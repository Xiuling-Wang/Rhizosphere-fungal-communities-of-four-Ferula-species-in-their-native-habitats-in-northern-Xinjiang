## Rebuild the script-generated publication figures from frozen processed data.
## This script does not recalculate the published statistical tables.
options(stringsAsFactors = FALSE)
set.seed(20260609)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(cowplot)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
  library(corrplot)
  library(RColorBrewer)
  library(ggplotify)
  library(ggrepel)
  library(ggtext)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
tb <- file.path(root, "data", "processed")
static_dir <- file.path(root, "data", "static")
fg <- file.path(root, "outputs", "figures", "components")
publication_dir <- file.path(root, "outputs", "figures", "publication")
dir.create(fg, recursive = TRUE, showWarnings = FALSE)
dir.create(publication_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  asv_counts = file.path(tb, "asv_table_fungi_rarefied_3819.tsv"),
  taxonomy = file.path(tb, "asv_taxonomy_fungi.tsv"),
  metadata = file.path(tb, "sample_metadata.tsv"),
  alpha_diversity = file.path(tb, "alpha_diversity_final.tsv"),
  biomarkers_traits = file.path(tb, "biomarkers_with_fungaltraits.tsv"),
  permanova = file.path(tb, "permanova_final.tsv"),
  hierarchical_partitioning = file.path(tb, "rdaccahp_final.tsv"),
  fungaltraits_groups = file.path(tb, "fungaltraits_guild_by_group.tsv"),
  core_genera = file.path(tb, "core_genera.tsv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing processed input files:\n- ", paste(missing_inputs, collapse = "\n- "), call. = FALSE)
}

rr <- read.delim(inputs[["asv_counts"]], check.names = FALSE)
rare <- as.matrix(rr[, -1]); rownames(rare) <- rr$sample
tax <- read.delim(inputs[["taxonomy"]], check.names = FALSE)
rownames(tax) <- tax$asv_id
meta <- read.delim(inputs[["metadata"]], check.names = FALSE)
meta <- meta[match(rownames(rare), meta$sample), ]
if (anyNA(meta$sample) || !all(colnames(rare) %in% tax$asv_id)) {
  stop("ASV counts, taxonomy, and metadata identifiers are not aligned.", call. = FALSE)
}
meta$site <- factor(meta$site, levels = c("HD", "XJ", "DS", "DG"))
meta$species_code <- meta$site
meta$depth_order <- factor(meta$depth_order, levels = c(1, 2, 3), labels = c("3 cm", "20 cm", "40 cm"))
go <- c("HD1","HD2","HD3","XJ1","XJ2","XJ3","DS1","DS2","DS3","DG1","DG2","DG3")
meta$group1 <- factor(meta$group1, levels = go)

site_cols  <- c(HD = "#4DBBD5", XJ = "#E64B35", DS = "#00A087", DG = "#925E9F")
## Sequential Blues: light (shallow) → dark (deep) — same hue, different intensity
depth_cols <- setNames(RColorBrewer::brewer.pal(9, "Blues")[c(3, 6, 9)],
                       c("3 cm", "20 cm", "40 cm"))

theme_paper <- function(base_size = 8, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid       = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 0.35),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      axis.text        = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", colour = "black"),
      legend.key       = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.title     = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0)
    )
}

save_plot <- function(name, plot, width, height) {
  ggsave(file.path(fg, paste0(name, ".pdf")), plot, width = width, height = height,
         units = "cm", bg = "white", device = cairo_pdf)
}

## ─────────────────────────────────────────────
## Figure 2: composition + ASV overlap
## ─────────────────────────────────────────────
long <- as.data.frame(rare) |>
  rownames_to_column("sample") |>
  pivot_longer(-sample, names_to = "asv_id", values_to = "abun") |>
  filter(abun > 0) |>
  left_join(tax[, c("asv_id", "Phylum", "Class", "Genus")], by = "asv_id") |>
  left_join(meta[, c("sample", "group1", "site")], by = "sample")

## Paired: 12 colours in light/dark pairs; drop #11 (pale yellow, invisible on white)
tax_pal <- RColorBrewer::brewer.pal(12, "Paired")[-11]

build_rank <- function(rank, req_class = FALSE, topn = NULL) {
  d <- long
  d$taxon <- d[[rank]]
  d$taxon[is.na(d$taxon) | d$taxon == ""] <- "Unassigned"
  if (req_class) d$taxon[d$taxon != "Unassigned" & (is.na(d$Class) | d$Class == "")] <- "Other"
  d <- d |>
    group_by(group1, taxon) |>
    summarise(ab = sum(abun), .groups = "drop")
  named <- d |>
    filter(!taxon %in% c("Unassigned", "Other")) |>
    group_by(taxon) |>
    summarise(total = sum(ab), .groups = "drop") |>
    arrange(desc(total)) |>
    pull(taxon)
  if (!is.null(topn) && length(named) > topn) {
    keep <- named[seq_len(topn)]
    d <- d |>
      mutate(taxon = ifelse(taxon %in% c(keep, "Unassigned"), taxon, "Other")) |>
      group_by(group1, taxon) |>
      summarise(ab = sum(ab), .groups = "drop")
    named <- keep
  }
  has_other <- "Other" %in% d$taxon
  d$taxon <- factor(d$taxon, levels = c("Unassigned", if (has_other) "Other", rev(named)))
  d$group1 <- factor(d$group1, levels = go)
  list(d = d, named = named, has_other = has_other)
}

composition_panel <- function(obj, title) {
  cols <- setNames(tax_pal[seq_along(obj$named)], obj$named)
  if (obj$has_other) cols["Other"]      <- "#C8C8C8"   # light grey
  cols["Unassigned"] <- "white"                         # white for unassigned
  ggplot(obj$d, aes(group1, ab, fill = taxon)) +
    geom_col(position = "fill", width = 0.74, colour = "grey70", linewidth = 0.18) +
    scale_fill_manual(values = cols,
                      breaks = c("Unassigned", if (obj$has_other) "Other", obj$named),
                      name = title) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x * 100)) +
    labs(x = NULL, y = "Relative abundance (%)") +
    theme_paper(7.4) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.8),
      axis.text.y  = element_text(size = 7),
      axis.title.y = element_text(size = 8, margin = margin(r = 3)),
      legend.position  = "right",
      legend.key.size  = unit(0.27, "cm"),
      legend.text      = element_text(size = 6.1),
      legend.title     = element_text(size = 7.2)
    )
}

bP <- build_rank("Phylum", topn = 10)
bC <- build_rank("Class", topn = 10)
bG <- build_rank("Genus", req_class = TRUE, topn = 10)

## Four-site ASV overlap (panel D) — UpSet-style panel drawn natively in ggplot2.
## Four-set Venn diagrams become visually crowded because all 15 regions must be
## labelled. The UpSet-style design keeps every count editable and makes the
## site-specific and shared ASV pools easier to compare.
pres <- rare > 0
upset_order <- c("HD", "XJ", "DS", "DG")
Pmat <- sapply(upset_order, function(s)
  colSums(pres[as.character(meta$site) == s, , drop = FALSE]) > 0)
.subsets <- unlist(lapply(1:4, function(k) combn(upset_order, k, simplify = FALSE)),
                   recursive = FALSE)
upset_df <- map_dfr(.subsets, function(S) {
  patt <- upset_order %in% S
  tibble(
    intersection = paste(S, collapse = "&"),
    degree = length(S),
    count = sum(apply(Pmat, 1, function(r) all(r == patt)))
  )
}) |>
  filter(count > 0) |>
  arrange(degree, desc(count), intersection) |>
  mutate(x = row_number(),
         intersection = factor(intersection, levels = intersection))

upset_matrix <- expand_grid(
  intersection = upset_df$intersection,
  set = factor(rev(upset_order), levels = rev(upset_order))
) |>
  left_join(upset_df[, c("intersection", "x")], by = "intersection") |>
  mutate(active = map2_lgl(as.character(intersection), as.character(set),
                           ~ .y %in% strsplit(.x, "&", fixed = TRUE)[[1]]),
         y = as.numeric(set))
upset_segments <- upset_matrix |>
  filter(active) |>
  group_by(intersection, x) |>
  summarise(ymin = min(as.numeric(set)), ymax = max(as.numeric(set)),
            degree = n(), .groups = "drop") |>
  filter(degree > 1)

matrix_y <- setNames(seq(-24, -84, length.out = length(upset_order)), upset_order)
upset_matrix$y_plot <- matrix_y[as.character(upset_matrix$set)]
upset_segments <- upset_matrix |>
  filter(active) |>
  group_by(intersection, x) |>
  summarise(ymin = min(y_plot), ymax = max(y_plot),
            degree = n(), .groups = "drop") |>
  filter(degree > 1)

pD <- ggplot() +
  geom_col(data = upset_df, aes(x = x, y = count, fill = factor(degree)),
           width = 0.70, colour = NA) +
  geom_text(data = upset_df, aes(x = x, y = count, label = count),
            vjust = -0.28, size = 2.45, family = "Helvetica",
            fontface = "bold", colour = "grey15") +
  geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.42) +
  geom_segment(data = upset_segments,
               aes(x = x, xend = x, y = ymin, yend = ymax),
               linewidth = 0.42, colour = "grey35", lineend = "round") +
  geom_point(data = upset_matrix,
             aes(x = x, y = y_plot, fill = ifelse(active, as.character(set), "inactive")),
             shape = 21, size = 2.35, stroke = 0.22, colour = "white") +
  geom_text(data = data.frame(set = upset_order, y = matrix_y[upset_order]),
            aes(x = 0.18, y = y, label = set, colour = set),
            inherit.aes = FALSE, hjust = 1, size = 2.75,
            family = "Helvetica", fontface = "bold") +
  scale_fill_manual(values = c(`1` = "#6E7781", `2` = "#8EA4BF",
                               `3` = "#B89AC8", `4` = "#E1B36A",
                               site_cols,
                               inactive = "#D8DCE2")) +
  scale_colour_manual(values = site_cols) +
  scale_x_continuous(breaks = upset_df$x, labels = rep("", nrow(upset_df)),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = c(0, 100, 200, 300),
                     labels = c("0", "100", "200", "300"),
                     expand = expansion(mult = c(0.02, 0.10))) +
  coord_cartesian(xlim = c(0.35, max(upset_df$x) + 0.60),
                  ylim = c(min(matrix_y) - 14, max(upset_df$count) * 1.16),
                  clip = "off") +
  labs(x = "Intersection", y = "ASVs") +
  theme_paper(7.2) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_text(size = 7.1, margin = margin(t = 2)),
    axis.title.y = element_text(size = 7.3, margin = margin(r = 2)),
    axis.text.y = element_text(size = 6.2),
    legend.position = "none",
    plot.margin = margin(4, 4, 2, 11)
  )

## Standalone component generated only temporarily for final assembly.
ggsave(file.path(fg, "Figure2D_overlap_upset_polished.pdf"), pD, width = 4.55, height = 4.45,
       device = cairo_pdf, bg = "white")

fig2 <- plot_grid(
  composition_panel(bP, "Phylum"),
  composition_panel(bC, "Class"),
  composition_panel(bG, "Genus"),
  pD,
  labels = c("A", "B", "C", "D"),
  label_fontface = "bold", label_size = 18,
  label_x = 0.006, label_y = 0.996,
  hjust = 0, vjust = 1,
  ncol = 2, align = "hv"
)
save_plot("Figure2_polished", fig2, 32, 20)
cat("Figure 2 done\n")

## ─────────────────────────────────────────────
## Figure 3AB: alpha diversity – QIIME2-style
## Facet by site × metric; depth on x-axis; depth_cols fill
## ─────────────────────────────────────────────
A_tbl <- read.delim(file.path(tb, "alpha_diversity_final.tsv"), check.names = FALSE)
A_tbl <- merge(A_tbl, meta[, c("sample", "site", "group1")], by = "sample")
A_tbl$site  <- factor(A_tbl$site,  levels = c("HD", "XJ", "DS", "DG"))
A_tbl$group1 <- factor(A_tbl$group1, levels = go)
A_tbl$depth_order <- factor(A_tbl$depth_order,
                             levels = c(1, 2, 3), labels = c("3 cm", "20 cm", "40 cm"))

letters_from_pairwise <- function(df, metric, group_var = "depth_order", alpha = 0.05) {
  groups <- levels(df[[group_var]])
  y <- df[[metric]]; g <- df[[group_var]]
  if (length(unique(g[!is.na(y)])) < 2) return(setNames(rep("a", length(groups)), groups))
  kw <- suppressWarnings(kruskal.test(y ~ g))
  ## Return empty strings when KW is not significant — no letters cluttering plot
  if (is.na(kw$p.value) || kw$p.value >= alpha) return(setNames(rep("", length(groups)), groups))
  pw <- suppressWarnings(
    pairwise.wilcox.test(y, g, p.adjust.method = "BH", exact = FALSE)$p.value)
  means <- tapply(y, g, median, na.rm = TRUE)
  ord   <- names(sort(means, decreasing = TRUE))
  let   <- setNames(rep("", length(groups)), groups)
  alphabet <- letters[seq_len(8)]
  for (grp in ord) {
    for (L in alphabet) {
      members <- names(let)[grepl(L, let, fixed = TRUE)]
      conflict <- FALSE
      if (length(members)) {
        for (m in members) {
          p <- if (grp %in% rownames(pw) && m %in% colnames(pw)) pw[grp, m]
               else if (m %in% rownames(pw) && grp %in% colnames(pw)) pw[m, grp]
               else NA
          if (!is.na(p) && p < alpha) conflict <- TRUE
        }
      }
      if (!conflict) { let[grp] <- paste0(let[grp], L); break }
    }
  }
  let[let == ""] <- "a"
  let
}

## Method note shown on the alpha-diversity depth panels (Fig.3A,B and FigS1):
## documents the test behind the compact-letter display so readers can tell
## that the absence of letters means "no significant within-site depth effect".
depth_method_caption <- "Kruskal-Wallis within each species/site: no significant depth effect."

alpha_long <- A_tbl |>
  select(sample, site, group1, depth_order, Shannon, Chao1) |>
  pivot_longer(c(Shannon, Chao1), names_to = "metric", values_to = "value")
alpha_long$metric <- factor(alpha_long$metric,
                             levels = c("Shannon", "Chao1"),
                             labels = c("Shannon diversity", "Chao1 richness"))

## significance letter positions – per (metric × site), label position slightly above max
letter_tbl <- alpha_long |>
  group_by(metric, site) |>
  group_modify(~{
    let <- letters_from_pairwise(.x, "value", group_var = "depth_order")
    data.frame(depth_order = names(let), label = unname(let))
  }) |>
  ungroup()
letter_tbl$depth_order <- factor(letter_tbl$depth_order, levels = c("3 cm", "20 cm", "40 cm"))

max_df <- alpha_long |>
  group_by(metric, site, depth_order) |>
  summarise(ymax = max(value, na.rm = TRUE), .groups = "drop")

## nudge y position: 4 % of each facet's range
nudge_df <- alpha_long |>
  group_by(metric) |>
  summarise(nudge = 0.04 * diff(range(value, na.rm = TRUE)), .groups = "drop")

lab_df <- left_join(max_df, letter_tbl, by = c("metric", "site", "depth_order")) |>
  left_join(nudge_df, by = "metric") |>
  mutate(y_label = ymax + nudge)

p3ab <- ggplot(alpha_long,
               aes(depth_order, value, fill = site, colour = site, alpha = depth_order)) +
  ## alternating background per site handled by facet border + white fill
  geom_boxplot(
    width = 0.55, outlier.shape = NA, linewidth = 0.4,
    fatten = 1.4
  ) +
  geom_jitter(
    width = 0.14, height = 0, size = 1.7, stroke = 0.3,
    shape = 21
  ) +
  geom_text(
    data = lab_df,
    aes(x = depth_order, y = y_label, label = label),
    inherit.aes = FALSE,
    size = 2.8, fontface = "bold", colour = "black", vjust = 0
  ) +
  facet_grid(metric ~ site, scales = "free_y", switch = "y") +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_colour_manual(values = site_cols, guide = "none") +
  scale_alpha_manual(values = c("3 cm" = 0.45, "20 cm" = 0.70, "40 cm" = 1.00),
                     name = "Depth") +
  labs(x = NULL, y = NULL) +
  theme_paper(8.5) +
  theme(
    panel.border    = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
    axis.line       = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 9, face = "bold"),
    strip.text.x    = element_text(size = 9.5, face = "bold"),
    panel.spacing   = unit(0.5, "lines"),
    axis.text.x     = element_text(size = 8),
    legend.position = "bottom",
    legend.title    = element_text(size = 7.2, face = "bold"),
    legend.text     = element_text(size = 6.8),
    legend.key.size = unit(0.40, "cm")
  ) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey35", colour = "grey35", size = 3)))
save_plot("Fig3AB_alpha_polished", p3ab, 22, 14)
cat("Fig3AB done\n")

## Single-metric helper (used in Figure 3 combined)
p_alpha_one <- function(metric_label, caption = NULL) {
  d      <- alpha_long |> filter(metric == metric_label)
  labs_one <- lab_df   |> filter(metric == metric_label)
  ggplot(d, aes(depth_order, value, fill = site, colour = site, alpha = depth_order)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.4, fatten = 1.4) +
    geom_jitter(width = 0.14, height = 0, size = 1.6, stroke = 0.28, shape = 21) +
    geom_text(data = labs_one,
              aes(x = depth_order, y = y_label, label = label),
              inherit.aes = FALSE, size = 2.5, fontface = "bold", colour = "black", vjust = 0) +
    facet_wrap(~site, nrow = 1) +
    scale_fill_manual(values = site_cols, guide = "none") +
    scale_colour_manual(values = site_cols, guide = "none") +
    scale_alpha_manual(values = c("3 cm" = 0.45, "20 cm" = 0.70, "40 cm" = 1.00),
                       name = "Depth") +
    labs(x = NULL, y = metric_label, caption = caption) +
    theme_paper(7.4) +
    theme(
      panel.border      = element_rect(colour = "grey40", fill = NA, linewidth = 0.45),
      axis.line         = element_blank(),
      strip.text        = element_text(size = 8.5, face = "bold"),
      panel.spacing     = unit(0.45, "lines"),
      axis.text.x       = element_text(size = 7),
      axis.title.y      = element_text(face = "bold", size = 8.2),
      legend.position   = "bottom",
      legend.title      = element_text(size = 6.8, face = "bold"),
      legend.text       = element_text(size = 6.5),
      legend.key.size   = unit(0.38, "cm"),
      plot.caption      = element_text(size = 8.4, hjust = 0, colour = "grey20",
                                       lineheight = 1.05, margin = margin(t = 5))
    ) +
    guides(alpha = guide_legend(override.aes = list(fill = "grey35", colour = "grey35", size = 3)))
}

## ─────────────────────────────────────────────
## Fig3C: alpha-soil-depth correlation (corrplot)
## ─────────────────────────────────────────────
alpha_all <- read.delim(file.path(tb, "alpha_diversity_final.tsv"), check.names = FALSE)
alpha_all <- merge(alpha_all, meta, by = "sample")
mets <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "Pielou")
soil_order <- c("SM","OM","TN","TP","TK","Nitrate_N","Ammonium_N","Olsen_P","AK","pH","EC","TDS")
cm_in <- cbind(alpha_all[, mets],
               meta[match(alpha_all$sample, meta$sample), soil_order],
               depth = as.numeric(meta$depth_order[match(alpha_all$sample, meta$sample)]))
pretty_names <- c(
  Observed = "Observed", Chao1 = "Chao1", ACE = "ACE", Shannon = "Shannon",
  Simpson = "Simpson", Pielou = "Pielou", SM = "SM", OM = "OM", TN = "TN",
  TP = "TP", TK = "TK", Nitrate_N = "NO3--N", Ammonium_N = "NH4+-N",
  Olsen_P = "Olsen P", AK = "AK", pH = "pH", EC = "EC", TDS = "TDS",
  depth = "Depth"
)
colnames(cm_in) <- unname(pretty_names[colnames(cm_in)])
cc <- cor(cm_in, method = "spearman", use = "pairwise")
pv <- matrix(NA_real_, ncol(cm_in), ncol(cm_in), dimnames = dimnames(cc))
for (i in seq_len(ncol(cm_in))) for (j in seq_len(ncol(cm_in)))
  pv[i, j] <- suppressWarnings(
    cor.test(cm_in[, i], cm_in[, j], method = "spearman", exact = FALSE)$p.value)
corr_cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200)

draw_fig3c <- function() {
  corrplot(cc, method = "color", type = "lower", order = "original",
           diag = FALSE, col = corr_cols, cl.ratio = 0.15, cl.align.text = "l",
           tl.col = "black", tl.cex = 0.82, tl.srt = 45, tl.pos = "ld",
           addgrid.col = "white", p.mat = pv,
           sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
           pch.col = "#222222", pch.cex = 1.18, mar = c(1.4, 0, 1, 0))
  mtext("Spearman's r", side = 1, line = 0.2, cex = 0.82, font = 2)
}
pdf(file.path(fg, "Fig3C_alpha_soil_corr_polished.pdf"), width = 7.4, height = 7.1, bg = "white")
draw_fig3c(); dev.off()

pC <- ggplotify::as.ggplot(~draw_fig3c()) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

fig3_combined <- plot_grid(
  plot_grid(
    p_alpha_one("Shannon diversity"),
    p_alpha_one("Chao1 richness", caption = depth_method_caption),
    labels = c("A", "B"), ncol = 1,
    label_fontface = "bold", label_size = 18,
    label_x = 0.006, label_y = 0.996,
    hjust = 0, vjust = 1, align = "v"
  ),
  pC,
  labels = c("", "C"), ncol = 2,
  rel_widths = c(1.0, 1.15),
  label_fontface = "bold", label_size = 18,
  label_x = c(0, 0.070), label_y = 0.990,
  hjust = 0, vjust = 1
)
save_plot("Figure3_polished", fig3_combined, 36, 18)
cat("Figure3 done\n")

## ─────────────────────────────────────────────
## Helper: aggregate by rank
## ─────────────────────────────────────────────
agg <- function(rank) {
  v    <- tax[colnames(rare), rank]
  keep <- !is.na(v) & v != ""
  cnt  <- t(rowsum(t(rare[, keep, drop = FALSE]), group = v[keep]))
  cnt / rowSums(rare)
}

## ─────────────────────────────────────────────
## Fig 4A: Class z-score heatmap (ComplexHeatmap)
##   – LEFT  annotation: Phylum colour bar  (mirrors old figure)
##   – RIGHT annotation: Mean RA (%) bar + Class name labels
##   – No row dendrogram; all classes with mean RA ≥ 0.01 %
## ─────────────────────────────────────────────
cl      <- agg("Class")
cls_site <- apply(cl, 2, function(x) tapply(x, meta$site, mean))

## Include all classes with mean RA ≥ 0.01 % (mirrors old figure's ~19 rows)
topc   <- names(sort(colMeans(cl), decreasing = TRUE))
topc   <- topc[colMeans(cl)[topc] * 100 >= 0.01]
n_cls  <- length(topc)

z      <- t(scale(cls_site[, topc]))   # rows = class, cols = site
z[is.na(z)] <- 0

## Mean RA across all samples for each class (%)
class_ra <- colMeans(cl)[topc] * 100
ra_max   <- ceiling(max(class_ra) / 5) * 5   # round up to nearest 5 (likely 20)

## ── Phylum for each class (majority of ASVs) ─────────────────────────────
class_phylum_raw <- sapply(topc, function(cls) {
  ph <- tax$Phylum[!is.na(tax$Class) & tax$Class == cls &
                   !is.na(tax$Phylum) & tax$Phylum != ""]
  if (!length(ph)) return("Other")
  names(which.max(table(ph)))
})
## Group minor phyla as "Other" to keep legend concise
main_ph <- c("Ascomycota","Basidiomycota","Glomeromycota",
             "Mortierellomycota","Mucoromycota","Chytridiomycota")
class_phylum_4a <- ifelse(class_phylum_raw %in% main_ph,
                          class_phylum_raw, "Other")
phylum_pal_full <- c(
  Ascomycota        = "#FC8D62",   # orange
  Basidiomycota     = "#8DA0CB",   # blue-purple
  Glomeromycota     = "#66C2A5",   # teal
  Mortierellomycota = "#E78AC3",   # pink
  Mucoromycota      = "#A6D854",   # yellow-green
  Chytridiomycota   = "#FFD92F",   # yellow
  Other             = "#B3B3B3"    # grey
)
## Keep only phyla actually present (Fig4A classes)
phylum_pal_4a <- phylum_pal_full[names(phylum_pal_full) %in% unique(class_phylum_4a)]

## ── Annotations ──────────────────────────────────────────────────────────
col_fun <- colorRamp2(seq(-2, 2, length.out = 101),
                      colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))

top_ann_4a <- HeatmapAnnotation(
  Site = factor(colnames(z), levels = c("HD","XJ","DS","DG")),
  col  = list(Site = site_cols),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size     = unit(0.35, "cm"),
  show_legend          = TRUE,
  annotation_legend_param = list(
    Site = list(title = "Site", title_gp = gpar(fontsize = 8, fontface = "bold"),
                labels_gp = gpar(fontsize = 7.5))
  )
)

## Left: thin Phylum colour bar (like old figure)
left_ann_4a <- rowAnnotation(
  Phylum = class_phylum_4a,
  col    = list(Phylum = phylum_pal_4a),
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontsize = 7.5, fontface = "bold"),
  simple_anno_size     = unit(0.35, "cm"),
  show_legend          = TRUE,
  annotation_legend_param = list(
    Phylum = list(title = "Phylum",
                  title_gp  = gpar(fontsize = 8, fontface = "bold"),
                  labels_gp = gpar(fontsize = 7.5),
                  ncol      = 1)
  )
)

## Right: RA bar (0 at heatmap side) + Class name labels
right_ann_4a <- rowAnnotation(
  "Mean RA (%)" = anno_barplot(
    class_ra,
    gp        = gpar(fill = "#4575b4", col = NA),
    bar_width = 0.72,
    axis_param = list(gp = gpar(fontsize = 6.5),
                      at     = c(0, ra_max/2, ra_max),
                      labels = c("0", as.character(ra_max/2), as.character(ra_max))),
    width = unit(2.2, "cm")
  ),
  "Class" = anno_text(
    names(class_ra),
    gp       = gpar(fontsize = 7.5),
    just     = "left",
    location = unit(1, "mm")
  ),
  annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
  annotation_name_side = c("Mean RA (%)" = "bottom"),
  show_annotation_name = c("Mean RA (%)" = TRUE, "Class" = FALSE),
  gap = unit(2, "mm")
)

## "Class" header text drawn manually above the Class-name column
## (annotation_name_side = "top" is unreliable for row annotations)
add_class_header_4a <- function() {
  decorate_annotation("Class", {
    grid.text("Class", x = unit(0, "npc"), y = unit(1, "npc") + unit(2, "mm"),
              just = c("left", "bottom"), gp = gpar(fontsize = 8, fontface = "bold"))
  })
}

ht4a <- Heatmap(
  z,
  name              = "z-score",
  col               = col_fun,
  top_annotation    = top_ann_4a,
  left_annotation   = left_ann_4a,
  right_annotation  = right_ann_4a,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_dend     = FALSE,
  show_column_dend  = TRUE,
  column_dend_height = unit(1, "cm"),
  show_row_names    = FALSE,
  column_names_gp   = gpar(fontsize = 8.5),
  column_names_rot  = 0,
  rect_gp           = gpar(col = "white", lwd = 0.5),
  border            = TRUE,
  heatmap_legend_param = list(
    title     = "z-score",
    title_gp  = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.5),
    at        = c(-2, -1, 0, 1, 2),
    direction = "vertical"
  ),
  width  = unit(11, "cm"),
  height = unit(n_cls * 0.72, "cm")   # dynamic: larger panel A heatmap body
)

## NOTE: `topc` currently keeps ALL classes with mean RA >= 0.01 % (n_cls classes,
## see cat() message below) -- not a fixed top-15/20 cutoff. This mirrors the old
## figure's ~19-row layout. Adjust the `>= 0.01` threshold above for a stricter cap.
cat(sprintf("Fig4A: showing %d classes with mean RA >= 0.01%%\n", n_cls))

fig4a_w <- 11.8                         # inches: tightened to content so panel A's enlarged body has no gap before B
fig4a_h <- (n_cls * 0.72 / 2.54) + 4.2   # inches: heatmap body + margins + column dendrogram
pdf(file.path(fg, "Fig4A_class_heatmap_polished.pdf"),
    width = fig4a_w, height = fig4a_h, bg = "white")
draw(ht4a, padding = unit(c(12, 3, 4, 3), "mm"),
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     merge_legend           = TRUE)
add_class_header_4a()
dev.off()
cat("Fig4A done\n")

## ─────────────────────────────────────────────
## Fig 4B: biomarker dotplot – genus label coloured by Class
## Recommendation: Class (5 classes, manageable; more informative than Phylum)
## ─────────────────────────────────────────────
bm_ft <- read.delim(file.path(tb, "biomarkers_with_fungaltraits.tsv"), check.names = FALSE)
bm_ft$enriched <- factor(bm_ft$enriched, levels = c("HD", "XJ", "DS", "DG"))

## Map genus → Class via taxonomy table
gen_tax <- unique(tax[, c("Genus", "Class", "Phylum")])
bm_ft <- merge(bm_ft, gen_tax, by.x = "genus", by.y = "Genus", all.x = TRUE)

## Sort: site first, then IndVal ascending within site
bm_ft <- bm_ft[order(bm_ft$enriched, bm_ft$IndVal), ]
bm_ft$genus_label <- factor(bm_ft$genus, levels = bm_ft$genus)

## Assign Class colours (5 classes: use a qualitative palette)
class_levels <- sort(unique(bm_ft$Class[!is.na(bm_ft$Class)]))
class_cols <- setNames(
  c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")[seq_along(class_levels)],
  class_levels
)

## Build italic genus + class-coloured label as markdown; append "(g)" = genus rank
bm_ft$md_label <- paste0(
  "<span style='color:", class_cols[bm_ft$Class], "'>*", bm_ft$genus, "*</span> (g)"
)
bm_ft$md_label <- factor(bm_ft$md_label, levels = bm_ft$md_label)

## Phylum group for left-hand colour strip (mirrors Fig4A left Phylum annotation)
bm_ft$phylum_grp <- ifelse(!is.na(bm_ft$Phylum) & bm_ft$Phylum %in% main_ph,
                            bm_ft$Phylum, "Other")
phylum_pal_4b <- phylum_pal_full[names(phylum_pal_full) %in% unique(bm_ft$phylum_grp)]

p4b <- ggplot(bm_ft, aes(IndVal, md_label, colour = enriched, size = mean_RA)) +
  geom_point(alpha = 0.90) +
  scale_colour_manual(values = site_cols, name = "Enriched site", drop = TRUE) +
  scale_size_continuous(name = "Mean RA (%)", range = c(1.8, 9),
                        breaks = c(0.5, 1, 2, 4)) +
  scale_x_continuous(limits = c(0.47, 0.86), breaks = seq(0.5, 0.8, 0.1),
                     expand = expansion(add = c(0.01, 0.06))) +
  labs(x = "IndVal (fidelity)", y = NULL) +
  theme_paper(8.5) +
  theme(
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.28),
    axis.line.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    legend.position   = "right",
    legend.key.size   = unit(0.38, "cm"),
    legend.text       = element_text(size = 7),
    legend.title      = element_text(size = 7.5, face = "bold"),
    axis.text.y       = element_markdown(size = 7.8)
  )

## ── Left-hand Phylum colour strip (one tile per genus row, aligned with p4b) ──
## height = 1 + no border -> adjacent rows of the same Phylum form one
## continuous coloured block, mirroring Fig4A's left Phylum bar.
phylum_strip_4b <- ggplot(bm_ft, aes(x = 1, y = md_label, fill = phylum_grp)) +
  geom_tile(width = 1, height = 1, colour = NA) +
  scale_fill_manual(values = phylum_pal_4b, drop = FALSE) +
  scale_y_discrete(limits = levels(bm_ft$md_label)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Phylum") +
  theme_void() +
  theme(
    plot.title  = element_text(size = 7, face = "bold", angle = 90,
                                hjust = 0.5, vjust = 0.5, margin = margin(r = 2)),
    legend.position = "none",
    plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 1)
  )

## ── Compact legends, stacked vertically ────────────────────────────────────
## 1) main legend (Enriched site colour + Mean RA size) from p4b itself
main_legend <- get_legend(p4b)

## 2) Phylum colour legend (single column, matches the left strip)
phylum_legend_plot <- ggplot(data.frame(Phylum = factor(names(phylum_pal_4b), levels = names(phylum_pal_4b))),
                              aes(x = 1, y = Phylum, fill = Phylum)) +
  geom_tile() +
  scale_fill_manual(values = phylum_pal_4b) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.34, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.box.margin = margin(0, 0, 0, 0)
  )
phylum_legend <- get_legend(phylum_legend_plot)

## 3) Class colour legend (single column, compact swatches matching label colours)
class_legend_plot <- ggplot(data.frame(Class = factor(class_levels, levels = class_levels)),
                             aes(x = 1, y = Class, colour = Class)) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = class_cols) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.34, "cm"),
    legend.text      = element_text(size = 7, face = "italic"),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.box.margin = margin(0, 0, 0, 0)
  )
class_legend <- get_legend(class_legend_plot)

legends_stack_4b <- plot_grid(main_legend, phylum_legend, class_legend,
                               ncol = 1, rel_heights = c(1, 0.45, 0.55),
                               align = "v", axis = "l")

p4b_main <- p4b + theme(legend.position = "none")

p4b_combined <- plot_grid(
  phylum_strip_4b, p4b_main, legends_stack_4b,
  ncol = 3, rel_widths = c(0.035, 1, 0.30),
  align = "h", axis = "tb"
)
save_plot("Fig4B_biomarkers_polished", p4b_combined, 19, 13)
cat("Fig4B done\n")

## ─────────────────────────────────────────────
## Fig 4 (combined): Fig4A (Class heatmap) + Fig4B (Genus biomarkers) side-by-side,
## drawn as panels A / B of a single Figure 4. Both 4A and 4B are still saved
## individually above for flexibility.
## Drawn directly into nested viewports (not grid.grabExpr) so ComplexHeatmap's
## legend layout is computed against the correct panel size and nothing clips.
## ─────────────────────────────────────────────
fig4b_w_in       <- 16.5 / 2.54
fig4b_h_in       <- 13 / 2.54
fig4_combined_w_in <- fig4a_w + fig4b_w_in
fig4_combined_h_in <- fig4a_h
panelA_frac <- fig4a_w / fig4_combined_w_in

draw_fig4_combined <- function() {
  pushViewport(viewport(x = 0, y = 0, width = panelA_frac, height = 1, just = c("left", "bottom")))
  draw(ht4a, padding = unit(c(12, 3, 4, 3), "mm"),
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legend           = TRUE,
       newpage                = FALSE)
  add_class_header_4a()
  upViewport()

  ## Panel B drawn at its native size (19 x 13 cm), vertically centred in its
  ## column -- avoids stretching B to A's (taller) height, which previously
  ## made B look oversized relative to A.
  pushViewport(viewport(
    x      = unit(fig4a_w + fig4b_w_in / 2, "in"),
    y      = unit(fig4_combined_h_in / 2, "in"),
    width  = unit(fig4b_w_in, "in"),
    height = unit(fig4b_h_in, "in"),
    just   = "center"
  ))
  print(p4b_combined, newpage = FALSE)
  upViewport()

  ## Panel labels A / B: placed close to each panel but away from rank labels.
  grid.text("A", x = unit(0.105, "npc"), y = unit(0.845, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text("B", x = unit(panelA_frac + 0.052, "npc"),
            y = unit((fig4_combined_h_in / 2 + fig4b_h_in / 2) / fig4_combined_h_in, "npc") - unit(4, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 20, fontface = "bold"))
}

pdf(file.path(fg, "Fig4_combined_polished.pdf"),
    width = fig4_combined_w_in, height = fig4_combined_h_in, bg = "white")
draw_fig4_combined()
dev.off()

cat("Fig4_combined done\n")

## ─────────────────────────────────────────────
## Fig 5A: dbRDA biplot (unchanged aesthetics)
## ─────────────────────────────────────────────
dbrda_sel <- c("SM", "OM", "TK", "Nitrate_N", "Ammonium_N", "Olsen_P", "AK", "EC")
all_env   <- c("SM","OM","TN","TP","TK","Nitrate_N","Ammonium_N","Olsen_P","AK","pH","EC","TDS")
envz_all  <- as.data.frame(scale(meta[, all_env]))
envz_sel  <- envz_all[, dbrda_sel]
bray5a    <- vegdist(rare, "bray")
db5a      <- capscale(as.formula(paste("bray5a ~", paste(dbrda_sel, collapse = "+"))),
                      data = envz_sel)
eig5a     <- eigenvals(db5a)
cap_eig5a <- eig5a[grepl("^CAP", names(eig5a))]
tot5a     <- sum(eig5a)
cap1_lab  <- paste0("CAP1 (", round(100 * cap_eig5a[1] / tot5a, 1), "%)")
cap2_lab  <- paste0("CAP2 (", round(100 * cap_eig5a[2] / tot5a, 1), "%)")
st5a <- as.data.frame(scores(db5a, display = "sites", choices = 1:2))
st5a$sample <- rownames(st5a)
st5a <- merge(meta[, c("sample", "site", "depth_order")], st5a, by = "sample")
st5a$site        <- factor(st5a$site, levels = c("HD", "XJ", "DS", "DG"))
st5a$depth_label <- factor(as.character(st5a$depth_order), levels = c("3 cm", "20 cm", "40 cm"))
bp5a <- as.data.frame(scores(db5a, display = "bp", choices = 1:2))
bp5a$var <- rownames(bp5a)
label_map <- c(Nitrate_N = "NO[3]^'-'*'-N'", Ammonium_N = "NH[4]^'+'*'-N'",
               Olsen_P = "Olsen~P", SM = "SM", OM = "OM", TK = "TK", AK = "AK", EC = "EC")
bp5a$var_label <- ifelse(bp5a$var %in% names(label_map), label_map[bp5a$var], bp5a$var)
sc5a <- 0.80 * max(abs(st5a[, c("CAP1", "CAP2")])) / max(abs(bp5a[, 1:2]))

p5a <- ggplot(st5a, aes(CAP1, CAP2, colour = site, alpha = depth_label)) +
  geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.28, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.28, linetype = "dashed") +
  geom_point(size = 2.6, stroke = 0.45, shape = 16) +
  geom_segment(data = bp5a, aes(x = 0, y = 0, xend = CAP1 * sc5a, yend = CAP2 * sc5a),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
               colour = "grey22", linewidth = 0.55) +
  geom_text_repel(data = bp5a, aes(x = CAP1 * sc5a, y = CAP2 * sc5a, label = var_label),
                  inherit.aes = FALSE, size = 2.65, colour = "grey10",
                  segment.size = 0.28, box.padding = 0.32,
                  force = 0.9, seed = 42, max.overlaps = 20,
                  parse = TRUE) +
  scale_colour_manual(values = site_cols, name = "Species/site") +
  scale_alpha_manual(values = c("3 cm" = 1, "20 cm" = 0.62, "40 cm" = 0.32), name = "Depth") +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1)),
    alpha  = guide_legend(override.aes = list(colour = "grey30"))
  ) +
  labs(x = cap1_lab, y = cap2_lab) +
  theme_paper(8.5) +
  theme(
    panel.border     = element_rect(colour = "grey55", fill = NA, linewidth = 0.4),
    axis.line        = element_blank(),
    legend.position  = "right",
    legend.key.size  = unit(0.38, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7.5, face = "bold")
  )
## Add PERMANOVA (adonis2) annotation ─────────────────────────────────────
perm5a <- read.delim(file.path(tb, "permanova_final.tsv"),
                     check.names = FALSE, row.names = 1)
perm_site  <- perm5a["site", ]
perm_depth <- perm5a["depth_order", ]
perm_label <- sprintf(
  "PERMANOVA (adonis2)\nSpecies/site (plant-level): R^2 = %.3f, p = %.3f%s\nDepth (within-plant): R^2 = %.3f, p = %.3f%s",
  perm_site[["R2"]], perm_site[["Pr(>F)"]],
  ifelse(perm_site[["Pr(>F)"]] <= 0.001, "***",
         ifelse(perm_site[["Pr(>F)"]] <= 0.01, "**",
                ifelse(perm_site[["Pr(>F)"]] <= 0.05, "*", " ns"))),
  perm_depth[["R2"]], perm_depth[["Pr(>F)"]],
  ifelse(perm_depth[["Pr(>F)"]] <= 0.001, "***",
         ifelse(perm_depth[["Pr(>F)"]] <= 0.01, "**",
                ifelse(perm_depth[["Pr(>F)"]] <= 0.05, "*", " ns")))
)
p5a <- p5a +
  annotate("text",
           x = Inf, y = -Inf,
           label = perm_label,
           hjust = 1.04, vjust = -0.18,
           size  = 2.4, colour = "grey20",
           family = "Helvetica", lineheight = 1.35,
           fontface = "plain")

save_plot("Fig5A_dbRDA_polished", p5a, 15, 12)
cat("Fig5A done\n")

## ─────────────────────────────────────────────
## Fig 5B: rdacca.hp bar chart
## ─────────────────────────────────────────────
hp_df <- read.delim(file.path(tb, "rdaccahp_final.tsv"), check.names = FALSE)
hp_label_map <- c(Nitrate_N = "NO[3]^'-'*'-N'", Ammonium_N = "NH[4]^'+'*'-N'",
                  Olsen_P = "Olsen~P", SM = "SM", OM = "OM",
                  TK = "TK", AK = "AK", EC = "EC")
hp_df$var_label <- ifelse(hp_df$var %in% names(hp_label_map),
                          hp_label_map[hp_df$var], hp_df$var)
hp_df$var_label <- factor(hp_df$var_label,
                           levels = hp_df$var_label[order(hp_df$ind)])

p5b <- ggplot(hp_df, aes(ind * 100, var_label)) +
  geom_col(fill = "#4575b4", width = 0.66) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, NA)) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(x = expression("Individual effect (% adj " * R^2 * ")"), y = NULL) +
  theme_paper(8.5) +
  theme(
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.28),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_text(size = 8)
  )
save_plot("Fig5B_varpart_polished", p5b, 11, 8)
cat("Fig5B done\n")

## ─────────────────────────────────────────────
## Fig 5C: Genus–soil Spearman heatmap (ComplexHeatmap)
##   – large significance symbols filling cell
##   – "Genus" row title; colour bar labelled "Spearman's r"
##   – column bar annotation for soil variable type
## ─────────────────────────────────────────────
g    <- agg("Genus")
topg <- names(sort(colMeans(g), decreasing = TRUE))[seq_len(min(30, ncol(g)))]
G    <- g[, topg, drop = FALSE]
cat("Fig5C: showing", length(topg), "of", ncol(g), "named genera (top by mean RA)\n")

env_vars <- c("SM","OM","TN","TP","TK","Nitrate_N","Ammonium_N","Olsen_P","AK","pH","EC","TDS")
env_label_expr <- expression(SM, OM, TN, TP, TK, NO[3]^"-"*"-N", NH[4]^"+"*"-N",
                             Olsen~P, AK, pH, EC, TDS)

rho <- matrix(NA_real_, length(topg), length(env_vars),
              dimnames = list(topg, env_vars))
pv5c <- rho
for (i in topg) for (j in env_vars) {
  ct <- suppressWarnings(
    cor.test(G[, i], meta[[j]], method = "spearman", exact = FALSE))
  rho[i, j]  <- ct$estimate
  pv5c[i, j] <- ct$p.value
}
colnames(rho) <- env_vars; colnames(pv5c) <- env_vars

## Cell size: aim for symbols to fill the cell
## Use a large fontsize so *** fills the box
cell_fn <- function(j, i, x, y, width, height, fill) {
  p <- pv5c[i, j]
  sym <- if (!is.na(p) && p < 0.001) "***"
         else if (!is.na(p) && p < 0.01)  "**"
         else if (!is.na(p) && p < 0.05)  "*"
         else ""
  if (sym != "") {
    grid.text(sym, x, y,
              gp = gpar(fontsize = if (sym == "***") 14.2 else if (sym == "**") 14.5 else 15,
                        col = "black", fontface = "bold"))
  }
}

col5c <- colorRamp2(seq(-1, 1, length.out = 101),
                    colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))

ht5c <- Heatmap(
  rho,
  name             = "Spearman's r",
  col              = col5c,
  cell_fun         = cell_fn,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_dend    = TRUE,
  show_column_dend = TRUE,
  row_title        = NULL,
  row_names_side   = "right",
  row_names_gp     = gpar(fontsize = 8.5, fontface = "italic"),
  column_labels    = env_label_expr,
  column_names_gp  = gpar(fontsize = 9),
  column_names_rot = 45,
  rect_gp          = gpar(col = "white", lwd = 0.5),
  border           = TRUE,
  heatmap_legend_param = list(
    title     = "Spearman's r",
    title_gp  = gpar(fontsize = 8.5, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    at        = c(-1, -0.5, 0, 0.5, 1),
    direction = "vertical",
    legend_height = unit(3.5, "cm")
  ),
  width  = unit(9, "cm"),
  height = unit(14, "cm")
)

fig5c_w_in <- 7.3
fig5c_h_in <- 7.7
fig5c_pad  <- unit(c(2, 2, 2, 2), "mm")

pdf(file.path(fg, "Fig5C_genus_soil_polished.pdf"), width = fig5c_w_in, height = fig5c_h_in, bg = "white")
draw(ht5c, padding = fig5c_pad)
grid.text("Genera", x = unit(1, "npc") - unit(34, "mm"),
          y = unit(1, "npc") - unit(11, "mm"),
          just = c("center", "bottom"),
          gp = gpar(fontsize = 9, fontface = "bold"))
dev.off()
cat("Fig5C done\n")

## ─────────────────────────────────────────────
## Fig 5 combined: A (dbRDA) + B (varpart) stacked on the left,
##                 C (genus–soil heatmap) on the right
## ─────────────────────────────────────────────
fig5a_w_in <- 15 / 2.54
fig5a_h_in <- 12 / 2.54
fig5b_w_in <- 15 / 2.54   # widen B to match A's width in the combined layout (avoids dead space)
fig5b_h_in <- 8  / 2.54

gap_ab  <- 0.30   # vertical gap between A and B
gap_lc  <- 0.30   # horizontal gap between left column and C

left_col_w <- fig5a_w_in
left_col_h <- fig5a_h_in + gap_ab + fig5b_h_in

fig5_w_in <- left_col_w + gap_lc + fig5c_w_in
fig5_h_in <- max(left_col_h, fig5c_h_in)

pdf(file.path(fg, "Fig5_combined_polished.pdf"), width = fig5_w_in, height = fig5_h_in, bg = "white")

## Panel A: top-left
pushViewport(viewport(x = unit(0, "in"), y = unit(fig5_h_in, "in"),
                       width = unit(fig5a_w_in, "in"), height = unit(fig5a_h_in, "in"),
                       just = c("left", "top")))
print(p5a, newpage = FALSE)
upViewport()

## Panel B: bottom-left
pushViewport(viewport(x = unit(0, "in"), y = unit(fig5_h_in - fig5a_h_in - gap_ab, "in"),
                       width = unit(fig5b_w_in, "in"), height = unit(fig5b_h_in, "in"),
                       just = c("left", "top")))
print(p5b, newpage = FALSE)
upViewport()

## Panel C: right column, top-aligned with A
pushViewport(viewport(x = unit(left_col_w + gap_lc, "in"), y = unit(fig5_h_in, "in"),
                       width = unit(fig5c_w_in, "in"), height = unit(fig5c_h_in, "in"),
                       just = c("left", "top")))
draw(ht5c, padding = fig5c_pad, newpage = FALSE)
grid.text("Genera", x = unit(1, "npc") - unit(34, "mm"),
          y = unit(1, "npc") - unit(11, "mm"),
          just = c("center", "bottom"),
          gp = gpar(fontsize = 9, fontface = "bold"))
upViewport()

## Panel labels
grid.text("A", x = unit(2, "mm"), y = unit(fig5_h_in, "in") - unit(2, "mm"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("B", x = unit(2, "mm"), y = unit(fig5_h_in - fig5a_h_in - gap_ab, "in") + unit(2, "mm"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("C", x = unit(left_col_w + gap_lc, "in") + unit(2, "mm"), y = unit(fig5_h_in, "in") - unit(2, "mm"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))

dev.off()
cat("Fig5_combined done\n")

## ─────────────────────────────────────────────
## FigS1: FungalTraits guild composition
##   Panel A – stacked bar per group1 (12 groups, site-labelled)
##   Panel B – key functional guilds compared across sites (dot plot)
## ─────────────────────────────────────────────
ft_gl <- read.delim(file.path(tb, "fungaltraits_guild_by_group.tsv"), check.names = FALSE)
ft_gl <- merge(ft_gl, unique(meta[, c("group1", "site")]), by = "group1")
ft_gl$group1 <- factor(ft_gl$group1, levels = go)
ft_gl$site   <- factor(ft_gl$site,   levels = c("HD", "XJ", "DS", "DG"))

## Top 9 guilds by total abundance (excl. Unassigned)
top9 <- ft_gl |>
  filter(life != "Unassigned/no-trait") |>
  group_by(life) |>
  summarise(s = sum(ab), .groups = "drop") |>
  arrange(desc(s)) |>
  slice_head(n = 9) |>
  pull(life)

ft_gl <- ft_gl |>
  mutate(lab = case_when(
    life == "Unassigned/no-trait" ~ "Unassigned",
    life %in% top9               ~ life,
    TRUE                         ~ "Other"
  ))

## ggplot2 stacking rule: FIRST factor level → TOP (y=1); LAST → BOTTOM (y=0)
## Mirrors Fig 2 logic: levels = c("Unassigned", "Other", rev(named))
##   → Unassigned at top (white), most abundant named guild at bottom
## ggplot2: FIRST level → TOP (y=1); LAST level → BOTTOM (y=0)
## → Unassigned first (white, top); soil_saprotroph last (most abundant, bottom)
stack_levels <- c(
  "Unassigned",             # first  → TOP  (white)
  "Other",                  #          2nd from top (grey)
  "dung_saprotroph",        # least abundant named
  "animal_parasite",
  "arbuscular_mycorrhizal",
  "ectomycorrhizal",
  "sooty_mold",
  "plant_pathogen",
  "litter_saprotroph",
  "wood_saprotroph",
  "unspecified_saprotroph",
  "soil_saprotroph"         # last   → BOTTOM (y=0, most abundant)
)
stack_levels <- stack_levels[stack_levels %in% unique(ft_gl$lab)]

guild_labels_s <- c(
  Unassigned             = "Unassigned",
  Other                  = "Other",
  dung_saprotroph        = "Dung saprotroph",
  animal_parasite        = "Animal parasite",
  arbuscular_mycorrhizal = "Arbuscular mycorrhizal",
  ectomycorrhizal        = "Ectomycorrhizal",
  sooty_mold             = "Sooty mold",
  plant_pathogen         = "Plant pathogen",
  litter_saprotroph      = "Litter saprotroph",
  wood_saprotroph        = "Wood saprotroph",
  unspecified_saprotroph = "Unspecified saprotroph",
  soil_saprotroph        = "Soil saprotroph"
)
guild_pal_s <- c(
  Unassigned             = "white",
  Other                  = "#c7c7c7",
  dung_saprotroph        = "#8c564b",
  animal_parasite        = "#a65628",
  arbuscular_mycorrhizal = "#238b45",
  ectomycorrhizal        = "#74c476",
  sooty_mold             = "#984ea3",
  plant_pathogen         = "#de2d26",
  litter_saprotroph      = "#fdbe85",
  wood_saprotroph        = "#fd8d3c",
  unspecified_saprotroph = "#9ecae1",
  soil_saprotroph        = "#4292c6"
)

ft_agg <- ft_gl |>
  group_by(group1, site, lab) |>
  summarise(rel = sum(rel), .groups = "drop") |>
  mutate(lab = factor(lab, levels = stack_levels))

## Panel A: stacked bar with site background shading + labels
pSA <- ggplot(ft_agg, aes(group1, rel, fill = lab)) +
  annotate("rect",
    xmin = c(0.5, 6.5), xmax = c(3.5, 9.5),
    ymin = -Inf, ymax = Inf,
    fill = "grey92", alpha = 0.55) +
  geom_col(position = "fill", width = 0.76, colour = "grey65", linewidth = 0.15) +
  annotate("text",
    x = c(2, 5, 8, 11), y = 1.045,
    label = c("HD", "XJ", "DS", "DG"),
    size = 2.9, fontface = "bold", vjust = 0,
    colour = unname(site_cols[c("HD", "XJ", "DS", "DG")])) +
  scale_fill_manual(
    values = guild_pal_s,
    labels = guild_labels_s[stack_levels],
    breaks = stack_levels,
    name   = "Primary lifestyle"
  ) +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) paste0(x * 100),
                     limits = c(0, 1.08)) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_paper(7.8) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
    legend.key.size = unit(0.28, "cm"),
    legend.text   = element_text(size = 6.5),
    legend.title  = element_text(size = 7.5, face = "bold"),
    plot.margin   = margin(t = 14, r = 5, b = 5, l = 5)
  )

## Panel B: key guild dot plot per site (ordered by mean RA ascending → AMF at top)
focus_life <- c(
  "arbuscular_mycorrhizal", "ectomycorrhizal", "plant_pathogen",
  "soil_saprotroph", "litter_saprotroph", "wood_saprotroph",
  "unspecified_saprotroph"
)

site_gl <- ft_gl |>
  group_by(site, life) |>
  summarise(ab = sum(ab), .groups = "drop") |>
  group_by(site) |>
  mutate(rel_pct = ab / sum(ab) * 100) |>
  ungroup() |>
  filter(life %in% focus_life)

## Dot-plot order: sort by mean RA ascending (least abundant at bottom)
dot_order <- site_gl |>
  group_by(life) |>
  summarise(m = mean(rel_pct), .groups = "drop") |>
  arrange(m) |>
  pull(life)

site_gl <- site_gl |>
  mutate(guild_disp = factor(
    unname(guild_labels_s[life]),
    levels = unname(guild_labels_s[dot_order])
  ))

pSB <- ggplot(site_gl, aes(rel_pct, guild_disp, colour = site)) +
  geom_line(aes(group = guild_disp), colour = "grey78", linewidth = 0.5) +
  geom_point(size = 3.5, alpha = 0.92) +
  scale_colour_manual(values = site_cols, name = "Site") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "Relative abundance (%)", y = NULL) +
  theme_paper(8.5) +
  theme(
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.28),
    axis.line.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.text.y   = element_text(size = 8),
    legend.key.size = unit(0.38, "cm"),
    legend.text   = element_text(size = 7.5),
    legend.title  = element_text(size = 8, face = "bold")
  )

figS1 <- plot_grid(pSA, pSB,
  labels = c("A", "B"), label_fontface = "bold", label_size = 18,
  label_x = 0.006, label_y = 0.996, hjust = 0, vjust = 1,
  ncol = 2, rel_widths = c(1.55, 1))
save_plot("FigS1_fungaltraits_polished", figS1, 34, 12)
cat("FigS1 done\n")

## ─────────────────────────────────────────────
## FigS2: Core fungal genera (34 genera shared by all 4 Ferula species)
##   – sorted by mean RA% descending (most abundant at bottom),
##     Unassigned at very top; full legend names; white for Unassigned
## ─────────────────────────────────────────────
cg <- read.delim(file.path(tb, "core_genera.tsv"), check.names = FALSE)
cg$primary_lifestyle[is.na(cg$primary_lifestyle) | cg$primary_lifestyle == ""] <- "Unassigned"

## Sort: primary = n_core_ASVs descending (most → BOTTOM of chart);
##       secondary = lifestyle group (same colour adjacent), Unassigned LAST (→ TOP);
##       tertiary = genus name alphabetically within each lifestyle group.
## Lifestyle rank order: saprotrophs → AMF → pathogen/sooty_mold → Unassigned
ls_sort_order <- c(
  "soil_saprotroph", "unspecified_saprotroph",
  "wood_saprotroph", "litter_saprotroph", "dung_saprotroph",
  "arbuscular_mycorrhizal",
  "plant_pathogen", "sooty_mold",
  "Unassigned"          # placed last → sorts to TOP of chart
)

cg <- cg |>
  mutate(ls_rank = match(primary_lifestyle, ls_sort_order)) |>
  arrange(desc(n_core_ASVs), ls_rank, genus) |>
  select(-ls_rank)
## In coord_flip bar charts: first factor level → BOTTOM (highest ASV count)
cg$genus_f <- factor(cg$genus, levels = cg$genus)

## Only lifestyles actually present in this dataset
present_ls <- unique(cg$primary_lifestyle)

## Legend order: ecologically meaningful (AMF first, then saprotrophs by
## descending community importance, then pathogen/sooty_mold, then Unassigned)
ls_order <- c(
  "arbuscular_mycorrhizal",
  "unspecified_saprotroph", "soil_saprotroph",
  "wood_saprotroph",        "litter_saprotroph",
  "dung_saprotroph",
  "plant_pathogen",         "sooty_mold",
  "Unassigned"
)
ls_order <- ls_order[ls_order %in% present_ls]   # keep only those present

## Full, unabbreviated names (no "Arb.", no "Unspec.")
core_labels_full <- c(
  arbuscular_mycorrhizal = "Arbuscular mycorrhizal",
  unspecified_saprotroph = "Unspecified saprotroph",
  soil_saprotroph        = "Soil saprotroph",
  wood_saprotroph        = "Wood saprotroph",
  litter_saprotroph      = "Litter saprotroph",
  dung_saprotroph        = "Dung saprotroph",
  plant_pathogen         = "Plant pathogen",
  sooty_mold             = "Sooty mold",
  Unassigned             = "Unassigned"
)[ls_order]

## Colours: match main composition bars; Unassigned = white (with border)
core_pal2 <- c(
  arbuscular_mycorrhizal = "#238b45",
  unspecified_saprotroph = "#9ecae1",
  soil_saprotroph        = "#4292c6",
  wood_saprotroph        = "#fd8d3c",
  litter_saprotroph      = "#fdbe85",
  dung_saprotroph        = "#8c564b",
  plant_pathogen         = "#de2d26",
  sooty_mold             = "#984ea3",
  Unassigned             = "white"
)[ls_order]

pS2 <- ggplot(cg, aes(n_core_ASVs, genus_f, fill = primary_lifestyle)) +
  geom_col(width = 0.72, colour = "grey65", linewidth = 0.18) +
  geom_text(aes(label = n_core_ASVs), hjust = -0.3, size = 2.6, colour = "grey30") +
  scale_fill_manual(
    values = core_pal2,
    labels = core_labels_full,
    breaks = ls_order,
    name   = "Primary lifestyle"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20)), breaks = 1:6) +
  labs(
    x     = "No. of core ASVs",
    y     = NULL,
    title = "Core fungal genera shared across all four Ferula species (n = 34)"
  ) +
  theme_paper(8) +
  theme(
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_text(face = "italic", size = 7.8),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.28),
    legend.key.size = unit(0.32, "cm"),
    legend.key      = element_rect(fill = NA, colour = "grey70", linewidth = 0.3),
    legend.text     = element_text(size = 7.2),
    legend.title    = element_text(size = 7.8, face = "bold"),
    plot.title      = element_text(size = 8.5, face = "plain",
                                   margin = margin(b = 6))
  )
save_plot("FigS2_core_genera_polished", pS2, 18, 20)
cat("FigS2 done\n")

## ─────────────────────────────────────────────
## FigS1_alpha_all: All 4 alpha-diversity metrics (redo of old FigS1)
##   ACE richness | Shannon diversity | Simpson (D) | Pielou's evenness
##   Same facet_grid style as Fig3AB; 4 rows × 4 site columns
## ─────────────────────────────────────────────
alpha_sup_long <- A_tbl |>
  select(sample, site, depth_order, ACE, Shannon, Simpson, Pielou) |>
  pivot_longer(c(ACE, Shannon, Simpson, Pielou),
               names_to = "metric", values_to = "value")

alpha_sup_long$metric <- factor(
  alpha_sup_long$metric,
  levels  = c("ACE", "Shannon", "Simpson", "Pielou"),
  labels  = c("ACE richness", "Shannon diversity", "Simpson (D)", "Pielou's evenness (J)")
)
alpha_sup_long$site       <- factor(alpha_sup_long$site, levels = c("HD","XJ","DS","DG"))
alpha_sup_long$depth_order <- factor(alpha_sup_long$depth_order,
                                      levels = c("3 cm","20 cm","40 cm"))

## Significance letters (reuse letters_from_pairwise defined above)
letter_tbl_s1 <- alpha_sup_long |>
  group_by(metric, site) |>
  group_modify(~{
    let <- letters_from_pairwise(.x, "value", group_var = "depth_order")
    data.frame(depth_order = names(let), label = unname(let))
  }) |>
  ungroup()
letter_tbl_s1$depth_order <- factor(letter_tbl_s1$depth_order,
                                     levels = c("3 cm","20 cm","40 cm"))

max_df_s1 <- alpha_sup_long |>
  group_by(metric, site, depth_order) |>
  summarise(ymax = max(value, na.rm = TRUE), .groups = "drop")

nudge_df_s1 <- alpha_sup_long |>
  group_by(metric) |>
  summarise(nudge = 0.04 * diff(range(value, na.rm = TRUE)), .groups = "drop")

lab_df_s1 <- left_join(max_df_s1, letter_tbl_s1, by = c("metric","site","depth_order")) |>
  left_join(nudge_df_s1, by = "metric") |>
  mutate(y_label = ymax + nudge)

pS1_alpha <- ggplot(alpha_sup_long,
    aes(depth_order, value, fill = depth_order, colour = depth_order)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.4,
               fatten = 1.4, alpha = 0.65) +
  geom_jitter(width = 0.14, height = 0, size = 1.7, stroke = 0.3,
              shape = 21, alpha = 0.85) +
  geom_text(data = lab_df_s1,
            aes(x = depth_order, y = y_label, label = label),
            inherit.aes = FALSE,
            size = 2.8, fontface = "bold", colour = "black", vjust = 0) +
  facet_grid(metric ~ site, scales = "free_y", switch = "y") +
  scale_fill_manual(values = depth_cols, name = "Depth") +
  scale_colour_manual(values = depth_cols, name = "Depth") +
  labs(x = NULL, y = NULL, caption = depth_method_caption) +
  theme_paper(8.5) +
  theme(
    panel.border      = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
    axis.line         = element_blank(),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 90, size = 8.5, face = "bold"),
    strip.text.x      = element_text(size = 9.5, face = "bold"),
    panel.spacing     = unit(0.5, "lines"),
    axis.text.x       = element_text(size = 7.5),
    legend.position   = "right",
    legend.key.size   = unit(0.38, "cm"),
    legend.text       = element_text(size = 7.5),
    legend.title      = element_text(size = 8, face = "bold"),
    plot.caption      = element_text(size = 6.5, hjust = 0, colour = "grey25",
                                     lineheight = 1.05, margin = margin(t = 5))
  )
save_plot("FigS1_alpha_all_polished", pS1_alpha, 22, 26)
cat("FigS1_alpha done\n")

## ─────────────────────────────────────────────
## FigS3_soil_heatmap: Soil properties × site-depth (redo of old FigS3)
##   Rows  = 12 soil variables (clustered with dendrogram)
##   Cols  = 12 site×depth mean values (site order fixed, no col clustering)
##   Top   = Site + Depth colour annotations
##   Polished: RdBu z-score, site_cols/depth_cols, white borders
## ─────────────────────────────────────────────
soil_vars_s3 <- c("SM","OM","TN","TP","TK","Nitrate_N","Ammonium_N",
                   "Olsen_P","AK","pH","EC","TDS")
soil_pretty_s3 <- c(
  SM         = "SM (%)",
  OM         = "OM (g/kg)",
  TN         = "TN (g/kg)",
  TP         = "TP (g/kg)",
  TK         = "TK (g/kg)",
  Nitrate_N  = "Nitrate-N (mg/kg)",
  Ammonium_N = "Ammonium-N (mg/kg)",
  Olsen_P    = "AP (mg/kg)",
  AK         = "AK (mg/kg)",
  pH         = "pH",
  EC         = "EC (mS/cm)",
  TDS        = "TDS (g/kg)"
)

## Mean per site × depth group
soil_grp <- meta |>
  mutate(depth_label = as.character(depth_order)) |>
  group_by(site, depth_label) |>
  summarise(across(all_of(soil_vars_s3), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  mutate(col_id = paste0(site, "\n", depth_label))

## Fixed column order: HD, XJ, DS, DG × 3cm, 20cm, 40cm
site_ord_s3  <- c("HD","XJ","DS","DG")
depth_ord_s3 <- c("3 cm","20 cm","40 cm")
col_ord_s3   <- as.vector(outer(site_ord_s3, depth_ord_s3,
                                 function(s,d) paste0(s,"\n",d)))
soil_grp <- soil_grp[match(col_ord_s3, soil_grp$col_id), ]

## Build z-score matrix: rows = soil variables, cols = site×depth
soil_mat_s3 <- as.matrix(soil_grp[, soil_vars_s3])  # [col × var]
soil_z_s3   <- scale(soil_mat_s3)                    # z-score per variable
soil_z_s3   <- t(soil_z_s3)                          # rows = var, cols = groups
rownames(soil_z_s3) <- soil_pretty_s3[soil_vars_s3]
colnames(soil_z_s3) <- soil_grp$col_id
soil_z_s3[is.na(soil_z_s3)] <- 0

## Chemical-symbol row labels (no units); order matches soil_vars_s3
row_labels_s3 <- expression(SM, OM, TN, TP, TK,
                             NO[3]^"-"*"-N", NH[4]^"+"*"-N",
                             AP, AK, pH, EC, TDS)

## Top annotations
ann_site_s3  <- gsub("\n.*", "", colnames(soil_z_s3))
ann_depth_s3 <- gsub(".*\n", "", colnames(soil_z_s3))

top_ann_s3 <- HeatmapAnnotation(
  Site  = factor(ann_site_s3,  levels = c("HD","XJ","DS","DG")),
  Depth = factor(ann_depth_s3, levels = c("3 cm","20 cm","40 cm")),
  col   = list(Site = site_cols, Depth = depth_cols),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size     = unit(0.35, "cm"),
  show_legend          = TRUE,
  annotation_legend_param = list(
    Site  = list(title = "Site",  title_gp = gpar(fontsize=8, fontface="bold"),
                 labels_gp = gpar(fontsize=7.5)),
    Depth = list(title = "Depth", title_gp = gpar(fontsize=8, fontface="bold"),
                 labels_gp = gpar(fontsize=7.5))
  )
)

col_s3 <- colorRamp2(seq(-2, 2, length.out = 101),
                     colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))

## Column names replaced by depth-only labels; site shown via column_split title
## Use bare numbers (3 / 20 / 40) to avoid label collisions; unit "(cm)" is
## appended once, after the last column-group, via decorate_column_names().
colnames(soil_z_s3) <- gsub(" cm", "", ann_depth_s3)   # "3", "20", "40" × 4 sites

ht_s3 <- Heatmap(
  soil_z_s3,
  name                   = "z-score",
  col                    = col_s3,
  top_annotation         = top_ann_s3,
  cluster_rows           = TRUE,
  cluster_columns        = FALSE,         # site order fixed
  clustering_method_rows = "ward.D2",
  show_row_dend          = TRUE,
  show_column_dend       = FALSE,
  row_dend_side          = "left",
  row_names_side         = "right",
  row_labels             = row_labels_s3,
  row_names_gp           = gpar(fontsize = 8),
  column_names_gp        = gpar(fontsize = 7.5),
  column_names_rot       = 0,
  column_split           = factor(ann_site_s3, levels = c("HD","XJ","DS","DG")),
  column_title_gp        = gpar(fontsize = 9, fontface = "bold"),
  column_gap             = unit(2, "mm"),
  rect_gp                = gpar(col = "white", lwd = 0.5),
  border                 = TRUE,
  heatmap_legend_param   = list(
    title     = "z-score",
    title_gp  = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.5),
    at        = c(-2, -1, 0, 1, 2),
    direction = "vertical"
  ),
  width  = unit(9, "cm"),
  height = unit(7, "cm")
)

pdf(file.path(fg, "FigS3_soil_heatmap_polished.pdf"),
    width = 7.7, height = 5.5, bg = "white")
draw(ht_s3, padding = unit(c(4, 8, 4, 4), "mm"),
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     merge_legend           = TRUE)
## Unit label "(cm)", once, to the right of the last (DG) column group
decorate_column_names("z-score", {
  grid.text("(cm)", x = unit(1, "npc") + unit(1, "mm"), y = unit(0.5, "npc"),
            just = "left", gp = gpar(fontsize = 7.5))
}, slice = 4)
dev.off()
cat("FigS3_soil_heatmap done\n")

## ─────────────────────────────────────────────
## Finalize: assemble the numbered publication set.
##   01_Fig1.pdf and 07_FigS2.pdf are static assets and are not distributed in
##   this sanitized code archive. If supplied under data/static/, they are
##   copied into the publication output alongside the regenerated figures.
##   All component PDFs are retained under outputs/figures/components/.
## ─────────────────────────────────────────────
## Optional static reference figures.
static_map <- c(
  "Fig1.pdf"  = "01_Fig1.pdf",
  "FigS2.pdf" = "07_FigS2.pdf"
)
for (src in names(static_map)) {
  sp <- file.path(static_dir, src)
  if (file.exists(sp)) {
    file.copy(sp, file.path(publication_dir, static_map[src]), overwrite = TRUE)
  } else {
    message("Optional static figure not supplied: ", sp)
  }
}

## Script-generated combined/final figures -> numbered submission names
## (S1->S1, S2->S2, S3->S3 unchanged; old "FigS1" fungaltraits and "FigS2"
##  core-genera figures are renumbered to S4/S5 in the new sequence)
final_map <- c(
  "Figure2_polished.pdf"            = "02_Fig2.pdf",
  "Figure3_polished.pdf"            = "03_Fig3.pdf",
  "Fig4_combined_polished.pdf"      = "04_Fig4.pdf",
  "Fig5_combined_polished.pdf"      = "05_Fig5.pdf",
  "FigS1_alpha_all_polished.pdf"    = "06_FigS1.pdf",
  "FigS3_soil_heatmap_polished.pdf" = "08_FigS3.pdf",
  "FigS1_fungaltraits_polished.pdf" = "09_FigS4.pdf",
  "FigS2_core_genera_polished.pdf"  = "10_FigS5.pdf"
)
for (src in names(final_map)) {
  copied <- file.copy(file.path(fg, src), file.path(publication_dir, final_map[src]), overwrite = TRUE)
  if (!copied) stop("Failed to assemble publication figure: ", src, call. = FALSE)
}

## Fig. S1 and Fig. S5 received a later publication-layout pass. Run that
## canonical finalizer automatically so the numbered output is not left at the
## earlier draft layout embedded above.
finalizer <- file.path(root, "scripts", "05_finalize_supplementary_figures.R")
if (!file.exists(finalizer)) stop("Missing supplementary-figure finalizer: ", finalizer, call. = FALSE)
source(finalizer, chdir = FALSE)

cat("\nScript-generated publication PDFs written to:\n", publication_dir, "\n")
cat("Component PDFs retained in:\n", fg, "\n")
