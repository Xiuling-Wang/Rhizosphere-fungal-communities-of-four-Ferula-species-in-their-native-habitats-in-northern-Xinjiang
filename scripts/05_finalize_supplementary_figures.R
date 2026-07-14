## Final publication layouts for Fig. S1 and Fig. S5.
## These two panels were revised after the main figure script was written.
options(stringsAsFactors = FALSE)
set.seed(20260615)

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggtext)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
tb <- file.path(root, "data", "processed")
publication_dir <- file.path(root, "outputs", "figures", "publication")
dir.create(publication_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  metadata = file.path(tb, "sample_metadata.tsv"),
  asv_counts = file.path(tb, "asv_table_fungi_rarefied_3819.tsv"),
  alpha_diversity = file.path(tb, "alpha_diversity_final.tsv"),
  stringent_core = file.path(tb, "core_genera_stringent_depth_site.tsv"),
  four_site_core = file.path(tb, "core_genera_four_site_shared.tsv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing processed input files:\n- ", paste(missing_inputs, collapse = "\n- "), call. = FALSE)
}

site_cols <- c(HD = "#4DBBD5", XJ = "#E64B35", DS = "#00A087", DG = "#925E9F")
depth_alpha <- c("3" = 0.45, "20" = 0.70, "40" = 1.00)

theme_paper <- function(base_size = 8, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.title = element_text(face = "bold")
    )
}

save_pdf <- function(plot, filename, width_cm, height_cm) {
  ggsave(filename, plot, width = width_cm, height = height_cm,
         units = "cm", bg = "white", device = cairo_pdf)
}

letters_from_pairwise <- function(df, metric, group_var = "depth_order", alpha = 0.05) {
  groups <- levels(df[[group_var]])
  y <- df[[metric]]
  g <- df[[group_var]]
  if (length(unique(g[!is.na(y)])) < 2) return(setNames(rep("a", length(groups)), groups))
  kw <- suppressWarnings(kruskal.test(y ~ g))
  if (is.na(kw$p.value) || kw$p.value >= alpha) return(setNames(rep("", length(groups)), groups))
  pw <- suppressWarnings(pairwise.wilcox.test(y, g, p.adjust.method = "BH", exact = FALSE)$p.value)
  means <- tapply(y, g, median, na.rm = TRUE)
  ordered_groups <- names(sort(means, decreasing = TRUE))
  labels <- setNames(rep("", length(groups)), groups)
  for (group in ordered_groups) {
    for (letter in letters[seq_len(8)]) {
      members <- names(labels)[grepl(letter, labels, fixed = TRUE)]
      conflict <- FALSE
      for (member in members) {
        p <- if (group %in% rownames(pw) && member %in% colnames(pw)) pw[group, member]
        else if (member %in% rownames(pw) && group %in% colnames(pw)) pw[member, group]
        else NA_real_
        if (!is.na(p) && p < alpha) conflict <- TRUE
      }
      if (!conflict) {
        labels[group] <- paste0(labels[group], letter)
        break
      }
    }
  }
  labels[labels == ""] <- "a"
  labels
}

## Fig. S1: additional alpha-diversity indices.
meta <- read.delim(inputs[["metadata"]], check.names = FALSE)
meta$site <- factor(meta$site, levels = c("HD", "XJ", "DS", "DG"))
alpha <- read.delim(inputs[["alpha_diversity"]], check.names = FALSE)
alpha <- merge(alpha, meta[, c("sample", "site")], by = "sample")
alpha$site <- factor(alpha$site, levels = c("HD", "XJ", "DS", "DG"))
alpha$depth_order <- factor(alpha$depth_order, levels = c(1, 2, 3), labels = c("3", "20", "40"))

alpha_long <- alpha |>
  select(sample, site, depth_order, ACE, Simpson, Pielou) |>
  pivot_longer(c(ACE, Simpson, Pielou), names_to = "metric", values_to = "value")
alpha_long$metric <- factor(
  alpha_long$metric,
  levels = c("ACE", "Simpson", "Pielou"),
  labels = c("ACE richness", "Simpson (D)", "Pielou's evenness (J)")
)

letter_tbl <- alpha_long |>
  group_by(metric, site) |>
  group_modify(~data.frame(
    depth_order = names(letters_from_pairwise(.x, "value")),
    label = unname(letters_from_pairwise(.x, "value"))
  )) |>
  ungroup()
letter_tbl$depth_order <- factor(letter_tbl$depth_order, levels = c("3", "20", "40"))
max_tbl <- alpha_long |>
  group_by(metric, site, depth_order) |>
  summarise(ymax = max(value, na.rm = TRUE), .groups = "drop")
nudge_tbl <- alpha_long |>
  group_by(metric) |>
  summarise(nudge = 0.04 * diff(range(value, na.rm = TRUE)), .groups = "drop")
label_tbl <- max_tbl |>
  left_join(letter_tbl, by = c("metric", "site", "depth_order")) |>
  left_join(nudge_tbl, by = "metric") |>
  mutate(y_label = ymax + nudge)

fig_s1 <- ggplot(alpha_long, aes(depth_order, value, fill = site, colour = site, alpha = depth_order)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.40, fatten = 1.35) +
  geom_jitter(width = 0.14, height = 0, size = 1.65, stroke = 0.30, shape = 21) +
  geom_text(data = label_tbl, aes(depth_order, y_label, label = label), inherit.aes = FALSE,
            size = 2.6, fontface = "bold", colour = "black", vjust = 0) +
  facet_grid(metric ~ site, scales = "free_y", switch = "y") +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_colour_manual(values = site_cols, guide = "none") +
  scale_alpha_manual(values = depth_alpha, name = "Depth") +
  labs(x = NULL, y = NULL,
       caption = "No significant differences among species or depths (Kruskal-Wallis, FDR-corrected).") +
  theme_paper(8.4) +
  theme(
    panel.border = element_rect(colour = "grey38", fill = NA, linewidth = 0.48),
    axis.line = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 8.6, face = "bold"),
    strip.text.x = element_text(size = 9.5, face = "bold"),
    panel.spacing = unit(0.52, "lines"),
    axis.text.x = element_text(size = 7.4),
    legend.position = "bottom",
    legend.key.size = unit(0.40, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7.4, face = "bold"),
    plot.caption = element_text(size = 9.5, hjust = 0, colour = "grey20",
                                lineheight = 1.05, margin = margin(t = 7))
  ) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey35", colour = "grey35", size = 3)))

save_pdf(fig_s1, file.path(publication_dir, "06_FigS1.pdf"), 22, 20)

## Fig. S5: strict and four-site core genera.
strict_core <- read.delim(inputs[["stringent_core"]], check.names = FALSE)
four_site_core <- read.delim(inputs[["four_site_core"]], check.names = FALSE)
rare_tbl <- read.delim(inputs[["asv_counts"]], check.names = FALSE)
rare <- as.matrix(rare_tbl[, -1, drop = FALSE])
rownames(rare) <- rare_tbl$sample
meta_core <- meta[match(rownames(rare), meta$sample), ]
present <- rare > 0
four_site_asvs <- Reduce(intersect, lapply(levels(meta_core$site), function(site) {
  colnames(present)[colSums(present[meta_core$site == site, , drop = FALSE]) > 0]
}))
if (length(four_site_asvs) != 68L) stop("Unexpected four-site core size: ", length(four_site_asvs), call. = FALSE)

lifestyle_labels <- c(
  arbuscular_mycorrhizal = "Arbuscular mycorrhizal",
  unspecified_saprotroph = "Unspecified saprotroph",
  soil_saprotroph = "Soil saprotroph",
  wood_saprotroph = "Wood saprotroph",
  litter_saprotroph = "Litter saprotroph",
  plant_pathogen = "Plant pathogen",
  sooty_mold = "Sooty mold",
  Unassigned = "Unassigned"
)
lifestyle_cols <- c(
  arbuscular_mycorrhizal = "#238B45",
  unspecified_saprotroph = "#86BBD8",
  soil_saprotroph = "#2F80B7",
  wood_saprotroph = "#E8A23C",
  litter_saprotroph = "#F0C27B",
  plant_pathogen = "#C44E52",
  sooty_mold = "#8E6BBE",
  Unassigned = "#D9D9D9"
)
lifestyle_order <- names(lifestyle_labels)

strict_plot <- strict_core |>
  mutate(primary_lifestyle = factor(primary_lifestyle, levels = lifestyle_order)) |>
  arrange(desc(n_core_ASVs), primary_lifestyle, genus) |>
  mutate(genus_label = ifelse(genus == "Unclassified", "Unclassified", paste0("<i>", genus, "</i>")),
         genus_f = factor(genus_label, levels = rev(genus_label)))
shared_plot <- four_site_core |>
  mutate(primary_lifestyle = factor(primary_lifestyle, levels = lifestyle_order)) |>
  arrange(desc(n_core_ASVs), primary_lifestyle, genus) |>
  mutate(genus_f = factor(genus, levels = rev(genus)))

p_strict <- ggplot(strict_plot, aes(n_core_ASVs, genus_f, colour = primary_lifestyle)) +
  geom_segment(aes(x = 0, xend = n_core_ASVs, yend = genus_f),
               linewidth = 0.55, alpha = 0.42, colour = "grey55") +
  geom_point(size = 3.3) +
  geom_text(aes(label = n_core_ASVs), nudge_x = 0.12, size = 2.9, colour = "grey20") +
  scale_colour_manual(values = lifestyle_cols, breaks = lifestyle_order,
                      labels = lifestyle_labels, drop = TRUE, name = "Primary lifestyle") +
  scale_x_continuous(limits = c(0, max(strict_plot$n_core_ASVs) + 0.45),
                     breaks = seq(0, max(strict_plot$n_core_ASVs), by = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = "No. of core ASVs", y = NULL, title = "Strict core",
       subtitle = paste0(sum(strict_plot$n_core_ASVs), " ASVs; ",
                         sum(strict_plot$genus != "Unclassified"), " named genera (+",
                         strict_plot$n_core_ASVs[strict_plot$genus == "Unclassified"], " unassigned)")) +
  theme_paper(8) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 9.4),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.28),
    plot.title = element_text(size = 8.6, face = "bold"),
    plot.subtitle = element_text(size = 7.2, colour = "grey30"),
    legend.position = "none"
  )

p_shared <- ggplot(shared_plot, aes(n_core_ASVs, genus_f, colour = primary_lifestyle)) +
  geom_segment(aes(x = 0, xend = n_core_ASVs, yend = genus_f),
               linewidth = 0.45, alpha = 0.38, colour = "grey58") +
  geom_point(aes(size = mean_relative_abundance), alpha = 0.95) +
  geom_text(aes(label = n_core_ASVs), nudge_x = 0.18, size = 2.25, colour = "grey20") +
  scale_colour_manual(values = lifestyle_cols, breaks = lifestyle_order,
                      labels = lifestyle_labels, drop = TRUE, name = "Primary lifestyle") +
  scale_size_continuous(range = c(1.8, 4.6), name = "Mean RA (%)") +
  scale_x_continuous(limits = c(0, max(shared_plot$n_core_ASVs) + 0.8),
                     breaks = seq(0, max(shared_plot$n_core_ASVs), by = 2),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "No. of shared ASVs", y = NULL, title = "Four-site shared core",
       subtitle = paste0(length(four_site_asvs), " ASVs; ", nrow(shared_plot), " named genera")) +
  theme_paper(8) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 7.8),
    panel.grid.major.x = element_line(colour = "grey91", linewidth = 0.24),
    plot.title = element_text(size = 8.6, face = "bold"),
    plot.subtitle = element_text(size = 7.2, colour = "grey30"),
    legend.position = "right", legend.key.size = unit(0.34, "cm"),
    legend.text = element_text(size = 6.7),
    legend.title = element_text(size = 7.2, face = "bold")
  ) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3)),
         size = guide_legend(order = 2))

fig_s5_main <- plot_grid(
  p_strict, p_shared, labels = c("A", "B"), ncol = 2,
  rel_widths = c(0.82, 1.55), label_fontface = "bold", label_size = 18,
  label_x = 0.006, label_y = 0.996, hjust = 0, vjust = 1,
  align = "h", axis = "tb"
)
caption_panel <- ggdraw() +
  draw_label(
    "Strict core, present in at least half of the samples at every site (matches Table S2);\nFour-site core, detected in all four sites.",
    x = 0.01, y = 0.70, hjust = 0, vjust = 1, size = 9.5,
    colour = "grey20", lineheight = 1.05
  )
fig_s5 <- plot_grid(fig_s5_main, caption_panel, ncol = 1, rel_heights = c(1, 0.075))
save_pdf(fig_s5, file.path(publication_dir, "10_FigS5.pdf"), 28, 20)

cat("Finalized supplementary PDFs:\n")
cat(file.path(publication_dir, "06_FigS1.pdf"), "\n")
cat(file.path(publication_dir, "10_FigS5.pdf"), "\n")
