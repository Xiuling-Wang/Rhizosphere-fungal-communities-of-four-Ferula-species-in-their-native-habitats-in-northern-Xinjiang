## Figure polishing entry script for the revised Rhizosphere submission.
## The statistical analyses are produced by run_final*.R; this script provides
## the shared colours, theme, and file checks used when regenerating polished figures.
options(stringsAsFactors = FALSE)
set.seed(20260609)

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
tb <- file.path(root, "tables")
fg <- file.path(root, "figures")
dir.create(fg, recursive = TRUE, showWarnings = FALSE)

site_cols <- c(HD = "#4DBBD5", XJ = "#E64B35", DS = "#00A087", DG = "#925E9F")
depth_cols <- setNames(RColorBrewer::brewer.pal(9, "Blues")[c(3, 6, 9)],
                       c("3 cm", "20 cm", "40 cm"))

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
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

save_plot <- function(name, plot, width, height) {
  ggsave(file.path(fg, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "cm",
         bg = "white", device = cairo_pdf)
}

required_tables <- c(
  "asv_table_fungi_rarefied_3819.tsv",
  "sample_metadata_downstream.tsv",
  "alpha_diversity_final.tsv",
  "genus_biomarkers_final.tsv",
  "core_genera.tsv"
)
missing_tables <- required_tables[!file.exists(file.path(tb, required_tables))]
if (length(missing_tables)) {
  stop("Missing required final tables: ", paste(missing_tables, collapse = ", "),
       "\nRun scripts/run_final.R, scripts/run_final_part2.R, and scripts/run_final_part3.R first.")
}

message("Figure polishing helpers loaded successfully.")
message("Use site_cols, depth_cols, theme_paper(), and save_plot() to regenerate journal-style panels from final tables.")
message("Final tables found in: ", tb)
message("Figures will be written to: ", fg)
