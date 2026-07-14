## Validate the sanitized publication-code release without recalculating results.
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
data_dir <- file.path(root, "data", "processed")

required <- file.path(data_dir, c(
  "asv_table_fungi_rarefied_3819.tsv",
  "asv_taxonomy_fungi.tsv",
  "sample_metadata.tsv",
  "alpha_diversity_final.tsv",
  "permanova_final.tsv",
  "core_genera.tsv",
  "core_genera_stringent_depth_site.tsv"
))
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing release files:\n- ", paste(missing, collapse = "\n- "), call. = FALSE)

rr <- read.delim(required[[1]], check.names = FALSE)
rare <- as.matrix(rr[, -1, drop = FALSE])
rownames(rare) <- rr$sample
storage.mode(rare) <- "numeric"

tax <- read.delim(required[[2]], check.names = FALSE, quote = "")
meta <- read.delim(required[[3]], check.names = FALSE, quote = "")
alpha <- read.delim(required[[4]], check.names = FALSE)
permanova <- read.delim(required[[5]], row.names = 1, check.names = FALSE)
core_named <- read.delim(required[[6]], check.names = FALSE)
core_all <- read.delim(required[[7]], check.names = FALSE)
present <- rare > 0
four_site_asvs <- Reduce(intersect, lapply(c("HD", "XJ", "DS", "DG"), function(site) {
  colnames(present)[colSums(present[meta$site == site, , drop = FALSE]) > 0]
}))

stopifnot(
  nrow(rare) == 35L,
  ncol(rare) == 1397L,
  all(rowSums(rare) == 3819),
  !anyDuplicated(rownames(rare)),
  nrow(tax) == 1761L,
  !anyDuplicated(tax$asv_id),
  all(colnames(rare) %in% tax$asv_id),
  nrow(meta) == 35L,
  identical(rownames(rare), meta$sample),
  !any(c("raw_R1", "raw_R2") %in% names(meta)),
  identical(as.integer(table(meta$site)[c("HD", "XJ", "DS", "DG")]), c(9L, 9L, 9L, 8L)),
  identical(as.integer(table(meta$depth_order)[c("1", "2", "3")]), c(12L, 12L, 11L)),
  nrow(alpha) == 35L,
  isTRUE(all.equal(permanova["site", "R2"], 0.36866, tolerance = 1e-8)),
  isTRUE(all.equal(permanova["site", "Pr(>F)"], 0.001, tolerance = 1e-8)),
  isTRUE(all.equal(permanova["depth_order", "R2"], 0.06096, tolerance = 1e-8)),
  isTRUE(all.equal(permanova["depth_order", "Pr(>F)"], 0.036, tolerance = 1e-8)),
  nrow(core_named) == 4L,
  sum(core_named$n_core_ASVs) == 4L,
  sum(core_all$n_core_ASVs) == 6L,
  core_all$n_core_ASVs[core_all$genus == "Unclassified"] == 2L,
  length(four_site_asvs) == 68L
)

script_files <- list.files(file.path(root, "scripts"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
for (script in script_files) parse(script)

text_files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
text_files <- text_files[
  !grepl("/(\\.git|outputs)/", text_files) &
  grepl("([.]R|[.]md|[.]tsv|[.]cff|[.]gitignore)$", text_files)
]
unsafe_patterns <- c(
  local_user_path = paste0("/", "Users", "/"),
  github_token = "gh[opusr]_[A-Za-z0-9_]+",
  fine_grained_token = "github_pat_[A-Za-z0-9_]+",
  private_key = "BEGIN (RSA |OPENSSH )?PRIVATE KEY",
  email_address = "[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+[.][A-Za-z]{2,}"
)
for (path in text_files) {
  content <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  hits <- names(unsafe_patterns)[vapply(unsafe_patterns, grepl, logical(1), x = content, perl = TRUE)]
  if (length(hits)) stop("Privacy scan failed for ", path, ": ", paste(hits, collapse = ", "), call. = FALSE)
}

cat("Release validation passed.\n")
cat(sprintf("- %d samples; %d rarefied ASVs; %d reads per sample\n", nrow(rare), ncol(rare), unique(rowSums(rare))))
cat(sprintf("- %d taxonomy records; %d R scripts parsed\n", nrow(tax), length(script_files)))
cat("- Published PERMANOVA, stringent-core, and four-site-core guard values matched\n")
cat(sprintf("- Privacy scan passed for %d tracked-text candidates\n", length(text_files)))
