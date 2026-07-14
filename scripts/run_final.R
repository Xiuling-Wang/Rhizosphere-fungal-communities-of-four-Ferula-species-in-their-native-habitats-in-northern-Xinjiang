## Compatibility entry point retained for links to the original repository layout.
## Canonical script: scripts/01_prepare_curated_asv.R
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(file_arg)) stop("Run this compatibility entry point with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
source(file.path(script_dir, "01_prepare_curated_asv.R"), chdir = FALSE)
