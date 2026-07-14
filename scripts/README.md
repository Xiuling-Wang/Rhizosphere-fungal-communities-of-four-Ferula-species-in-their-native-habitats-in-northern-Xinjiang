# Script order

Run commands from the repository root. Each script also accepts the repository
root as its first argument.

| Stage | Canonical script | Default status |
| --- | --- | --- |
| Validate | `validate_release.R` | Runnable with base R |
| 1 | `01_prepare_curated_asv.R` | Requires non-distributed DADA2/BLAST intermediates |
| 2 | `02_run_community_statistics.R` | Runnable from included processed data |
| 3 | `03_run_traits_core_permdisp.R` | Requires the third-party FungalTraits source table |
| 4 | `04_make_publication_figures.R` | Runnable from included processed data; calls stage 5 |
| 5 | `05_finalize_supplementary_figures.R` | Final Fig. S1/Fig. S5 layouts; also runnable alone |

Recommended public-release checks:

```bash
Rscript scripts/validate_release.R
Rscript scripts/02_run_community_statistics.R
Rscript scripts/04_make_publication_figures.R
```

Stages 1–3 write only under `outputs/reanalysis/`. Stages 4–5 write only under
`outputs/figures/`. Frozen files under `data/processed/` are never overwritten.

The legacy filenames retained at this level are compatibility wrappers for the
canonical numbered scripts. New work should cite and use the numbered files.
