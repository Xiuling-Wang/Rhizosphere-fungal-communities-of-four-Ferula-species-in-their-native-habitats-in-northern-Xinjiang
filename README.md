# Rhizosphere fungal communities of four *Ferula* species

Publication-code archive for:

> Wang X, Liang G, Zhuang L (2026). **Rhizosphere fungal communities of
> four *Ferula* species in their native habitats in northern Xinjiang.**
> *Rhizosphere*. Manuscript identifier: RHISPH-D-26-00651R1.

This release synchronizes the repository with the final accepted analysis and
the retained publication-figure source. It replaces machine-specific paths and
revision-stage entry points with portable scripts, frozen processed inputs, and
release checks. No statistical method or published result was changed during
this cleanup.

## Quick start

Run from the repository root:

```bash
Rscript scripts/validate_release.R
Rscript scripts/04_make_publication_figures.R
```

The figure script writes PDFs to:

```text
outputs/figures/
  components/   # retained component PDFs
  publication/  # numbered publication PDFs generated from code
```

It regenerates Figs. 2–5 and Figs. S1, S3–S5 from the included processed data.
Fig. 1 and Fig. S2 are static assets and are intentionally not redistributed in
this sanitized code archive.

## Repository layout

```text
data/
  processed/    frozen non-identifying inputs and published result tables
  README.md     data provenance, exclusions, and optional-input instructions
scripts/
  01_prepare_curated_asv.R
  02_run_community_statistics.R
  03_run_traits_core_permdisp.R
  04_make_publication_figures.R
  05_finalize_supplementary_figures.R
  validate_release.R
outputs/        generated locally; ignored by Git
```

The former filenames (`run_final*.R` and `polish_figures.R`) remain as small
compatibility entry points so existing links do not break.

## Workflow and reproducibility boundary

1. `01_prepare_curated_asv.R` documents the BLAST-curated ASV filtering,
   rarefaction, composition summary, and genus-indicator workflow. It requires
   two non-distributed DADA2/BLAST intermediate files listed in
   `data/README.md`.
2. `02_run_community_statistics.R` is runnable from the included processed
   data. It reproduces alpha-diversity tests, PERMANOVA, dbRDA, hierarchical
   partitioning, and genus–soil correlations without overwriting the frozen
   publication tables.
3. `03_run_traits_core_permdisp.R` documents FungalTraits annotation, core
   genera, and PERMDISP. It requires the third-party FungalTraits source table;
   the derived publication tables are included.
4. `04_make_publication_figures.R` is the canonical final figure script and is
   runnable from the included processed data. It automatically runs stage 5.
5. `05_finalize_supplementary_figures.R` applies the later publication layouts
   for Fig. S1 and Fig. S5 using only frozen processed tables.

The raw-read-to-ASV DADA2 preprocessing workspace was not retained in the final
local archive. This repository therefore supports a verified processed-data
rerun and preserves the upstream provenance code, but it does not claim a
one-command reconstruction from FASTQ files.

Permutation p-values are Monte Carlo estimates. The script resets a documented
seed immediately before each permutation test, but the last reported digit can
still vary across R, `vegan`, or `permute` versions. The frozen values used in
the article remain in `data/processed/`; reruns should preserve effect sizes and
inferential classifications rather than be expected to reproduce every sampled
permutation p-value byte-for-byte.

The retained accepted PDFs also received a final export/crop pass whose command
was not preserved. Regenerated vector PDFs reproduce the plotted data and
labels, but some page boxes and low-level PDF metadata are therefore not
expected to be byte-identical to the accepted files.

## Statistical design retained in the code

- Thirty-five rhizosphere samples were retained; planned sample `DGAW3.2` was
  unavailable.
- The three depths were sampled within each plant. Species/site was tested at
  the plant level after aggregating depths, while depth was tested using
  within-plant restricted permutations.
- Species and site are confounded in this native-habitat survey. Results are
  interpreted as species/site associations, not isolated host-species effects.
- The frozen rarefied table contains 1,397 ASVs at 3,819 reads per sample.

## Data availability

Raw ITS1 reads are publicly available from CNCB-NGDC:

- [GSA CRA015212](https://ngdc.cncb.ac.cn/gsa/browse/CRA015212)
- [BioProject PRJCA023762](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA023762)

The repository includes only compact processed inputs and result tables needed
for code inspection and figure reproduction. Raw reads, exact coordinates,
author/contact files, manuscripts, submission correspondence, Zotero audit
logs, provider documents, and local machine paths are excluded.

## R dependencies

The scripts use R plus the following packages where applicable:

`tidyverse`, `vegan`, `permute`, `ape`, `cowplot`, `RColorBrewer`, `labdsv`,
`VennDiagram`, `futile.logger`, `ggplot2`, `ggrepel`, `pheatmap`, `rdacca.hp`,
`reshape2`, `ComplexHeatmap`, `circlize`, `corrplot`, `ggplotify`, and `ggtext`.

PDF export uses Cairo and Helvetica. On systems without Helvetica, configure an
equivalent installed sans-serif font before comparing typography.

## Release validation

`scripts/validate_release.R` checks:

- sample, ASV, taxonomy, and rarefaction dimensions;
- identifier alignment and absence of the removed raw-file path columns;
- published PERMANOVA and stringent-core guard values;
- R syntax for every released script; and
- local-path, email, private-key, and common GitHub-token patterns.
