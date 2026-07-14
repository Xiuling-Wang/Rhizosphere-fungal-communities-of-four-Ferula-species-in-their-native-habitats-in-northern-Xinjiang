# Data included in this code archive

`processed/` contains the frozen, non-identifying inputs and result tables used
by the publication scripts. These files are derived from the accepted analysis
and supplementary workbook; they do not contain author contact details,
submission correspondence, exact sampling coordinates, local filesystem paths,
or raw FASTQ files.

Key files:

- `asv_table_fungi_rarefied_3819.tsv`: 35 samples by 1,397 retained ASVs,
  rarefied to 3,819 reads per sample.
- `asv_taxonomy_fungi.tsv`: UNITE v10.0 taxonomy for the 1,761 retained fungal
  ASVs; unassigned ranks are represented by empty fields.
- `sample_metadata.tsv`: sample codes, species/site and depth design, and soil
  measurements. The original `raw_R1` and `raw_R2` machine paths were removed.
- `*_final.tsv` and the other small tables: frozen published statistics and
  processed inputs used by the publication-figure script.

One planned sample (`DGAW3.2`) was unavailable and is not present in the
35-sample downstream dataset.

Soil-variable units follow the article: `SM` (%); `OM`, `TN`, `TP`, and `TK`
(g/kg); `Nitrate_N`, `Ammonium_N`, `Olsen_P`, and `AK` (mg/kg); `EC` (mS/cm);
`TDS` (g/kg); and `pH` (unitless).

## Data not redistributed here

Raw ITS1 reads are available from the Genome Sequence Archive at CNCB-NGDC:

- GSA accession: `CRA015212`
- BioProject: `PRJCA023762`

To run `scripts/01_prepare_curated_asv.R`, create `data/intermediate/` and add:

- `asv_table_nonchim.tsv`
- `fungal_set_BLASTrebuilt.txt`

These DADA2/BLAST intermediate files were not retained in the final local
archive and are therefore not claimed to be downloadable from this repository.

To run `scripts/03_run_traits_core_permdisp.R`, obtain the FungalTraits source
table described by Polme et al. (2020) and place `polme2020_genera.csv` in
`data/external/`. The third-party source table is not redistributed here; the
derived publication tables are included in `processed/`.

`scripts/04_make_publication_figures.R` regenerates the code-produced figures
without either of those optional inputs. Figure 1 and Fig. S2 are static assets;
if needed for a complete numbered set, place `Fig1.pdf` and `FigS2.pdf` in
`data/static/` before running the script.
