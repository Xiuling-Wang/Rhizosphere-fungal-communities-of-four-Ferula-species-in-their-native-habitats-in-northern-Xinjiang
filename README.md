# Rhizosphere fungal communities of four Ferula species

This repository contains the final R scripts used for the revised downstream analysis and figure polishing of the manuscript:

**Rhizosphere fungal communities of four Ferula species in their native habitats in northern Xinjiang**

The uploaded scripts correspond to the BLAST-curated fungal ASV workflow used for the revision. Legacy exploratory scripts based on earlier non-final filtering, Mantel tests, Levins niche breadth, and generalist/specialist analyses were intentionally excluded to avoid confusion.

## Workflow

Run the scripts from a project directory that contains the expected `analysis/`, `tables/`, and `figures/` folders:

```bash
Rscript scripts/run_final.R /path/to/project_root
Rscript scripts/run_final_part2.R /path/to/project_root
Rscript scripts/run_final_part3.R /path/to/project_root
Rscript scripts/polish_figures.R /path/to/project_root
```

If no path is supplied, each script treats the current working directory as the project root.

## Script roles

- `scripts/run_final.R`: builds the BLAST-curated fungal ASV table, rarefies reads, generates Figure 2 components, and calculates genus-level biomarkers using Kruskal-Wallis/BH plus IndVal.
- `scripts/run_final_part2.R`: reproduces alpha diversity, class heatmap, dbRDA/PERMANOVA, variation partitioning, and genus-soil correlation analyses for Figures 3-5.
- `scripts/run_final_part3.R`: adds FungalTraits functional annotation, core genera shared by all four Ferula species, and PERMDISP checks.
- `scripts/polish_figures.R`: creates the final polished publication figures from the final analysis tables.

## Expected inputs

The scripts expect the corrected local analysis outputs to be arranged as follows:

```text
project_root/
  analysis/
    asv_trial/
      outputs/
        asv_table_nonchim.tsv
      taxonomy/
        fungal_set_BLASTrebuilt.txt
        asv_taxonomy_fungiBLASTset_UNITE.tsv
        fungaltraits/
          polme2020_genera.csv
    downstream_fungi_final/
      tables/
        sample_metadata_downstream.tsv
  tables/
    asv_table_fungi_rarefied_3819.tsv
    sample_metadata_downstream.tsv
    genus_biomarkers_final.tsv
    ...
  figures/
```

Raw FASTQ files and large intermediate files are not committed here.

## Main revision choices

- Used DADA2-style ASV-level downstream tables rather than the older OTU table.
- Rebuilt the fungal set using BLAST-informed curation because both EUKARYOME and UNITE showed errors near the plant/fungal boundary.
- Removed confirmed host/plant-derived ASVs before downstream analysis.
- Retained only assigned fungal phyla for composition summaries.
- Used UNITE taxonomy for fungal annotation.
- Used FungalTraits (Polme et al. 2020) for genus-level lifestyle annotation, including arbuscular mycorrhizal fungi, saprotrophs, and plant pathogens.
- Added core genera and PERMDISP as concise reviewer-response analyses.

## R packages

The scripts use: `tidyverse`, `vegan`, `cowplot`, `RColorBrewer`, `labdsv`, `VennDiagram`, `ape`, `ggplot2`, `ggrepel`, `pheatmap`, `rdacca.hp`, `reshape2`, `ComplexHeatmap`, `circlize`, `corrplot`, `ggplotify`, `ggsci`, `ggtext`, `png`, and `grid`.

## References for analysis resources

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581-583.
- Wang Q, Garrity GM, Tiedje JM, Cole JR. 2007. Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied and Environmental Microbiology* 73:5261-5267.
- Abarenkov K et al. 2024. UNITE general FASTA release for fungi.
- Polme S et al. 2020. FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. *Fungal Diversity* 105:1-16.
