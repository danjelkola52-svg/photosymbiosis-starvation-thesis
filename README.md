# Photosymbiosis under food deprivation: *Exaiptasia* vs. *Waminoa*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23264022.svg)](https://doi.org/10.5281/zenodo.23264022)

Code, data and supplementary tables for the Master's thesis **"Resolving Starvation Responses of Photosymbiotic Holobionts"**
(Danjel Kola, M.Sc. Microbiology, University of Oldenburg; thesis work at HIFMB, supervisor Dr. Nils Rädecker).

The thesis compares how fed and starved *Exaiptasia diaphana* (sea anemone) and *Waminoa* sp. (acoel flatworm) regulate their
algal symbiosis (symbiont density, oxygen production, photosynthetic efficiency, ammonium uptake, host protein) and
characterises their bacterial communities by 16S rRNA amplicon sequencing (V5–V7, primers 799F/1193R).

Everything in the thesis can be reproduced from **two R scripts**:

| Script | What it does | Run time |
|---|---|---|
| [`scripts/physiology_analysis.R`](scripts/physiology_analysis.R) | physiology: assumptions, Welch t-tests (Holm), two-way ANOVA, Table A1, Figures 4–11 | < 1 min |
| [`scripts/16S_amplicon_analysis.R`](scripts/16S_amplicon_analysis.R) | 16S: DADA2 → SILVA → contaminant removal → diversity, PERMANOVA, ANCOM-BC2, Figure 12, Tables B1/B2 | 10–20 min (+ 2–3 h if you start from raw reads) |

## Quick start

```r
# from the repository root
install.packages(c("readxl","dplyr","tidyr","purrr","tibble","stringr","readr","ggplot2","patchwork","car","openxlsx",
                   "vegan","ggrepel","ragg"))
BiocManager::install(c("dada2","phyloseq","Biostrings","ShortRead","ANCOMBC","decontam"))

source("scripts/physiology_analysis.R")      # -> results/physiology/{tables,figures}
source("scripts/16S_amplicon_analysis.R")    # -> 16S/{tables,figures}   (skips DADA2 if 16S/dada2_out/ps_raw.rds exists)
```

Software used for the thesis: R 4.4.3, DADA2 1.34, ANCOMBC 2.8.0, SILVA 138.2 (`16S/R_session_info.txt`,
`16S/dada2_out/dada2_settings_and_session.txt`).

## Repository layout

```
scripts/
  physiology_analysis.R          physiology statistics, Table A1, Figures 4-11
  16S_amplicon_analysis.R        16S pipeline (Part 1 DADA2, Part 2 analysis)
  figure_scripts_python/         scripts that drew the printed Figures 3-11 (matplotlib); same statistics as the R script
data/physiology/                 master_data_all_parameters.xlsx   (22 animals x 13 variables, sheet "Master")
16S/dada2_out/                   ps_raw.rds (phyloseq object = input of Part 2), ASV taxonomy, read tracking, primer-dimer QC
16S/tables/16S_all_tables.xlsx   every 16S table (36 sheets), thesis analysis n = 17 animals
16S/figures/                     16S figures (PDF, vector)
16S/reference_db/                put the SILVA files here (see below)
results/physiology/              output of physiology_analysis.R
supplement/                      Supplementary_Tables_and_Figure_Data.xlsx  (what the thesis cites as "Supplementary workbook")
thesis/                          the thesis as submitted (Word)
```

## What reproduces what

| Thesis item | Produced by | Where |
|---|---|---|
| Table A1 (means ± SD, change, Holm p, interaction) | `physiology_analysis.R` | `results/physiology/tables/table_a1.csv`, sheet `Table_A1` of the workbook |
| Figures 4–11 (Fed vs Starved) | `physiology_analysis.R` (R version) / `figure_scripts_python/thesis_physiology_figures.py` (printed version) | `results/physiology/figures/` |
| Assumption tests (Shapiro–Wilk, Levene) | `physiology_analysis.R` | `assumptions.csv` |
| Contamination removal (75.8 % of reads) | `16S_amplicon_analysis.R`, Part 2 | `16S/tables/16S_all_tables.xlsx` (`ASV_filter_summary`, `Contamination_per_sample`, `ASV_all_with_status`) |
| Table B1 (PERMANOVA) | Part 2 | sheet `PERMANOVA` |
| Table B2 (PERMANOVA with primer-dimer load first) | Part 2 | sheet `PERMANOVA_dimer_adj` |
| Figure 12 A/B/C (families, PCoA, ANCOM-BC2) | Part 2 | `16S/figures/main/Figure_main_ABC_BrayCurtis.*` |

Thesis analysis = **17 animals**: five libraries (AF5, AS1, AS2, WF2, WS4; 870–1,717 reads after cleaning, 73–93 % primer dimer)
are excluded. `EXCLUDE_ANIMALS=""` runs all 20 animals (sensitivity analysis). The permutation tests are seeded, so the
PERMANOVA numbers reproduce exactly.

## Contents of the supplementary workbook

`supplement/Supplementary_Tables_and_Figure_Data.xlsx` has 37 sheets (one per table or figure, plus an index; the `README` sheet links to every sheet). A figure sheet shows the figure as printed together with its adjusted p-values and the per-animal values.

| Sheets | Content | Used in |
|---|---|---|
| `Table_A1`, `Phys_Contrasts`, `Phys_ANOVA`, `Phys_Assumptions`, `Phys_Descriptives` | Table A1; fed vs. starved tests with exact permutation p-values and Hedges' g; two-way ANOVA; Shapiro-Wilk and Levene tests; group means | Sections 2.9, 3.1; Appendix A |
| `Phys_Data`, `Data_dictionary` | Values for every animal and parameter; meaning of each column of the physiology workbook | Figures 4-11 |
| `Fig04_TotalProtein` to `Fig11_O2perCell` (8 sheets) | Each figure as printed, with its adjusted p-values and the per-animal data | Figures 4-11 |
| `Fig12_16S`, `Fig12A_Family_bars`, `Fig12B_PCoA_BrayCurtis`, `Fig12C_ANCOMBC2_family` | Figure 12 and the data behind panels A-C | Figure 12 |
| `Table_B1_PERMANOVA`, `Table_B1b_Pairwise`, `Table_B2_PERMANOVA_dimer` | PERMANOVA of community composition; pairwise tests; model with primer-dimer load as first term | Section 3.2; Appendix B |
| `16S_Samples`, `16S_ReadTracking`, `16S_Cleaning_Steps`, `16S_NEGCON1`, `16S_Contam_per_sample` | Reads per library through DADA2 and cleaning; primer-dimer load; reads removed per contaminant rule; role of the extraction blank | Section 2.8; Appendix B.2-B.4 |
| `16S_ASV_status_log`, `16S_ASV_counts_clean` | Status of every ASV (retained, or removed by which rule); decontaminated ASV count table with taxonomy | Appendix B.3 |
| `16S_Alpha_diversity`, `16S_ANCOMBC2_all_levels`, `16S_Top15_Genus_after`, `16S_Top15_Family_before`, `16S_Phylum_percent`, `16S_Family_percent`, `16S_Genus_percent` | Alpha diversity; complete ANCOM-BC2 results; taxonomic composition per animal | Sections 3.2, 4.5 |

## Raw data and reference databases

* **Raw 16S reads** (25 libraries, paired-end 2 × 250 bp, BMK/Novogene run CS705-003): ENA/SRA accession **[ACCESSION TO BE ADDED]**.
  Place the 50 `.fq` files in `16S/raw_fastq/` only if you want to rerun Part 1 (DADA2).
* **SILVA 138.2** (`silva_nr99_v138.2_toGenus_trainset.fa.gz`, `silva_v138.2_assignSpecies.fa.gz`): <https://doi.org/10.5281/zenodo.14169026>.
  Optional, for the 138.1 comparison sheet: `silva_nr99_v138.1_train_set.fa.gz`, `silva_species_assignment_v138.1.fa.gz`, <https://doi.org/10.5281/zenodo.4587955>.
  Put them in `16S/reference_db/`.

## Notes

* The workbook still labels *Exaiptasia* as "Aiptasia" (original sample sheets); the scripts relabel it.
* The in-silico primer-coverage check (thesis Appendix B.1: TestPrime against SILVA and a 100,000-sequence sample of the SILVA 138.1 training set) was run outside this repository.
* Intermediate DADA2 checkpoints and the 2.5 GB of raw reads are not in the repository.

## License and citation

Code: MIT (`LICENSE`). Data, tables and figures: CC BY 4.0.
Please cite the archived repository: Kola, D. (2026). *Resolving Starvation Responses of Photosymbiotic Holobionts* (Master's thesis, Carl von Ossietzky Universität Oldenburg): code, data and supplementary material. Zenodo. <https://doi.org/10.5281/zenodo.23264022> (concept DOI; always resolves to the latest version). See also `CITATION.cff`.
