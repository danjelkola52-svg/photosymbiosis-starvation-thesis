# Starvation responses of photosymbiotic holobionts: *Exaiptasia* vs. *Waminoa*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23264022.svg)](https://doi.org/10.5281/zenodo.23264022)

Code, data and supplementary tables for the Master's thesis **"Resolving Starvation Responses of Photosymbiotic Holobionts"**
(Danjel Kola, M.Sc. Microbiology, Carl von Ossietzky Universität Oldenburg; thesis work at HIFMB, supervisor Dr. Nils Rädecker).
In the thesis this repository is cited as **Kola (2026)**; supplementary tables are cited by the sheet name in the supplementary workbook,
for example *(Kola, 2026; sheet Table_A1)*. The thesis cites **version 1.1.0**: <https://doi.org/10.5281/zenodo.23337121>.

The thesis compares fed and starved *Exaiptasia diaphana* (sea anemone) and *Waminoa* sp. (acoel flatworm):
symbiont density, oxygen production, photosynthetic efficiency (Fv/Fm), ammonium uptake and host protein, plus the
bacterial communities of both hosts (16S rRNA gene, V5–V7, primers 799F/1193R).

---

## What is in the repository

```
supplement/Supplementary_Tables_and_Figure_Data.xlsx   the supplementary workbook cited in the thesis (40 sheets, see below)

data/physiology/master_data_all_parameters.xlsx         physiology raw data: 22 animals x 13 variables (sheet "Master")
16S/dada2_out/ps_raw.rds                                16S input: ASV table + SILVA 138.2 taxonomy (phyloseq object)
16S/dada2_out/*.csv, *.txt                              DADA2 read tracking, primer-dimer QC, ASV taxonomy, DADA2 settings

scripts/physiology_analysis.R                           1  physiology statistics, Table A1
scripts/16S_amplicon_analysis.R                         2  16S: (DADA2) -> contaminant removal -> diversity, PERMANOVA, ANCOM-BC2
scripts/16S_figure12_and_ancombc2_tables.R              3  Figure 12B/C and the full ANCOM-BC2 tables (+ ALDEx2, interaction model)
scripts/figure_scripts_python/                          4  the printed Figures 3-11 (matplotlib)

results/physiology/tables/                              output of script 1 (CSV + physiology_tables.xlsx)
16S/tables/16S_all_tables.xlsx                          output of script 2 (all 16S tables, n = 17 animals)
16S/tables/16S_ANCOMBC2_full_results.xlsx               output of script 3
16S/tables/cleaning_state.rds                           which ASVs were kept and which animals were analysed (written by script 2, read by script 3)
16S/figures/                                            Figure 12 panels B and C (script 3)
figures/                                                Figures 3-12 exactly as printed in the thesis
```

## How to reproduce the thesis

Everything runs from the **repository root** (in RStudio: open the folder, `Session > Set Working Directory`, then `source()` each script).
No raw reads are needed: script 2 starts from `16S/dada2_out/ps_raw.rds`.

### 1. Install R packages (once)

R 4.4 or newer is recommended (the thesis used R 4.4.3).

```r
install.packages(c("readxl", "dplyr", "tidyr", "purrr", "tibble", "stringr", "readr", "ggplot2", "patchwork",
                   "car", "openxlsx", "vegan", "ggrepel", "ragg", "BiocManager"))
BiocManager::install(c("phyloseq", "ANCOMBC", "ALDEx2"))
# only needed to rerun DADA2 from the raw reads:
BiocManager::install(c("dada2", "Biostrings", "ShortRead"))
```

### 2. Run the scripts in this order

| Step | Command | Produces | Time |
|---|---|---|---|
| 1 | `source("scripts/physiology_analysis.R")` | assumption tests, Welch t-tests (Holm), two-way ANOVA, Table A1 → `results/physiology/tables/` | < 1 min |
| 2 | `source("scripts/16S_amplicon_analysis.R")` | contaminant removal (75.8 % of reads), alpha diversity, PERMANOVA (Tables B1, B2), ANCOM-BC2 → `16S/tables/16S_all_tables.xlsx` | 10–20 min |
| 3 | `source("scripts/16S_figure12_and_ancombc2_tables.R")` | Figure 12B (PCoA with the *Artemia* food library), Figure 12C, full ANCOM-BC2 results with Holm and BH values, ALDEx2, interaction model → `16S/tables/16S_ANCOMBC2_full_results.xlsx`, `16S/figures/` | ~5 min |
| 4 | `python scripts/figure_scripts_python/thesis_physiology_figures.py data/physiology/master_data_all_parameters.xlsx figures` and `python scripts/figure_scripts_python/phylo.py figures` | Figures 3–11 as printed | < 1 min |

Script 3 needs the output of script 2 (`16S/tables/cleaning_state.rds`; a copy from the thesis run is included, so script 3 also runs on its own).
The Python figure scripts need `numpy`, `pandas`, `scipy`, `matplotlib` and `openpyxl`; they use the same statistics as script 1.

### 3. Check that you get the thesis numbers

* **Physiology:** `results/physiology/tables/table_a1.csv` = Table A1 (e.g. symbiont density in *Exaiptasia*: Holm-adjusted p = 0.048).
* **16S:** the console summary of script 2 must report 17 animals, 75.8 % of reads removed, and PERMANOVA host R² = 0.171 (Bray–Curtis) / 0.123 (Aitchison);
  *Exaiptasia* fed vs starved R² = 0.288 / 0.326, p = 0.029; *Waminoa* p = 0.241 / 0.200. The permutations are seeded, so these reproduce exactly.
* **ANCOM-BC2:** no family or genus differs between fed and starved animals in either host after Holm adjustment; between hosts, Lachnospiraceae and
  Acetobacteraceae are higher in *Waminoa* (Holm-adjusted p = 0.015 and 0.021).

## Analysis choices you can change

* **Animals:** the thesis analysis uses 17 animals. Five libraries (AF5, AS1, AS2, WF2, WS4; 870–1,717 reads after cleaning, 73–93 % primer dimer)
  are excluded. Set the environment variable `EXCLUDE_ANIMALS=""` before script 2 to analyse all 20 animals.
* **DADA2 from raw reads:** script 2 runs DADA2 only if `16S/dada2_out/ps_raw.rds` is missing (or `RUN_DADA2=true`). Then put the 50 raw `.fq` files in
  `16S/raw_fastq/` and the SILVA 138.2 files (`silva_nr99_v138.2_toGenus_trainset.fa.gz`, `silva_v138.2_assignSpecies.fa.gz`,
  <https://doi.org/10.5281/zenodo.14169026>) in `16S/reference_db/`. DADA2 takes 2–3 h.

## The supplementary workbook

`supplement/Supplementary_Tables_and_Figure_Data.xlsx` has 40 sheets. The `README` sheet lists every sheet with a link and the thesis section it belongs to.
Tab colours: blue = physiology, green = 16S, orange = figure sheets (the figure as printed + its data).

| Sheets | Content | Thesis |
|---|---|---|
| `Table_A1`, `Phys_Contrasts`, `Phys_ANOVA`, `Phys_Assumptions`, `Phys_Descriptives` | Table A1; fed vs starved tests (Welch, Holm, exact permutation p, Hedges' g); two-way ANOVA; Shapiro–Wilk and Levene; group means | Sections 2.9, 3.1 |
| `Phys_Data`, `Data_dictionary` | every animal and parameter; meaning of each column | Figures 4–11 |
| `Fig04_TotalProtein` … `Fig11_O2perCell` | each figure as printed, with its adjusted p-values and the per-animal values | Figures 4–11 |
| `Fig12_16S`, `Fig12A_Family_bars`, `Fig12B_PCoA_BrayCurtis`, `Fig12C_ANCOMBC2_family` | Figure 12 and the data behind panels A–C (incl. distances to the *Artemia* food library) | Figure 12 |
| `Table_B1_PERMANOVA`, `Table_B1b_Pairwise`, `Table_B2_PERMANOVA_dimer` | PERMANOVA; pairwise tests; model with primer-dimer load as first term | Sections 3.2, 4.5 |
| `16S_Samples`, `16S_ReadTracking`, `16S_Cleaning_Steps`, `16S_NEGCON1`, `16S_Contam_per_sample`, `16S_ASV_status_log`, `16S_ASV_counts_clean` | reads per library through DADA2 and cleaning; primer-dimer load; reads removed per contaminant rule; status of every ASV; decontaminated count table | Section 2.8 |
| `16S_Alpha_diversity` | alpha diversity and tests | Section 3.2 |
| `16S_ANCOMBC2_all_levels`, `16S_ANCOMBC2_summary` | every tested taxon (family and genus, three comparisons) with group means, 95 % CI, Holm and BH values; summary per comparison | Section 3.2, Figure 12C |
| `16S_ANCOMBC2_interaction`, `16S_ALDEx2_within_host` | two checks of the fed vs starved result: host × treatment model and ALDEx2 | Section 4.5 |
| `16S_Top15_*`, `16S_Phylum_percent`, `16S_Family_percent`, `16S_Genus_percent` | taxonomic composition per animal | Figure 12A |

## Notes

* The raw data sheets label *Exaiptasia* as "Aiptasia" (original sample sheets); the scripts relabel it.
* ANCOM-BC2 log fold changes are on the natural-log scale. ANCOM-BC2 reports no confidence limits; the 95 % intervals use the t quantile that reproduces
  each reported p-value, so an interval excludes zero exactly when the unadjusted p < 0.05.
* The in-silico primer coverage check (TestPrime against SILVA: 46 % of bacterial sequences matched) was run online and is not part of the scripts.
* Raw 16S reads (25 libraries, 2 × 250 bp, BMK run CS705-003): ENA/SRA accession **[to be added]**.

## Citation and licence

Kola, D. (2026). *Resolving starvation responses of photosymbiotic holobionts: Code, data and supplementary material* [Data set and software].
Zenodo. <https://doi.org/10.5281/zenodo.23264022> (concept DOI, always resolves to the latest version). See also `CITATION.cff`.

The thesis cites version 1.1.0: <https://doi.org/10.5281/zenodo.23337121>.

Code: MIT (`LICENSE`). Data, tables and figures: CC BY 4.0.
