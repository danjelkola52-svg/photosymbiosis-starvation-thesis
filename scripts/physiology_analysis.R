# =============================================================================
#  PHYSIOLOGY ANALYSIS  -  Exaiptasia and Waminoa, fed vs starved
#  Master's thesis, Danjel Kola (Univ. Oldenburg / HIFMB)
#
#  WHAT THIS SCRIPT DOES (one script, top to bottom)
#    1  reads the physiology workbook (22 animals x 13 measured variables)
#    2  checks the assumptions        Shapiro-Wilk per group, Levene across the four groups
#    3  tests fed vs starved          Welch's t-test within each host, Holm-adjusted across the two hosts,
#                                     plus an exact permutation test and Hedges' g as a robustness check
#    4  tests host x treatment        two-way ANOVA, Type III sums of squares (the interaction term)
#    5  writes Table A1               means +/- SD, change, adjusted p, interaction F and p
#    6  draws Figures 4-11            boxplot + individual animals, Holm-adjusted p above each pair
#    7  saves everything              CSV files, one Excel workbook with a sheet per table, figure files
#
#  HOW TO RUN   from the repository root:   Rscript scripts/physiology_analysis.R
#               or open in RStudio, set the working directory to the repository root and Source it.
#  INPUT        data/physiology/master_data_all_parameters.xlsx   (sheet "Master")
#  OUTPUT       results/physiology/tables/    CSV files + physiology_tables.xlsx
#               results/physiology/figures/   Fig04 ... Fig11 as PDF (vector) and PNG (600 dpi)
#
#  PACKAGES     CRAN: readxl, dplyr, tidyr, purrr, tibble, stringr, readr, ggplot2, patchwork, car, openxlsx
#  Versions used for the thesis: R 4.4.x (see results/physiology/R_session_info.txt after a run).
#
#  NOTE ON THE FIGURES  The figures printed in the thesis (Figures 4-11) were drawn with the Python script
#  scripts/figure_scripts_python/thesis_physiology_figures.py from the same workbook. This script draws the
#  same panels in R with identical statistics (the p-values in figure_pvalues.csv agree), so every number in
#  the thesis can be reproduced in R alone. Figure 3 (phylogeny) is drawn by phylo.py in the same folder.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(stringr); library(readr); library(ggplot2); library(patchwork); library(car); library(openxlsx)
})

# ---- SETTINGS ------------------------------------------------------------------
REPO      <- Sys.getenv("REPO_DIR", if (basename(getwd()) == "scripts") dirname(getwd()) else getwd())
DATA_FILE <- Sys.getenv("PHYSIO_FILE", file.path(REPO, "data", "physiology", "master_data_all_parameters.xlsx"))
OUT_DIR   <- Sys.getenv("OUT_DIR",     file.path(REPO, "results", "physiology"))
P_ADJUST  <- "holm"        # correction across the two hosts, within each parameter (as in the thesis)
ALPHA     <- 0.05
FIG_DPI   <- 600
COL_TREAT <- c(Fed = "#5B8FC7", Starved = "#E08A62")     # thesis colours
HOSTS     <- c("Exaiptasia", "Waminoa")
GROUP_LEVELS <- c("Fed", "Starved")

dir_tab <- file.path(OUT_DIR, "tables");  dir_fig <- file.path(OUT_DIR, "figures")
for (d in c(dir_tab, dir_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

# Type III sums of squares need sum-to-zero contrasts, otherwise the main-effect rows are not meaningful
options(contrasts = c("contr.sum", "contr.poly"))


# =============================================================================
# 1. LOAD
# =============================================================================
say("1. Loading %s", DATA_FILE)
stopifnot(file.exists(DATA_FILE))
raw <- read_excel(DATA_FILE, sheet = "Master")

# The workbook has 16 columns in a fixed order (the headers contain Greek letters and sub/superscripts that
# do not survive every operating system), so they are renamed by position. The order is checked below.
COLS <- c("Organism", "Sample", "Treatment",
          "protein_ug",            # Protein (ug / animal)
          "symbionts_animal",      # Symbiont cells (per animal)
          "symbionts_per_ug",      # Cells / ug protein
          "nh4_uM",                # NH4+ uptake (uM)
          "nh4_per_protein",       # NH4+ / ug protein (uM ug-1)
          "nh4_per_cell",          # NH4+ / cell (uM cell-1)
          "ddO2_mgL_h",            # delta-delta O2 / time (mg L-1 h-1)
          "o2_rate_protein_raw",   # O2 rate / ug protein
          "o2_rate_cell_raw",      # O2 rate / cell
          "fvfm",                  # PAM Fv/Fm
          "o2_ug_h",               # delta O2 absolute (ug O2 h-1)
          "o2_nmol_protein",       # delta O2 (nmol O2 h-1 ug protein-1)
          "o2_nmol_cell")          # delta O2 (nmol O2 h-1 cell-1)
stopifnot(ncol(raw) == length(COLS), tolower(names(raw)[1:3]) == c("organism", "sample", "treatment"))
d <- raw; names(d) <- COLS
d <- d %>%
  mutate(across(-c(Organism, Sample, Treatment), as.numeric),
         Organism  = if_else(Organism == "Aiptasia", "Exaiptasia", Organism),   # the workbook still says Aiptasia
         Organism  = factor(Organism,  levels = HOSTS),
         Treatment = factor(Treatment, levels = GROUP_LEVELS)) %>%
  filter(!is.na(Organism), !is.na(Treatment))
say("   %d animals", nrow(d)); print(table(d$Organism, d$Treatment))
stopifnot(!anyNA(d[, COLS[4:16]]))


# =============================================================================
# 2. PARAMETERS (what is tested, how it is labelled, in which unit it is reported)
#    scale converts the workbook column into the reporting unit used in the thesis
# =============================================================================
params <- tribble(
  ~label,                  ~col,                 ~scale, ~unit,                            ~digits, ~fig, ~stub,
  "Total protein",         "protein_ug",         1,      "ug animal-1",                    1,       4,    "Fig04_TotalProtein",
  "Symbionts per animal",  "symbionts_animal",   1e-3,   "x 10^3 cells animal-1",          1,       5,    "Fig05_SymbiontsPerAnimal",
  "Symbiont density",      "symbionts_per_ug",   1,      "cells ug-1 protein",             1,       6,    "Fig06_SymbiontDensity",
  "Net O2 per animal",     "o2_ug_h",            1,      "ug O2 h-1 animal-1",             2,       7,    "Fig07_NetO2",
  "O2 per protein",        "o2_nmol_protein",    1,      "nmol O2 h-1 ug-1 protein",       2,       8,    "Fig08_O2perProtein",
  "Fv/Fm",                 "fvfm",               1,      "dimensionless",                  3,       9,    "Fig09_FvFm",
  "NH4 per protein",       "nh4_per_protein",    1,      "uM ug-1 protein",                3,       10,   "Fig10_NH4perProtein",
  "O2 per cell",           "o2_nmol_cell",       1e3,    "pmol O2 h-1 cell-1",             1,       11,   "Fig11_O2perCell"
)
# y-axis labels (plotmath) for the figures
YLAB <- c(
  "Total protein"        = "'Total protein ('*mu*'g animal'^-1*')'",
  "Symbionts per animal" = "'Symbionts per animal ('*'\u00d710'^3*')'",
  "Symbiont density"     = "atop('Symbiont density', '(cells '*mu*'g'^-1*' protein)')",
  "Net O2 per animal"    = "atop('Net O'[2]*' production', '('*mu*'g O'[2]*' h'^-1*' animal'^-1*')')",
  "O2 per protein"       = "atop('Net O'[2]*' per protein', '(nmol O'[2]*' h'^-1*' '*mu*'g'^-1*')')",
  "Fv/Fm"                = "'Maximum quantum yield (F'[v]*'/F'[m]*')'",
  "NH4 per protein"      = "'NH'[4]^'+'*' uptake ('*mu*'M '*mu*'g'^-1*' protein)'",
  "O2 per cell"          = "atop('Net O'[2]*' per symbiont cell', '(pmol O'[2]*' h'^-1*' cell'^-1*')')")

# long table: one row per animal x parameter, already in the reporting unit
long <- map_dfr(seq_len(nrow(params)), function(i)
  tibble(Parameter = params$label[i], Organism = d$Organism, Treatment = d$Treatment, Sample = d$Sample,
         value = d[[params$col[i]]] * params$scale[i]))
long$Parameter <- factor(long$Parameter, levels = params$label)


# =============================================================================
# 3. DESCRIPTIVES AND ASSUMPTION CHECKS
# =============================================================================
say("2. Descriptives and assumptions")
descriptives <- long %>%
  group_by(Parameter, Organism, Treatment) %>%
  summarise(n = n(), Mean = mean(value), SD = sd(value), Median = median(value),
            Min = min(value), Max = max(value), .groups = "drop")

levene <- long %>% group_by(Parameter) %>%
  group_modify(~ tibble(Levene_p = leveneTest(value ~ interaction(Organism, Treatment), data = .x)$`Pr(>F)`[1])) %>%
  ungroup()
assumptions <- long %>%
  group_by(Parameter, Organism, Treatment) %>%
  summarise(n = n(), Shapiro_p = shapiro.test(value)$p.value, .groups = "drop") %>%
  left_join(levene, by = "Parameter") %>%
  mutate(Normality_ok = Shapiro_p >= ALPHA, Equal_variance_ok = Levene_p >= ALPHA)
say("   Shapiro-Wilk p < %.2f in %d of %d groups; Levene p < %.2f for: %s", ALPHA,
    sum(!assumptions$Normality_ok), nrow(assumptions), ALPHA,
    paste(unique(as.character(assumptions$Parameter[!assumptions$Equal_variance_ok])), collapse = ", "))


# =============================================================================
# 4. FED vs STARVED WITHIN EACH HOST   (Welch + Holm; exact permutation and Hedges' g as checks)
#    Welch's t-test does not assume equal variances. Holm is applied across the two hosts of each parameter.
# =============================================================================
say("3. Within-host contrasts (Welch, %s across the two hosts)", P_ADJUST)
exact_perm_p <- function(f, s) {           # exact two-sided permutation test of the difference in means
  all_v <- c(f, s); n1 <- length(f); obs <- abs(mean(f) - mean(s))
  idx <- combn(length(all_v), n1)
  mean(apply(idx, 2, function(k) abs(mean(all_v[k]) - mean(all_v[-k]))) >= obs - 1e-12)
}
hedges_g <- function(f, s) {               # bias-corrected standardised mean difference (pooled SD)
  n1 <- length(f); n2 <- length(s)
  sp <- sqrt(((n1 - 1) * var(f) + (n2 - 1) * var(s)) / (n1 + n2 - 2))
  (mean(f) - mean(s)) / sp * (1 - 3 / (4 * (n1 + n2 - 2) - 1))
}
contrasts <- map_dfr(levels(long$Parameter), function(par) {
  map_dfr(HOSTS, function(h) {
    x <- long %>% filter(Parameter == par, Organism == h)
    f <- x$value[x$Treatment == "Fed"]; s <- x$value[x$Treatment == "Starved"]
    tt <- t.test(f, s)                                              # Welch is R's default
    tibble(Parameter = par, Organism = h, n_Fed = length(f), n_Starved = length(s),
           Fed_mean = mean(f), Fed_SD = sd(f), Starved_mean = mean(s), Starved_SD = sd(s),
           Change_pct = 100 * (mean(s) - mean(f)) / mean(f),        # starved relative to fed
           t = unname(tt$statistic), df = unname(tt$parameter), p_Welch = tt$p.value,
           p_permutation = exact_perm_p(f, s), Hedges_g = hedges_g(f, s))
  }) %>% mutate(p_adj = p.adjust(p_Welch, method = P_ADJUST))
}) %>% mutate(Parameter = factor(Parameter, levels = params$label)) %>% arrange(Parameter, Organism)


# =============================================================================
# 5. HOST x TREATMENT   (two-way ANOVA, Type III)
#    The interaction asks whether the two hosts respond differently to starvation. It is kept although
#    some groups depart from normality, because no non-parametric test addresses an interaction directly.
# =============================================================================
say("4. Two-way ANOVA (Type III)")
anova_tab <- map_dfr(levels(long$Parameter), function(par) {
  x <- long %>% filter(Parameter == par)
  m <- lm(value ~ Organism * Treatment, data = x)
  a <- car::Anova(m, type = "III")
  tibble(Parameter = par, Term = rownames(a), Df = a$Df, F = a$`F value`, p = a$`Pr(>F)`, Df_residual = df.residual(m)) %>%
    filter(Term %in% c("Organism", "Treatment", "Organism:Treatment"))
}) %>% mutate(Parameter = factor(Parameter, levels = params$label))


# =============================================================================
# 6. TABLE A1 (as printed in Appendix A / the supplement)
# =============================================================================
fmt_num <- function(x, par) {          # number of decimals follows the reporting unit (params$digits)
  dig <- params$digits[match(as.character(par), params$label)]
  mapply(function(v, k) formatC(v, format = "f", digits = k), x, dig)
}
fmt_change <- function(x) {            # rounded % change with an explicit sign; "0" instead of "-0"
  r <- round(x); ifelse(r == 0, "0", sprintf("%+d", as.integer(r)))
}
fmt_p3 <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
int_row <- anova_tab %>% filter(Term == "Organism:Treatment") %>%
  transmute(Parameter, Interaction = sprintf("F = %.2f, p %s", F, ifelse(p < 0.001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3)))))
table_A1 <- contrasts %>%
  left_join(int_row, by = "Parameter") %>%
  group_by(Parameter) %>%
  mutate(Interaction = if_else(row_number() == 1, Interaction, "")) %>% ungroup() %>%
  transmute(Parameter = as.character(Parameter), Unit = params$unit[match(Parameter, params$label)], Host = as.character(Organism),
            Fed = sprintf("%s ± %s", fmt_num(Fed_mean, Parameter), fmt_num(Fed_SD, Parameter)),
            Starved = sprintf("%s ± %s", fmt_num(Starved_mean, Parameter), fmt_num(Starved_SD, Parameter)),
            `Change (%)` = fmt_change(Change_pct), `p (fed vs starved, Holm)` = fmt_p3(p_adj),
            `Organism x treatment` = Interaction)


# =============================================================================
# 7. FIGURES 4-11
# =============================================================================
say("5. Figures")
theme_thesis <- theme_classic(base_size = 10) +
  theme(axis.line = element_line(linewidth = .4), axis.ticks = element_line(linewidth = .4),
        axis.text = element_text(colour = "black"), axis.title = element_text(colour = "black", size = 10),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3),
        plot.tag = element_text(face = "bold", size = 13), legend.position = "none")
fmt_p_fig <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

make_panel <- function(par, host, ylim, p_adj) {
  x <- long %>% filter(Parameter == par, Organism == host)
  n <- x %>% count(Treatment)
  lab <- setNames(sprintf("%s\n(n = %d)", n$Treatment, n$n), as.character(n$Treatment))
  span <- diff(ylim); yb <- ylim[2] - .10 * span
  ggplot(x, aes(Treatment, value)) +
    geom_boxplot(aes(fill = Treatment), width = .55, alpha = .85, outlier.shape = NA, colour = "black", linewidth = .4) +
    geom_jitter(aes(fill = Treatment), width = .10, height = 0, size = 2, shape = 21, colour = "black", stroke = .35) +
    scale_fill_manual(values = COL_TREAT) +
    scale_x_discrete(labels = lab) +
    scale_y_continuous(labels = scales::label_comma()) +
    coord_cartesian(ylim = ylim) +
    annotate("segment", x = 1, xend = 2, y = yb, yend = yb, linewidth = .4) +
    annotate("segment", x = 1, xend = 1, y = yb, yend = yb - .022 * span, linewidth = .4) +
    annotate("segment", x = 2, xend = 2, y = yb, yend = yb - .022 * span, linewidth = .4) +
    annotate("text", x = 1.5, y = yb + .03 * span, label = fmt_p_fig(p_adj), size = 3.1,
             fontface = if (p_adj < ALPHA) "bold" else "plain") +
    labs(x = NULL, y = NULL) + theme_thesis
}
for (i in seq_len(nrow(params))) {
  par  <- params$label[i]
  vals <- long$value[long$Parameter == par]
  rng  <- range(vals); span <- diff(rng)
  ylim <- c(rng[1] - .08 * span, rng[2] + .26 * span)             # same y-axis in both panels, room for the bracket
  pa <- make_panel(par, "Exaiptasia", ylim, contrasts$p_adj[contrasts$Parameter == par & contrasts$Organism == "Exaiptasia"]) +
    labs(y = parse(text = YLAB[[par]])[[1]])
  pb <- make_panel(par, "Waminoa", ylim, contrasts$p_adj[contrasts$Parameter == par & contrasts$Organism == "Waminoa"]) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  fig <- (pa | pb) + plot_annotation(tag_levels = "A")
  ggsave(file.path(dir_fig, paste0(params$stub[i], ".pdf")), fig, width = 6.3, height = 3.4, device = cairo_pdf)
  ggsave(file.path(dir_fig, paste0(params$stub[i], ".png")), fig, width = 6.3, height = 3.4, dpi = FIG_DPI)
  say("   Figure %d  %-22s Exaiptasia %-10s Waminoa %-10s", params$fig[i], par,
      fmt_p_fig(contrasts$p_adj[contrasts$Parameter == par & contrasts$Organism == "Exaiptasia"]),
      fmt_p_fig(contrasts$p_adj[contrasts$Parameter == par & contrasts$Organism == "Waminoa"]))
}


# =============================================================================
# 8. SAVE TABLES  (CSV + one Excel workbook, one sheet per table)
# =============================================================================
say("6. Writing tables")
figure_pvalues <- contrasts %>% transmute(Figure = params$fig[match(Parameter, params$label)], Parameter = as.character(Parameter),
                                          Host = as.character(Organism), n_Fed, n_Starved, p_Welch, p_Holm = p_adj) %>% arrange(Figure, Host)
data_dictionary <- tibble(
  Column = COLS,
  Meaning = c("Host species (the workbook says 'Aiptasia'; relabelled Exaiptasia in the analysis)", "Animal ID (A = Exaiptasia, W = Waminoa; F = fed, S = starved)", "Feeding treatment",
              "Soluble host protein per animal (ug)", "Symbiont cells per animal", "Symbiont cells per ug host protein",
              "Net NH4+ uptake over the 3 h incubation (uM)", "NH4+ uptake per ug host protein (uM ug-1)", "NH4+ uptake per symbiont cell (uM cell-1)",
              "Blank-corrected O2 change per time (mg L-1 h-1)", "O2 rate per ug protein (workbook value)", "O2 rate per cell (workbook value)",
              "Maximum quantum yield of PSII after dark adaptation", "Net O2 production per animal (ug O2 h-1)",
              "Net O2 production per ug host protein (nmol O2 h-1 ug-1)", "Net O2 production per symbiont cell (nmol O2 h-1 cell-1; x 1000 = pmol in the thesis)"))
tables <- list(README = tibble(Sheet = c("Table_A1", "Contrasts", "ANOVA", "Descriptives", "Assumptions", "Figure_pvalues", "Data_long", "Data_dictionary"),
                               Content = c("Table A1 as printed in the thesis (means +/- SD, change, Holm-adjusted p, interaction)",
                                           "Welch t-test, Holm-adjusted p, exact permutation p and Hedges' g for every parameter and host",
                                           "Two-way ANOVA, Type III: host, treatment and host x treatment",
                                           "n, mean, SD, median, min, max per group", "Shapiro-Wilk per group and Levene's test per parameter",
                                           "p-values printed on Figures 4-11", "Every animal x parameter value in the reporting unit", "Meaning of each workbook column")),
               Table_A1 = table_A1, Contrasts = contrasts, ANOVA = anova_tab, Descriptives = descriptives, Assumptions = assumptions,
               Figure_pvalues = figure_pvalues, Data_long = long, Data_dictionary = data_dictionary)
for (nm in setdiff(names(tables), "README")) write_csv(tables[[nm]], file.path(dir_tab, paste0(tolower(nm), ".csv")))
wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#E8EEF6", border = "bottom")
for (nm in names(tables)) {
  addWorksheet(wb, nm); writeData(wb, nm, as.data.frame(tables[[nm]]), headerStyle = hdr)
  setColWidths(wb, nm, cols = seq_len(ncol(tables[[nm]])), widths = "auto"); freezePane(wb, nm, firstRow = TRUE)
}
saveWorkbook(wb, file.path(dir_tab, "physiology_tables.xlsx"), overwrite = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "R_session_info.txt"))

cat("\n=== TABLE A1 ===\n"); print(as.data.frame(table_A1), row.names = FALSE, right = FALSE)
say("DONE. Tables: %s | Figures: %s", dir_tab, dir_fig)
