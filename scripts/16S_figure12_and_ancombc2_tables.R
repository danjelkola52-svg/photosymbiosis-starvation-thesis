# =============================================================================
#  16S FIGURE 12 (panels B and C) AND FULL ANCOM-BC2 TABLES
#  Master's thesis, Danjel Kola (Univ. Oldenburg / HIFMB)
#
#  RUN AFTER  scripts/16S_amplicon_analysis.R   (it writes 16S/tables/cleaning_state.rds, which this script reads,
#                                                so that exactly the same ASVs and the same 17 animals are used)
#  HOW TO RUN from the repository root:   Rscript scripts/16S_figure12_and_ancombc2_tables.R
#
#  WHAT IT DOES
#    1  rebuilds the decontaminated count table of the 17 analysed animals and the Artemia food library (ARTSEP)
#    2  Figure 12B   Bray-Curtis PCoA fitted jointly to the 17 animals and the Artemia library (axes 22.6 % / 15.9 %);
#                    the food library is shown for comparison only and is not part of any test
#    3  ANCOM-BC2    fed vs starved within Exaiptasia, within Waminoa, and Waminoa vs Exaiptasia (adjusted for treatment),
#                    at family and genus level, same settings as 16S_amplicon_analysis.R
#                    (prv_cut 0.30, struc_zero TRUE, neg_lb FALSE, pseudo-count sensitivity analysis, Holm adjustment).
#                    Benjamini-Hochberg values are added for comparison only.
#    4  Figure 12C   natural-log fold changes of the six classified families with the smallest p per comparison.
#                    ANCOM-BC2 gives no confidence limits and its p-values come from a t reference distribution,
#                    so each 95 % interval uses the t quantile that reproduces the reported p (it excludes 0 exactly
#                    when the unadjusted p < 0.05).
#    5  two checks of the within-host result: an organism x treatment interaction model (ANCOM-BC2, all 17 animals)
#                    and ALDEx2 (CLR, Welch t-test, BH) as an independent second method
#    6  writes 16S/tables/16S_ANCOMBC2_full_results.xlsx and the panel figures to 16S/figures/
#
#  PACKAGES  CRAN: vegan, ggplot2, openxlsx, patchwork, ragg      Bioconductor: phyloseq, ANCOMBC, ALDEx2
#  Run time  about 5 minutes.
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(ANCOMBC); library(ALDEx2); library(vegan)
  library(ggplot2); library(patchwork); library(openxlsx)
})

REPO    <- Sys.getenv("REPO_DIR", if (basename(getwd()) == "scripts") dirname(getwd()) else getwd())
PROJ    <- file.path(REPO, "16S")
PS_RAW  <- file.path(PROJ, "dada2_out", "ps_raw.rds")
STATE   <- file.path(PROJ, "tables", "cleaning_state.rds")
OUT_TAB <- file.path(PROJ, "tables");  OUT_FIG <- file.path(PROJ, "figures")
ARTEMIA <- "ARTSEP"
PREV    <- 0.30; ALPHA <- 0.05; SEED <- 123
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
if (!file.exists(STATE)) stop("Run scripts/16S_amplicon_analysis.R first: ", STATE, " is missing.")

# ---- 1. decontaminated counts of the 17 animals (+ Artemia) ----------------------------------------------
st <- readRDS(STATE)
ps <- readRDS(PS_RAW)
if (taxa_are_rows(ps)) ps <- phyloseq(t(otu_table(ps)), tax_table(ps), sample_data(ps))
stopifnot(ntaxa(ps) == length(st$status))
keep <- st$status == levels(st$status)[length(levels(st$status))]           # last level = "kept"
ps_all <- prune_taxa(taxa_names(ps)[keep], prune_samples(c(st$ani, ARTEMIA), ps))
host  <- ifelse(grepl("^A[FS]", sample_names(ps_all)), "Exaiptasia", ifelse(grepl("^W", sample_names(ps_all)), "Waminoa", "Artemia"))
treat <- ifelse(grepl("^.F", sample_names(ps_all)), "Fed", ifelse(grepl("^.S", sample_names(ps_all)), "Starved", "Food"))
sample_data(ps_all) <- sample_data(data.frame(Organism = host, Treatment = treat, row.names = sample_names(ps_all)))
ps17 <- prune_samples(st$ani, ps_all); ps17 <- prune_taxa(taxa_sums(ps17) > 0, ps17)
sd17 <- data.frame(sample_data(ps17))
sd17$Treatment <- factor(sd17$Treatment, levels = c("Fed", "Starved"))
sd17$Organism  <- factor(sd17$Organism,  levels = c("Exaiptasia", "Waminoa"))
sample_data(ps17) <- sample_data(sd17)
say("animals: %d (%s) | ASVs: %d", nsamples(ps17), paste(table(paste(sd17$Organism, sd17$Treatment)), collapse = "/"), ntaxa(ps17))

# ---- 2. Figure 12B: joint PCoA of animals + Artemia ----------------------------------------------------
cnt <- as(otu_table(ps_all), "matrix"); cnt <- cnt[, colSums(cnt) > 0]
bc  <- vegdist(cnt / rowSums(cnt), "bray")
pc  <- cmdscale(bc, k = 2, eig = TRUE); ev <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
# the sign of a PCoA axis is arbitrary; orient PCo2 so that the Artemia library is at the top, as in the printed figure
if (pc$points[ARTEMIA, 2] < 0) pc$points[, 2] <- -pc$points[, 2]
md  <- data.frame(Sample = rownames(cnt), Organism = host[match(rownames(cnt), sample_names(ps_all))],
                  Treatment = treat[match(rownames(cnt), sample_names(ps_all))], PCo1 = pc$points[, 1], PCo2 = pc$points[, 2])
md$Group <- ifelse(md$Organism == "Artemia", "Artemia (food)", paste0(md$Organism, ", ", tolower(md$Treatment)))
say("PCoA axes: %.1f %% / %.1f %%", ev[1], ev[2])
dm <- as.matrix(bc); an <- st$ani
art_dist <- data.frame(Sample = an, Group = md$Group[match(an, md$Sample)], BrayCurtis_to_Artemia = round(dm[an, ARTEMIA], 3))
art_summary <- data.frame(Measure = c("Median Bray-Curtis distance, animals to Artemia", "Median Bray-Curtis distance among animals"),
                          Value = round(c(median(dm[an, ARTEMIA]), median(dm[an, an][upper.tri(dm[an, an])])), 2))
hull <- do.call(rbind, lapply(split(md[md$Organism != "Artemia", ], md$Group[md$Organism != "Artemia"]), function(x) x[chull(x$PCo1, x$PCo2), ]))
COL <- c("Exaiptasia, fed" = "#5B8FC7", "Exaiptasia, starved" = "#E8826A", "Waminoa, fed" = "#5B8FC7", "Waminoa, starved" = "#E8826A", "Artemia (food)" = "#2BA37E")
SHP <- c("Exaiptasia, fed" = 21, "Exaiptasia, starved" = 21, "Waminoa, fed" = 22, "Waminoa, starved" = 22, "Artemia (food)" = 23)
LTY <- c(Exaiptasia = "solid", Waminoa = "dashed")
pB <- ggplot(md, aes(PCo1, PCo2)) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
  geom_polygon(data = hull, aes(group = Group, fill = Group, colour = Group, linetype = Organism), alpha = 0.15, linewidth = 0.45, show.legend = FALSE) +
  geom_point(aes(fill = Group, shape = Group), size = 2.8, colour = "grey15", stroke = 0.4) +
  geom_text(aes(label = ifelse(Organism == "Artemia", "Artemia food", Sample), fontface = ifelse(Organism == "Artemia", "italic", "plain")),
            size = 2.3, colour = "grey30", nudge_y = 0.03) +
  scale_fill_manual(values = COL, name = NULL) + scale_colour_manual(values = COL, guide = "none") +
  scale_shape_manual(values = SHP, name = NULL) + scale_linetype_manual(values = LTY, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.05))) +
  labs(x = sprintf("PCo1 (%.1f %%)", ev[1]), y = sprintf("PCo2 (%.1f %%)", ev[2])) +
  guides(fill = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  theme_classic(base_size = 9) + theme(legend.position = "bottom", legend.text = element_text(face = "italic"))
ggsave(file.path(OUT_FIG, "Figure12B_PCoA_BrayCurtis_with_Artemia.pdf"), pB, width = 4.2, height = 4.6, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "Figure12B_PCoA_BrayCurtis_with_Artemia.png"), pB, width = 4.2, height = 4.6, dpi = 600, device = ragg::agg_png)

# ---- 3. ANCOM-BC2 -------------------------------------------------------------------------------------
glom <- function(p, rank) {                      # collapse to a rank; short IDs (long names crash ancombc2)
  p <- prune_taxa(taxa_sums(p) > 0, p); g <- tax_glom(p, rank, NArm = FALSE)
  tax <- as.data.frame(as(tax_table(g), "matrix"), stringsAsFactors = FALSE)
  lab <- tax[[rank]]; par <- if (rank == "Genus") tax$Family else tax$Order
  na  <- is.na(lab) | lab == ""; lab[na] <- paste("unclassified", ifelse(is.na(par[na]), "Bacteria", par[na]))
  lab <- sub("([a-z])1$", "\\1", lab)          # SILVA suffix, e.g. "Alcanivoracaceae1"
  ids <- sprintf("T%03d", seq_len(ntaxa(g))); taxa_names(g) <- ids; list(g = g, key = setNames(lab, ids)) }
fit_ancombc2 <- function(g, ff, grp, struc_zero = TRUE) {   # a taxon found in one group only stops ancombc2: set aside, refit
  drop <- character(0)
  for (i in 1:10) {
    gt <- prune_taxa(setdiff(taxa_names(g), drop), g); set.seed(SEED)
    f <- tryCatch(ancombc2(data = gt, fix_formula = ff, group = grp, p_adj_method = "holm", prv_cut = PREV, lib_cut = 0,
                           struc_zero = struc_zero, neg_lb = FALSE, alpha = ALPHA, pseudo_sens = TRUE, n_cl = 1, verbose = FALSE),
                  error = function(e) e)
    if (!inherits(f, "error")) return(f$res)
    bad <- intersect(unlist(strsplit(conditionMessage(f), "[^A-Za-z0-9_]+")), taxa_names(gt))
    if (!length(bad)) stop(conditionMessage(f)); drop <- union(drop, bad) } }
t_ci <- function(lfc, se, p) {                   # t quantile that reproduces the reported two-sided p
  W <- abs(lfc / se)
  crit <- mapply(function(w, pp) {
    if (!is.finite(w) || is.na(pp) || pp <= 2 * pnorm(-w)) return(qnorm(0.975))
    d <- tryCatch(uniroot(function(d) 2 * pt(-w, d) - pp, c(0.02, 1e7))$root, error = function(e) 1e7)
    qt(0.975, d) }, W, p)
  cbind(lfc - crit * se, lfc + crit * se) }
grp_lab <- paste(sd17$Organism, tolower(sd17$Treatment))
jobs <- list(
  list(name = "Exaiptasia: fed vs starved", s = rownames(sd17)[sd17$Organism == "Exaiptasia"], ff = "Treatment", grp = "Treatment", coef = "TreatmentStarved",
       cols = c("Exaiptasia fed", "Exaiptasia starved"), sheet = "Exaip_fed_vs_starved"),
  list(name = "Waminoa: fed vs starved", s = rownames(sd17)[sd17$Organism == "Waminoa"], ff = "Treatment", grp = "Treatment", coef = "TreatmentStarved",
       cols = c("Waminoa fed", "Waminoa starved"), sheet = "Wam_fed_vs_starved"),
  list(name = "Waminoa vs Exaiptasia (adjusted for treatment)", s = rownames(sd17), ff = "Organism + Treatment", grp = "Organism", coef = "OrganismWaminoa",
       cols = sort(unique(grp_lab)), sheet = "Exaip_vs_Waminoa"))
res <- list()
for (j in jobs) for (rk in c("Family", "Genus")) {
  say("ANCOM-BC2 %s, %s", j$name, rk)
  p  <- prune_samples(j$s, ps17); gg <- glom(p, rk); r <- fit_ancombc2(gg$g, j$ff, j$grp); s <- j$coef
  rel <- as(otu_table(gg$g), "matrix"); if (taxa_are_rows(gg$g)) rel <- t(rel); rel <- 100 * rel / rowSums(rel)
  gl  <- grp_lab[match(rownames(rel), rownames(sd17))]
  mean_pct <- sapply(j$cols, function(cc) colMeans(rel[gl == cc, r$taxon, drop = FALSE]))
  present  <- sapply(j$cols, function(cc) sprintf("%d of %d", colSums(rel[gl == cc, r$taxon, drop = FALSE] > 0), sum(gl == cc)))
  if (is.null(dim(mean_pct))) { mean_pct <- t(mean_pct); present <- t(present) }
  lfc <- r[[paste0("lfc_", s)]]; se <- r[[paste0("se_", s)]]; pv <- r[[paste0("p_", s)]]; ci <- t_ci(lfc, se, pv)
  ss  <- r[[paste0("passed_ss_", s)]]; qh <- p.adjust(pv, "holm")
  d <- data.frame(Comparison = j$name, Level = rk, Taxon = unname(gg$key[r$taxon]),
                  setNames(as.data.frame(round(mean_pct, 3)), paste0("Mean_pct_", gsub(" ", "_", j$cols))),
                  setNames(as.data.frame(present), paste0("Present_", gsub(" ", "_", j$cols))),
                  LFC_ln = round(lfc, 3), Fold_change = round(exp(lfc), 2), SE = round(se, 3),
                  CI95_low = round(ci[, 1], 3), CI95_high = round(ci[, 2], 3),
                  p_value = signif(pv, 3), p_adj_Holm = signif(qh, 3), q_BH = signif(p.adjust(pv, "BH"), 3),
                  Passed_pseudocount_sensitivity = ss,
                  Verdict = ifelse(qh < ALPHA & ss, "different (robust)", ifelse(qh < ALPHA, "significant, not robust",
                            ifelse(pv < ALPHA, "unadjusted p < 0.05 only", "no difference"))),
                  check.names = FALSE, row.names = NULL)
  res[[paste(j$sheet, rk, sep = "_")]] <- d[order(d$p_value), ] }
summ <- do.call(rbind, lapply(res, function(d) data.frame(Comparison = d$Comparison[1], Level = d$Level[1], Taxa_tested = nrow(d),
  Unadjusted_p_below_0.05 = sum(d$p_value < ALPHA), Expected_by_chance = round(ALPHA * nrow(d), 1), Smallest_p = min(d$p_value),
  Holm_below_0.05 = sum(d$p_adj_Holm < ALPHA), Holm_below_0.05_and_robust = sum(d$Verdict == "different (robust)"),
  BH_below_0.05 = sum(d$q_BH < ALPHA), Smallest_BH = min(d$q_BH))))
rownames(summ) <- NULL; print(summ)

# ---- 4. Figure 12C ------------------------------------------------------------------------------------
keep_cols <- c("Comparison", "Taxon", "LFC_ln", "CI95_low", "CI95_high", "p_value", "p_adj_Holm", "Passed_pseudocount_sensitivity")
fc <- do.call(rbind, lapply(res[grepl("_Family$", names(res))], function(d) head(d[!grepl("^unclassified", d$Taxon), keep_cols], 6)))
fc$Comparison <- factor(fc$Comparison, levels = sapply(jobs, `[[`, "name"),
                        labels = c("Exaiptasia: starved vs fed", "Waminoa: starved vs fed", "Waminoa vs Exaiptasia"))
fc$Class <- factor(ifelse(fc$p_adj_Holm < ALPHA & fc$Passed_pseudocount_sensitivity, "Holm-adjusted p < 0.05, robust to pseudo-count",
                   ifelse(fc$p_adj_Holm < ALPHA, "Holm-adjusted p < 0.05, not robust to pseudo-count", "Unadjusted p < 0.05, Holm-adjusted p ≥ 0.05")),
                   levels = c("Holm-adjusted p < 0.05, robust to pseudo-count", "Holm-adjusted p < 0.05, not robust to pseudo-count",
                              "Unadjusted p < 0.05, Holm-adjusted p ≥ 0.05"))
fc$key <- paste(fc$Comparison, fc$Taxon); fc$key <- factor(fc$key, levels = rev(fc$key))
fc$lab <- ifelse(fc$p_adj_Holm < 0.1, sprintf("p[adj] == %.3f", fc$p_adj_Holm), sprintf("p[adj] == %.2f", fc$p_adj_Holm))
pC <- ggplot(fc, aes(LFC_ln, key)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = CI95_low, xmax = CI95_high, colour = Class), height = 0, linewidth = 0.5) +
  geom_point(aes(colour = Class, fill = Class), shape = 21, size = 2.3, stroke = 0.7) +
  geom_text(aes(x = 7.4, label = lab), parse = TRUE, hjust = 1, size = 2.2, colour = "grey30") +
  scale_colour_manual(values = c("#B5182B", "#B5182B", "grey35"), name = NULL, drop = FALSE) +
  scale_fill_manual(values = c("#B5182B", "white", "grey55"), name = NULL, drop = FALSE) +
  scale_y_discrete(labels = function(x) sub("^.*? (?=[A-Z][a-z]+aceae$|unclassified)", "", x, perl = TRUE)) +
  scale_x_continuous(limits = c(-7, 7.5), breaks = seq(-6, 6, 2)) +
  facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
  labs(x = "Natural-log fold change (95 % CI)", y = NULL) +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "italic", size = 9),
                                       legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(OUT_FIG, "Figure12C_ANCOMBC2_family.pdf"), pC, width = 4.8, height = 6.2, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "Figure12C_ANCOMBC2_family.png"), pC, width = 4.8, height = 6.2, dpi = 600, device = ragg::agg_png)

# ---- 5. checks: interaction model and ALDEx2 ----------------------------------------------------------
say("ANCOM-BC2 interaction model (family)")
gf <- glom(ps17, "Family"); r <- fit_ancombc2(gf$g, "Organism * Treatment", NULL, struc_zero = FALSE)
s  <- "OrganismWaminoa:TreatmentStarved"
inter <- data.frame(Taxon = unname(gf$key[r$taxon]), LFC_interaction_ln = round(r[[paste0("lfc_", s)]], 3), SE = round(r[[paste0("se_", s)]], 3),
                    p_value = signif(r[[paste0("p_", s)]], 3), p_adj_Holm = signif(p.adjust(r[[paste0("p_", s)]], "holm"), 3),
                    q_BH = signif(p.adjust(r[[paste0("p_", s)]], "BH"), 3), Passed_pseudocount_sensitivity = r[[paste0("passed_ss_", s)]],
                    LFC_starved_in_Exaiptasia = round(r$lfc_TreatmentStarved, 3), p_starved_in_Exaiptasia = signif(r$p_TreatmentStarved, 3), row.names = NULL)
inter <- inter[order(inter$p_value), ]
say("ALDEx2 within hosts")
ald <- do.call(rbind, lapply(c("Exaiptasia", "Waminoa"), function(h) do.call(rbind, lapply(c("Family", "Genus"), function(rk) {
  p <- prune_samples(rownames(sd17)[sd17$Organism == h], ps17); g <- glom(p, rk)
  m <- as(otu_table(g$g), "matrix"); if (!taxa_are_rows(g$g)) m <- t(m); rownames(m) <- make.unique(unname(g$key[rownames(m)]))
  m <- m[rowSums(m > 0) >= ceiling(PREV * ncol(m)), ]
  set.seed(1); x <- aldex(m, as.character(sample_data(g$g)$Treatment), mc.samples = 256, test = "t", effect = TRUE, denom = "all", verbose = FALSE)
  x <- x[order(x$we.ep), ]
  data.frame(Host = h, Level = rk, Taxon = rownames(x), CLR_difference_starved_minus_fed = signif(x$diff.btw, 3),
             Effect_size = signif(x$effect, 3), Welch_p = signif(x$we.ep, 3), Welch_BH_q = signif(x$we.eBH, 3), row.names = NULL) }))))

# ---- 6. Excel workbook ----------------------------------------------------------------------------------
readme <- data.frame(Item = c("Content", "Data", "ANCOM-BC2", "Adjustment", "Robust", "Confidence interval", "Log scale", "Sign",
                              "Small samples", "Interaction", "ALDEx2", "Figure 12", "Script"),
  Description = c(
  "Complete ANCOM-BC2 results for every tested taxon in the three comparisons of Figure 12C, at family and genus level, with two checks.",
  "17 animals after contaminant removal: Exaiptasia fed 4 (AF1-4), starved 3 (AS3-5); Waminoa fed 5 (WF1,3,4,5,6), starved 5 (WS1,2,3,5,6). Same ASVs as 16S_all_tables.xlsx.",
  "ANCOMBC 2.8, ancombc2(prv_cut = 0.30, lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, pseudo_sens = TRUE, alpha = 0.05). Taxa present in < 30 % of samples were not tested.",
  "p_adj_Holm: Holm adjustment across all taxa tested in a comparison and level (reported in the thesis). q_BH: Benjamini-Hochberg, given for comparison only.",
  "Passed_pseudocount_sensitivity: ANCOM-BC2's check that the result holds when the pseudo-count added to zeros changes. A difference is accepted only if Holm-adjusted p < 0.05 and this is TRUE.",
  "ANCOM-BC2 reports no confidence limits; its p-values use a t reference distribution. CI95 uses the t quantile that reproduces the reported p, so it excludes 0 exactly when p < 0.05.",
  "Fold changes are natural log (ln); Fold_change = exp(LFC_ln).",
  "Positive LFC: higher in starved animals (within-host sheets) or in Waminoa (host sheets).",
  "ANCOM-BC2 warns that variance estimates are unstable below 5 animals per group (Exaiptasia: 4 vs 3). An exact rank test with 4 vs 3 animals cannot reach p < 0.057.",
  "Exploratory: Organism * Treatment model on all 17 animals (family level, struc_zero = FALSE). Hits that fail the pseudo-count check do not meet the thesis criterion.",
  "ALDEx2 (CLR, 256 Monte-Carlo instances, Welch t-test, Benjamini-Hochberg) within each host, as an independent second method.",
  "Panel B: PCoA (Bray-Curtis) fitted jointly to the 17 animals and the Artemia library. Panel C: six classified families with the smallest p per comparison.",
  "scripts/16S_figure12_and_ancombc2_tables.R"))
wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", fgFill = "#DDE6F0", wrapText = TRUE, valign = "top")
add <- function(nm, d, w = NULL) { addWorksheet(wb, nm); writeData(wb, nm, d, headerStyle = hs); freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, seq_len(ncol(d)), if (is.null(w)) "auto" else w) }
add("README", readme, c(22, 150)); add("Summary", summ)
for (nm in names(res)) add(nm, res[[nm]])
add("Interaction_Family", inter); add("ALDEx2_within_host", ald)
add("Fig12B_PCoA_coordinates", data.frame(md[, c("Sample", "Group", "PCo1", "PCo2")], PCo1_pct = round(ev[1], 1), PCo2_pct = round(ev[2], 1)))
add("Fig12B_Artemia_distance", rbind(data.frame(Sample = art_summary$Measure, Group = "", BrayCurtis_to_Artemia = art_summary$Value), art_dist))
add("Fig12C_plotted", fc[, keep_cols])
saveWorkbook(wb, file.path(OUT_TAB, "16S_ANCOMBC2_full_results.xlsx"), overwrite = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT_TAB, "R_session_info_figure12.txt"))
say("DONE: %s", file.path(OUT_TAB, "16S_ANCOMBC2_full_results.xlsx"))
