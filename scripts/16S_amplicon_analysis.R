# =============================================================================
#  16S rRNA AMPLICON ANALYSIS  -  Exaiptasia and Waminoa, fed vs starved
#  Master's thesis, Danjel Kola (Univ. Oldenburg / HIFMB). V5-V7 region, primers 799F / 1193R.
#
#  ONE SCRIPT, TWO PARTS
#    PART 1  DADA2      raw paired-end fastq  ->  ASV table  ->  SILVA taxonomy  ->  dada2_out/ps_raw.rds
#                       (about 2-3 h; every stage is checkpointed, so a stopped run resumes)
#    PART 2  ANALYSIS   ps_raw.rds  ->  contaminant removal  ->  bar plots, alpha/beta diversity,
#                       PERMANOVA, ANCOM-BC2  ->  tables/16S_all_tables.xlsx + figures/   (10-20 min)
#
#  Part 1 is skipped automatically when dada2_out/ps_raw.rds already exists (it is in the repository),
#  so the thesis results can be reproduced without the raw reads.
#
#  HOW TO RUN   from the repository root:     Rscript scripts/16S_amplicon_analysis.R
#               or open it in RStudio, set the working directory to the repository root and Source it.
#               Everything is written to 16S/ (dada2_out/, tables/, figures/).
#
#  INPUT   16S/dada2_out/ps_raw.rds     (supplied; or built by Part 1 from the raw reads)
#          raw reads (only Part 1)       ENA/SRA accession: see README.md, put the 50 fastq files in 16S/raw_fastq/
#          SILVA (only Part 1)           silva_nr99_v138.2_toGenus_trainset.fa.gz and silva_v138.2_assignSpecies.fa.gz
#                                        (Zenodo, doi:10.5281/zenodo.14169026) in 16S/reference_db/
#                                        optional for the 138.1 comparison sheet: silva_nr99_v138.1_train_set.fa.gz,
#                                        silva_species_assignment_v138.1.fa.gz (doi:10.5281/zenodo.4587955)
#
#  THESIS ANALYSIS = 17 animals. Five libraries (AF5, AS1, AS2, WF2, WS4) with 870-1,717 reads after cleaning
#  and 73-93 % primer-dimer reads are excluded (EXCLUDE_ANIMALS, Part 2). Run with EXCLUDE_ANIMALS="" for
#  all 20 animals (the sensitivity analysis).
#
#  PACKAGES  CRAN:         ggplot2, dplyr, tidyr, openxlsx, ggrepel, patchwork, ragg, vegan
#            Bioconductor: dada2, phyloseq, Biostrings, ShortRead, ANCOMBC, decontam
#            BiocManager::install(c("dada2","phyloseq","Biostrings","ShortRead","ANCOMBC","decontam"))
#  Versions used for the thesis: R 4.4.3, DADA2 1.34, ANCOMBC 2.8.0 (full list: 16S/R_session_info.txt).
# =============================================================================

# ---- CONFIG SHARED BY BOTH PARTS ------------------------------------------------
REPO    <- Sys.getenv("REPO_DIR", if (basename(getwd()) == "scripts") dirname(getwd()) else getwd())
PROJ    <- Sys.getenv("PROJ_DIR", file.path(REPO, "16S"))                # working folder for everything below
RAW_DIR <- Sys.getenv("RAW_DIR",  file.path(PROJ, "raw_fastq"))          # the 50 raw fastq files (Part 1 only)
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
RUN_DADA2 <- if (nzchar(Sys.getenv("RUN_DADA2"))) as.logical(Sys.getenv("RUN_DADA2")) else
             !file.exists(file.path(PROJ, "dada2_out", "ps_raw.rds"))
dir.create(PROJ, recursive = TRUE, showWarnings = FALSE)


# #############################################################################
#  PART 1   DADA2   raw fastq  ->  ASV table  ->  SILVA taxonomy  ->  phyloseq
#  Follows the official DADA2 tutorial (benjjneb.github.io/dada2): quality profiles -> filterAndTrim ->
#  learnErrors -> dada -> mergePairs -> makeSequenceTable -> removeBimeraDenovo -> track reads -> assignTaxonomy.
# #############################################################################
if (RUN_DADA2) {
cat("\n######## PART 1: DADA2 ########\n")
# Primers 799F / 1193R (V5-V7) are the first 19 / 18 bases of each read -> cut with trimLeft.
TRIM_LEFT   <- c(19, 18)
# truncLen is the position in the ORIGINAL read where it is cut (trimLeft is included):
# 240/200 keeps 221 + 182 bases. The V5-V7 insert is ~376 bp, so the pairs overlap by ~27 bp.
# DO NOT use 210/180: it keeps only 191 + 162 = 353 bp < 376 bp, the pairs cannot overlap and
# ~0 % of the reads merge (tested on 120,000 read pairs).
TRUNC_LEN   <- c(240, 200)
MAX_EE      <- c(2, 5)       # maximum expected errors (forward, reverse)
TRUNC_Q     <- 2
MIN_OVERLAP <- 12
INSERT_BP   <- 376           # modal insert length after primer removal (used for the safety check only)
LEN_WINDOW  <- 350:400       # keep merged sequences of this length
BINNED_Q    <- c(2, 11, 25, 37)   # the only quality values on these NovaSeq fastqs -> binned error model
NBASES      <- as.numeric(Sys.getenv("DADA2_NBASES", "1e9"))   # bases used to learn the error model (1e9 = all reads;
                                                               # restart_lighter.bat sets 3e8 = a random 30 %, still 3x DADA2's default)
THREADS     <- as.integer(Sys.getenv("THREADS", "8"))   # CPU threads (on Windows filterAndTrim cannot use threads)
RESUME      <- TRUE
N_QC_READS  <- 200000        # reads per library used to estimate the primer-dimer load

# Taxonomy: main = newest SILVA (138.2, Nov 2024). The old 138.1 is run as well, only for the
# comparison table "SILVA_comparison" in the Excel file. Set old to NULL to skip it.
SILVA_MAIN <- list(name  = "SILVA_138.2",
                   train = file.path(PROJ, "reference_db/silva_nr99_v138.2_toGenus_trainset.fa.gz"),
                   species = file.path(PROJ, "reference_db/silva_v138.2_assignSpecies.fa.gz"))
old_db_dir <- file.path(PROJ, "reference_db")
SILVA_OLD  <- list(name  = "SILVA_138.1",
                   train = file.path(old_db_dir, "silva_nr99_v138.1_train_set.fa.gz"),
                   species = file.path(old_db_dir, "silva_species_assignment_v138.1.fa.gz"))
if (!file.exists(SILVA_OLD$train)) SILVA_OLD <- NULL     # comparison is optional
# -----------------------------------------------------------------------------

suppressMessages({ library(dada2); library(phyloseq); library(Biostrings); library(ggplot2); library(ShortRead) })
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
out  <- file.path(PROJ, "dada2_out");        ckp <- file.path(out, "checkpoints")
filt <- file.path(out, "filtered")
for (d in c(out, ckp, filt)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time(); set.seed(100)

# run a stage once; on re-run load the saved result instead
stage <- function(name, expr) {
  f <- file.path(ckp, paste0(name, ".rds"))
  if (RESUME && file.exists(f)) { say("  [resume] %s loaded from checkpoint", name); return(readRDS(f)) }
  x <- expr; saveRDS(x, f); x
}

# ---- safety check: can the two reads still overlap? ---------------------------
len1 <- TRUNC_LEN[1] - TRIM_LEFT[1]; len2 <- TRUNC_LEN[2] - TRIM_LEFT[2]
overlap <- len1 + len2 - INSERT_BP
say("dada2 %s | kept read length %d + %d = %d bp | insert ~%d bp | expected overlap %d bp",
    packageVersion("dada2"), len1, len2, len1 + len2, INSERT_BP, overlap)
if (overlap < MIN_OVERLAP + 10)
  stop("truncLen leaves only ", overlap, " bp of overlap (< ", MIN_OVERLAP + 10, "): the read pairs will not merge. ",
       "Increase TRUNC_LEN (see the comment next to it).")

# ---- samples ------------------------------------------------------------------
{                                  # build the sample table from the BMK library ids
  smap <- c(M0001 = "AF1", M0002 = "AF2", M0003 = "AF3", M0004 = "AF4", M0005 = "AF5", M0006 = "AS1", M0007 = "AS2",
            M0008 = "AS3", M0009 = "AS4", M0010 = "AS5", M0011 = "ARTSEP", M0012 = "NEGCON1", M0013 = "WF1", M0014 = "WF2",
            M0015 = "WF3", M0016 = "WF4", M0017 = "WF5", M0018 = "WF6", M0019 = "WS1", M0020 = "WS2", M0021 = "WS3",
            M0022 = "WS4", M0023 = "WS5", M0024 = "WS6", M0025 = "NEGCON2")
  fq <- function(id, r) { f <- file.path(RAW_DIR, sprintf("Unknown_CS705-003%s_%d.fq", id, r))
    if (file.exists(f)) f else paste0(f, ".gz") }
  man <- data.frame(BMK_id = names(smap), Name = unname(smap),
                    source = ifelse(names(smap) %in% c("M0015", "M0016", "M0020", "M0021"), "resequenced_2026-09-22", "run2_original"),
                    R1 = vapply(names(smap), fq, "", r = 1), R2 = vapply(names(smap), fq, "", r = 2), reads_R1 = NA,
                    stringsAsFactors = FALSE, row.names = NULL)
}
stopifnot(nrow(man) == 25, !any(duplicated(man$Name)), all(file.exists(man$R1, man$R2)),
          file.exists(SILVA_MAIN$train), file.exists(SILVA_MAIN$species))
nm <- man$Name; fnFs <- setNames(man$R1, nm); fnRs <- setNames(man$R2, nm)
write.csv(man[, c("BMK_id", "Name", "source", "R1", "R2", "reads_R1")],
          file.path(out, "samples_used.csv"), row.names = FALSE)

# ---- 1. read quality (tutorial: plotQualityProfile) + primer-dimer load ------------
say("[1/8] quality profiles and primer-dimer load")
if (!file.exists(file.path(out, "quality_profiles_R1.pdf")))
  tryCatch({
    pdf(file.path(out, "quality_profiles_R1.pdf"), width = 11, height = 8); print(plotQualityProfile(unname(fnFs), aggregate = FALSE)); invisible(dev.off())
    pdf(file.path(out, "quality_profiles_R2.pdf"), width = 11, height = 8); print(plotQualityProfile(unname(fnRs), aggregate = FALSE)); invisible(dev.off())
  }, error = function(e) say("  ! quality-profile PDFs skipped: %s", conditionMessage(e)))

# Primer dimer = 799F joined directly to the reverse complement of 1193R. These reads are full length, so
# filterAndTrim() cannot see them; they are lost silently at mergePairs(). The share of forward reads that
# contain the reverse complement of 1193R within their first 70 bases is the "primer-dimer load" used in
# Appendix B.2 (Table B2). It is computed on the first N_QC_READS reads of each library.
DIMER <- "GGAAGGTGGGGATGACGT"     # reverse complement of 1193R (5'-ACGTCATCCCCACCTTCC-3')
qc_file <- file.path(out, "qc_summary.csv")
if (!file.exists(qc_file)) {
  qc_tab <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
    s1 <- ShortRead::sread(ShortRead::yield(ShortRead::FastqStreamer(man$R1[i], n = N_QC_READS)))
    data.frame(Name = man$Name[i], BMK_id = man$BMK_id[i], source = man$source[i], reads_scanned = length(s1),
               pct_primer_dimer = round(100 * mean(vcountPattern(DIMER, subseq(s1, 1, pmin(70, width(s1)))) > 0), 1))
  }))
  write.csv(qc_tab, qc_file, row.names = FALSE)
}

# ---- 2. filter and trim -------------------------------------------------------
say("[2/8] filterAndTrim  (trimLeft %s, truncLen %s, maxEE %s)",
    paste(TRIM_LEFT, collapse = "/"), paste(TRUNC_LEN, collapse = "/"), paste(MAX_EE, collapse = "/"))
filtFs <- setNames(file.path(filt, paste0(nm, "_F_filt.fastq.gz")), nm)
filtRs <- setNames(file.path(filt, paste0(nm, "_R_filt.fastq.gz")), nm)
out_filt <- stage("01_filter", {
  o <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = TRIM_LEFT, truncLen = TRUNC_LEN,
                     maxN = 0, maxEE = MAX_EE, truncQ = TRUNC_Q, rm.phix = TRUE,
                     compress = TRUE, multithread = if (.Platform$OS.type == "windows") FALSE else THREADS)
  rownames(o) <- nm; o })
print(out_filt)
stopifnot(all(file.exists(filtFs), file.exists(filtRs)))

# ---- 3. error model (binned quality scores, monotone) --------------------------
say("[3/8] learnErrors  (binned Q %s)", paste(BINNED_Q, collapse = "/"))
# NovaSeq reads carry only 4 quality values, so the default loess fit is unreliable. dada2 >= 1.34 ships
# makeBinnedQualErrfun(); the installed dada2 (1.32) does not, so the same idea is written out here:
# estimate the error rate at each binned Q, interpolate linearly (log10 scale) between the bins, then force
# every rate to be non-increasing with Q. Checked against dada2's own function on run-3 data: identical at the
# four binned Q values (difference 0), correlation 0.98-0.99 over all Q.
nts <- c("A", "C", "G", "T")
monotone_errfun <- function(trans) {
  qq <- as.numeric(colnames(trans))
  err <- matrix(0, 16, length(qq), dimnames = list(paste0(rep(nts, each = 4), "2", rep(nts, 4)), colnames(trans)))
  for (b in nts) {
    tot <- colSums(trans[paste0(b, "2", nts), , drop = FALSE])
    for (ntj in setdiff(nts, b)) {
      r <- paste0(b, "2", ntj); ok <- qq %in% BINNED_Q & tot > 0
      lp <- log10((trans[r, ok] + 1) / (tot[ok] + 1))                 # rate at each binned Q (1 pseudo-count)
      err[r, ] <- if (sum(ok) >= 2) 10^approx(qq[ok], lp, xout = qq, rule = 2)$y else rep(10^mean(lp), length(qq))
    }
  }
  err <- pmin(pmax(err, 1e-7), 0.25)
  for (b in nts) {
    subs <- paste0(b, "2", setdiff(nts, b))
    for (r in subs) err[r, ] <- rev(cummax(rev(err[r, ])))            # never rise with quality
    err[paste0(b, "2", b), ] <- 1 - colSums(err[subs, , drop = FALSE])
  }
  err
}
errs <- stage("02_errors", list(
  F = learnErrors(filtFs, nbases = NBASES, errorEstimationFunction = monotone_errfun,
                  multithread = THREADS, randomize = TRUE),
  R = learnErrors(filtRs, nbases = NBASES, errorEstimationFunction = monotone_errfun,
                  multithread = THREADS, randomize = TRUE)))
pdf(file.path(out, "error_plots.pdf"), width = 9, height = 7)
print(plotErrors(errs$F, nominalQ = TRUE) + ggtitle("R1 error model"))
print(plotErrors(errs$R, nominalQ = TRUE) + ggtitle("R2 error model"))
invisible(dev.off())

# ---- 4. sample inference --------------------------------------------------------
say("[4/8] dada (sample by sample, as in the tutorial)")
dd <- stage("03_dada", list(F = dada(filtFs, err = errs$F, multithread = THREADS),
                            R = dada(filtRs, err = errs$R, multithread = THREADS)))

# ---- 5. merge pairs, sequence table, length window --------------------------------
say("[5/8] mergePairs (minOverlap %d)", MIN_OVERLAP)
mg <- stage("04_merged", mergePairs(dd$F, filtFs, dd$R, filtRs, minOverlap = MIN_OVERLAP, verbose = TRUE))
st_all <- makeSequenceTable(mg)
len_tab <- as.data.frame(table(length = nchar(getSequences(st_all))))
len_tab$reads <- tapply(colSums(st_all), nchar(getSequences(st_all)), sum)[as.character(len_tab$length)]
write.csv(len_tab, file.path(out, "merged_length_distribution.csv"), row.names = FALSE)
st <- st_all[, nchar(colnames(st_all)) %in% LEN_WINDOW, drop = FALSE]
say("  length window %d-%d: %d of %d ASVs, %.1f %% of merged reads kept",
    min(LEN_WINDOW), max(LEN_WINDOW), ncol(st), ncol(st_all), 100 * sum(st) / sum(st_all))

# ---- 6. chimeras ------------------------------------------------------------------
say("[6/8] removeBimeraDenovo")
seqtab_chim <- stage("05_nochim", removeBimeraDenovo(st, method = "consensus", multithread = THREADS, verbose = TRUE))
say("  non-chimeric: %d ASVs, %.1f %% of reads kept", ncol(seqtab_chim), 100 * sum(seqtab_chim) / sum(st))

# sequencing artefacts: Illumina adapter + poly-G / poly-C reads that survive merging
adapters <- c(P5 = "AATGATACGGCGACCACCGA", P7 = "ATCTCGTATGCCGTCTTCTGCTTG",
              TruSeq_R1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT", TruSeq_trim = "AGATCGGAAGAGC",
              Nextera = "CTGTCTCTTATACACATCT")
adapters <- c(adapters, setNames(as.character(reverseComplement(DNAStringSet(adapters))), paste0(names(adapters), "_rc")))
sq <- colnames(seqtab_chim)
hit_ad   <- vapply(sq, function(s) paste(names(adapters)[vapply(adapters, grepl, logical(1), x = s, fixed = TRUE)], collapse = ";"), "")
hit_homo <- grepl("G{12,}|C{12,}|A{12,}|T{12,}", sq)
artefact <- nzchar(hit_ad) | hit_homo
write.csv(data.frame(sequence = sq[artefact], length = nchar(sq[artefact]), reads = colSums(seqtab_chim)[artefact],
                     adapter_hits = hit_ad[artefact], homopolymer_ge12 = hit_homo[artefact]),
          file.path(out, "removed_artefact_ASVs.csv"), row.names = FALSE)
say("  adapter / homopolymer artefact ASVs removed: %d ASVs, %d reads", sum(artefact), sum(seqtab_chim[, artefact]))
seqtab <- seqtab_chim[, !artefact, drop = FALSE]
seqtab <- seqtab[, order(colSums(seqtab), decreasing = TRUE), drop = FALSE]   # ASV0001 = most abundant

# ---- 7. read tracking ---------------------------------------------------------------
say("[7/8] read tracking")
gN <- function(x) sum(getUniques(x))
track <- data.frame(Name = nm, input = out_filt[nm, 1], filtered = out_filt[nm, 2],
                    denoisedF = sapply(dd$F[nm], gN), denoisedR = sapply(dd$R[nm], gN),
                    merged = sapply(mg[nm], gN), len_filt = rowSums(st)[nm],
                    nonchim = rowSums(seqtab_chim)[nm], final = rowSums(seqtab)[nm])
track$pct_kept <- round(100 * track$final / track$input, 1)
track$source <- man$source[match(nm, man$Name)]
write.csv(track, file.path(out, "read_tracking.csv"), row.names = FALSE)
print(track, row.names = FALSE)

# ---- 8. taxonomy (main database + old one for the comparison) ---------------------
classify <- function(db) {
  say("[8/8] assignTaxonomy + addSpecies with %s", db$name)
  stage(paste0("06_taxa_", db$name), {
    tx <- assignTaxonomy(seqtab, db$train, multithread = THREADS, outputBootstraps = TRUE)
    list(tax = addSpecies(tx$tax, db$species), boot = tx$boot)
  })
}
tx_main <- classify(SILVA_MAIN)
tx_old  <- if (!is.null(SILVA_OLD)) classify(SILVA_OLD) else NULL

asv_id <- sprintf("ASV%04d", seq_len(ncol(seqtab)))
dna <- DNAStringSet(colnames(seqtab)); names(dna) <- asv_id
writeXStringSet(dna, file.path(out, "ASV_sequences.fasta"))
rename_tax <- function(tx) { tx$tax <- tx$tax[colnames(seqtab), , drop = FALSE]; rownames(tx$tax) <- asv_id
  tx$boot <- tx$boot[colnames(seqtab), , drop = FALSE]; rownames(tx$boot) <- asv_id; tx }
tx_main <- rename_tax(tx_main)
saveRDS(tx_main$tax, file.path(out, "taxa_main.rds"))
if (!is.null(tx_old)) { tx_old <- rename_tax(tx_old); saveRDS(tx_old$tax, file.path(out, "taxa_old_db.rds")) }
write.csv(data.frame(ASV = asv_id, tx_main$tax, boot = tx_main$boot, sequence = as.character(dna)),
          file.path(out, "ASV_taxonomy.csv"), row.names = FALSE)
seq_mat <- seqtab; colnames(seq_mat) <- asv_id
saveRDS(seq_mat, file.path(out, "seqtab_nochim.rds"))

# ---- phyloseq object -----------------------------------------------------------------
sdf <- data.frame(Name = nm, BMK_id = man$BMK_id, source = man$source, row.names = nm)
ps <- phyloseq(otu_table(seq_mat, taxa_are_rows = FALSE), sample_data(sdf[rownames(seq_mat), ]),
               tax_table(tx_main$tax), refseq(dna))
saveRDS(ps, file.path(out, "ps_raw.rds"))
writeLines(c(sprintf("trimLeft %s | truncLen %s | maxEE %s | minOverlap %d | length window %d-%d",
                     paste(TRIM_LEFT, collapse = "/"), paste(TRUNC_LEN, collapse = "/"), paste(MAX_EE, collapse = "/"),
                     MIN_OVERLAP, min(LEN_WINDOW), max(LEN_WINDOW)),
             sprintf("taxonomy: %s (+ %s for comparison)", SILVA_MAIN$name, if (is.null(SILVA_OLD)) "none" else SILVA_OLD$name),
             capture.output(sessionInfo())), file.path(out, "dada2_settings_and_session.txt"))
say("DONE in %.1f h  ->  %s", as.numeric(difftime(Sys.time(), t0, units = "hours")), file.path(out, "ps_raw.rds"))

}   # end of PART 1
rm(list = setdiff(ls(), c("REPO", "PROJ", "RAW_DIR", "THREADS", "RUN_DADA2")))   # PART 2 starts from a clean workspace


# #############################################################################
#  PART 2   ANALYSIS   phyloseq object  ->  contaminant removal  ->  every table (Excel)
#                      + bar plots, PCoA, PERMANOVA, ANCOM-BC2, alpha diversity, panel A/B/C figure
#
#  THE CONTAMINANT RULES (in this order; every ASV gets exactly one status)
#    1   not bacterial / no phylum / chloroplast / mitochondria
#    3   MANUAL BLOCKLIST: human-handling genera + families (blocklist section below)
#    3b  human oral / gut genera (habitat list in the blocklist section)
#    4   reagent / water-system genera - removed ONLY if the same genus is in the blank NEGCON1
#    5   "blank-dominant" ASVs - more abundant (relative) in NEGCON1 than in ANY animal
#        (unless Artemia carries it even higher: then it is cross-talk from the food library)
#  Only NEGCON1 is used; NEGCON2 is left out (too few usable reads). The sheet "NEGCON2_sensitivity"
#  shows what would change if it were used as well. The sheet "Why_NEGCON1" explains the blank.
# #############################################################################
cat("\n######## PART 2: ANALYSIS ########\n")

# =============================================================================
#  CONTAMINANT BLOCKLISTS  -  edit here to change what counts as a contaminant
#  Names must match the Genus / Family column of the SILVA taxonomy. A family-only assignment is
#  written "NA_<Family>" (e.g. "NA_Microbacteriaceae"). SILVA 138.2 renamed some phyla
#  (Proteobacteria -> Pseudomonadota, Firmicutes -> Bacillota ...), but genus names are mostly unchanged.
#  The "Blocklist_coverage" sheet shows which names were / were not found in the data.
# =============================================================================
# ---- STEP 3: manual "human handling" blocklist (skin, mouth, gut, kit/lab) -----
MANUAL_GENERA <- c(
  "Cutibacterium","Propionibacterium","Staphylococcus","Corynebacterium",
  "Lawsonella","Micrococcus","Kocuria","Dermabacter","Brevibacterium",
  "Streptococcus","Rothia","Neisseria","Haemophilus","Veillonella",
  "Fusobacterium","Porphyromonas","Prevotella","Prevotella_7","Gemella",
  "Actinomyces","Granulicatella","Leptotrichia","Anaerococcus",
  "Peptoniphilus","Finegoldia","Faecalibacterium","Bacteroides",
  "Blautia","Bifidobacterium","Lactobacillus","Enterococcus",
  "Ralstonia","Bradyrhizobium","Burkholderia","Delftia","Herbaspirillum",
  "Methylobacterium","Cupriavidus",
  "Acinetobacter","Stenotrophomonas","Arthrobacter","Rhodococcus",
  "Microbacterium","Curtobacterium","Aeromicrobium","Janibacter",
  "Dietzia","Tsukamurella","Microlunatus","Patulibacter","Beutenbergia",
  "Paenibacillus","Brevibacillus","Facklamia","Brochothrix","Abiotrophia",
  # SILVA 138 names of some of the above
  "Burkholderia-Caballeronia-Paraburkholderia","Methylobacterium-Methylorubrum",
  "Prevotella_9","Alloprevotella")

MANUAL_FAMILIES <- c(
  "Staphylococcaceae","Corynebacteriaceae","Propionibacteriaceae",
  "Micrococcaceae","Bifidobacteriaceae","Pasteurellaceae",
  "Bacteroidaceae","Actinomycetaceae","Veillonellaceae",
  "Streptococcaceae","Lactobacillaceae","Family XI",
  "Dermacoccaceae","Bogoriellaceae")

# ---- habitat lists (literature based; Salter 2014, Eisenhofer 2019 and others) --
# MARINE        protected: never removed by the reagent step (real sea-water bacteria)
# HUMAN_ORAL_GUT -> STEP 3b: removed (oral-plaque / gut genera have no source in an anemone)
# REAGENT       -> STEP 4: removed ONLY if the genus is also found in the extraction blank
# SOIL          informational only (not removed)
MARINE <- c(
    "Alteromonas", "Pseudoalteromonas", "Vibrio", "Marinobacter", "Oceanibulbus", "Alcanivorax",
    "Halioxenophilus", "Erythrobacter", "Halomonas", "JTB255_marine_benthic_group",
    "NA_UBA10353_marine_group", "Aurantivirga", "Thalassospira", "Donghicola", "Aestuariibacter",
    "Lewinella", "Spongiibacter", "Halobacteriovorax", "NA_Hyphomonadaceae", "Photobacterium",
    "Sphingorhabdus", "Marinovum", "NA_Halieaceae", "NA_Dadabacteriales", "AqS1", "SM1A02",
    "NA_Alteromonadaceae", "Fodinicurvata", "NA_EPR3968-O8a-Bc78", "Ruegeria", "Endozoicomonas",
    "Nautella", "Tenacibaculum", "Roseovarius", "Sulfitobacter", "Marivita", "Pseudophaeobacter",
    "Phaeobacter", "Maribacter", "Winogradskyella", "Kordiimonas", "Neptuniibacter", "Oleiphilus",
    "Oleibacter", "Thalassotalea", "Colwellia", "Shewanella", "Algibacter", "Polaribacter",
    "Glaciecola", "Owenweeksia", "Fulvivirga", "Rubritalea", "Pelagibius", "Labrenzia", "Stappia",
    "Maritimibacter", "Litoreibacter", "Lutibacter", "Ahrensia", "Filomicrobium", "Arenibacter",
    "Muricauda", "Aliiroseovarius", "Tenuibacillus", "Marinicella", "Pseudohongiella", "Amphritea",
    "Nitrosopumilus", "Psychrosphaera", "NA_Rhodobacteraceae", "Thalassobius", "Oceanicaulis",
    "Hyphomonas", "Robiginitomaculum", "Maricaulis", "Aliivibrio", "Idiomarina",
    "Pseudospongiibacter", "Porticoccus", "NA_Flavobacteriaceae", "Croceitalea", "Aquimarina",
    "Cohaesibacter", "Shimia", "NA_Saprospiraceae", "Portibacter", "Haliea", "Congregibacter",
    "NA_Nitrincolaceae", "Marinomonas", "Oceanospirillum", "Sediminicola", "Salinirepens",
    "Candidatus_Endoecteinascidia", "Leisingera", "Loktanella", "Octadecabacter", "Planktotalea",
    "Pseudomonas_marine", "Woeseia", "NA_Woeseiaceae", "Bdellovibrio_marine", "NA_Arenicellaceae",
    "Arenicella", "Halocynthiibacter", "Kiloniella", "NA_Kiloniellaceae", "Pelagibaca", "Roseibium"
  )

HUMAN_ORAL_GUT <- c(
    "Treponema", "Catonella", "Filifactor", "Prevotella_9", "Alloprevotella", "Tannerella",
    "Lautropia", "F0058", "Peptostreptococcus", "[Eubacterium]_brachy_group", "Capnocytophaga",
    "Eikenella", "Johnsonella", "Moryella", "Peptoanaerobacter", "Alysiella", "Bergeyella",
    "Caviibacter", "Hafnia-Obesumbacterium", "Enterobacter", "Phocaeicola",
    "Rikenellaceae_RC9_gut_group", "NA_Ruminococcaceae", "NA_Lachnospiraceae",
    "Clostridium_sensu_stricto_1", "NA_[Eubacterium]_coprostanoligenes_group",
    "Defluviitaleaceae_UCG-011", "NA_Clostridia_vadinBB60_group", "Aerococcus", "NA_Neisseriaceae",
    "NA_Prevotellaceae", "NA_Lactobacillales", "Brachybacterium", "Ethanoligenens",
    "Escherichia-Shigella", "Parvimonas", "Solobacterium", "Oribacterium", "Lachnoanaerobaculum",
    "Stomatobaculum", "Megasphaera", "Dialister", "Selenomonas", "Campylobacter", "Mogibacterium",
    "Atopobium", "Kingella", "Cardiobacterium", "Aggregatibacter", "Tannerellaceae", "Bacteroides",
    "Parabacteroides", "Alistipes", "Ruminococcus", "Subdoligranulum", "Dorea", "Anaerostipes",
    "Roseburia", "Collinsella", "Akkermansia", "Muribaculaceae", "NA_Muribaculaceae",
    "Lachnoclostridium", "Butyricicoccus", "Erysipelotrichaceae_UCG-003", "Helicobacter",
    "Proteus", "Klebsiella", "Serratia", "Citrobacter", "Morganella", "Candidatus_Saccharimonas",
    "Leptotrichia_like", "Gemella_like", "Actinomyces_like", "Dermacoccus", "Enhydrobacter"
  )

REAGENT <- c(
    "Pseudomonas", "Aquabacterium", "Pelomonas", "Sphingomonas", "Cloacibacterium", "Bacillus",
    "Curvibacter", "Methylibium", "Paracoccus", "Deinococcus", "Hyphomicrobium",
    "Chryseobacterium", "Mycobacterium", "Massilia", "Limnobacter", "Brevundimonas", "Tepidimonas",
    "Pseudoxanthomonas", "Methylotenera", "Roseomonas", "Leptothrix", "Flavobacterium",
    "Noviherbaspirillum", "Sphingobacterium", "Variovorax", "Methylophilus",
    "Methylobacterium-Methylorubrum", "Afipia", "Sphingobium", "Methyloversatilis",
    "Undibacterium", "Meiothermus", "Thermomonas", "Tepidiphilus", "NA_Microbacteriaceae",
    "Xanthomonas", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "NA_Caulobacteraceae",
    "NA_Comamonadaceae", "Piscinibacter", "Acidiphilium", "Aureimonas", "Reyranella",
    "Rhizobacter", "Spirosoma", "1174-901-12", "Nakamurella", "Rubellimicrobium", "Amaricoccus",
    "Tabrizicola", "Rhodovarius", "Anaerobacillus", "NA_Bacillales", "NA_Bacillaceae",
    "Caldalkalibacillus", "Truepera", "Novosphingobium", "Sphingopyxis", "Bosea", "Caulobacter",
    "Devosia", "Mesorhizobium", "Ochrobactrum", "Phyllobacterium", "Acidovorax", "Comamonas",
    "Duganella", "Herbaspirillum", "Janthinobacterium", "Polaromonas", "Hydrogenophaga",
    "Ideonella", "Rhodoferax", "Schlegelella", "Sulfuritalea", "Psychrobacter", "Legionella",
    "Dyadobacter", "Pedobacter", "Hydrotalea", "Aquicella", "Obscuribacter", "Lysobacter",
    "Achromobacter", "Pseudorhodoferax", "Rhizorhapis", "Delftia_like", "Bradyrhizobium_like",
    "Cutibacterium_like", "Rheinheimera", "Rhodobacter", "Gemmobacter", "Hymenobacter",
    "Porphyrobacter", "Altererythrobacter", "Erythromicrobium", "Blastomonas", "Sphingosinicella",
    "Asticcacaulis"
  )

SOIL <- c(
    "Nocardioides", "Gaiella", "NA_Gaiellales", "NA_67-14", "Solirubrobacter",
    "NA_Pedosphaeraceae", "Ellin6055", "Ellin6067", "MND1", "NA_Gemmatimonadaceae", "Subgroup_10",
    "NA_Subgroup_7", "Blastococcus", "Pseudonocardia", "Actinomycetospora", "Rubrobacter",
    "Rhodoplanes", "Marmoricola", "Lapillicoccus", "Conexibacter", "Streptomyces",
    "NA_Xanthobacteraceae", "NA_TRA3-20", "NA_IMCC26256", "NA_Microtrichales", "NA_Frankiales",
    "RB41", "wb1-P19", "Nordella", "Arsenicicoccus", "NA_Micrococcales", "Bryobacter",
    "Candidatus_Udaeobacter", "Haliangium_soil", "NA_Vicinamibacteraceae", "Vicinamibacteraceae",
    "NA_Solirubrobacteraceae", "Mycobacterium_soil", "Geodermatophilus", "Modestobacter",
    "Kribbella", "Microvirga", "Skermanella", "Dongia", "Pirellula", "Terrimonas",
    "Flavisolibacter", "Niastella", "Ferruginibacter", "Chitinophaga", "Segetibacter",
    "NA_Chitinophagaceae", "Candidatus_Nitrosotalea", "Nitrospira", "Nitrosospira",
    "Steroidobacter", "Ramlibacter", "Microlunatus_like", "Kineosporia", "Friedmanniella", "Iamia",
    "NA_Iamiaceae", "Crossiella", "Amycolatopsis", "Nonomuraea", "Actinoplanes",
    "Dactylosporangium", "Catellatospora", "Mycobacterium_soil2", "Pseudarthrobacter",
    "Paenarthrobacter"
  )



# ---- SETTINGS ----------------------------------------------------------------
PS_RAW      <- Sys.getenv("PS_RAW",  file.path(PROJ, "dada2_out/ps_raw.rds"))
OUT         <- Sys.getenv("OUT_DIR", PROJ)

MIN_READS   <- 1000          # animals with fewer reads AFTER cleaning are removed from the analysis
# extra animals to leave out whatever their read count, e.g. "AF5,AS1,AS2,WF2,WS4" (the 5 high-primer-dimer libraries
# excluded in the thesis analysis, n = 17). Can also be set from outside: environment variable EXCLUDE_ANIMALS.
EXCLUDE_ANIMALS <- strsplit(Sys.getenv("EXCLUDE_ANIMALS", "AF5,AS1,AS2,WF2,WS4"), ",")[[1]]   # thesis analysis (n = 17). Use EXCLUDE_ANIMALS="" for all 20 animals.
BLANKS      <- "NEGCON1"     # extraction blank(s) used for contaminant rules 4 and 5
DROP_ALWAYS <- "NEGCON2"     # library left out completely
ARTEMIA     <- "ARTSEP"      # the food library (drawn at the side of the bar plots)
TOP_N       <- 15            # families / genera shown in bar plots
SHOW_BLANK_IN_BEFORE <- TRUE # also draw NEGCON1 next to Artemia in the "before removal" bars
ANCOM_RANKS <- c("Genus", "Family")
ANCOM_PREV  <- 0.30          # taxon must be present in >= 30 % of samples to be tested
ALPHA       <- 0.05          # significance level
N_PERM      <- 9999          # PERMANOVA permutations
FIG_FORMATS <- c("pdf", "png")   # add "tiff" if a journal asks for it (large files)
FIG_DPI     <- 600

# optional inputs for comparison tables (skipped automatically if missing)
QC_CSV     <- Sys.getenv("QC_CSV",     file.path(PROJ, "dada2_out/qc_summary.csv"))   # primer-dimer % per library (written in PART 1)
PREV_TRACK <- Sys.getenv("PREV_TRACK", "")   # optional: read tracking of an earlier DADA2 run, for the run-comparison sheet

GROUP_COL <- c("Aiptasia Fed" = "#1F5FA8", "Aiptasia Starved" = "#8FB8E6",
               "Waminoa Fed"  = "#C4501B", "Waminoa Starved"  = "#F2AD85")
SIDE_COL  <- c("Artemia (food)" = "#1B9E77", "Extraction blank" = "#6B6B66")
# bar-plot colours: 16 distinct colours first, then
# extra colours for taxa that are only in the "before" top list. "Other / Unclassified" is one grey category.
TAXON_PALETTE <- c("#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#4a3aa7","#e34948",
                   "#00a3c4","#9c5518","#b83dbb","#6f9c00","#d4006e","#0f5fa8","#c97f2e","#7a3fd6",
                   "#5f8f8f","#b5a642","#8c564b","#2f4b7c","#ff9896","#17becf","#98df8a","#c5b0d5",
                   "#bcbd22","#aec7e8","#ffbb78","#9467bd","#393b79","#e7969c")
OTHER_COL   <- "#c9c8c2"; OTHER_LABEL <- "Other / Unclassified"
STATUS_COL <- c("Non-bacterial / chloroplast / mitochondria" = "#7F7F7F",
                "Manual blocklist (human handling)"          = "#D4006E",
                "Human oral / gut genera"                    = "#8E44AD",
                "Reagent genera found in NEGCON1"            = "#E0A33A",
                "Blank-dominant ASVs (NEGCON1)"              = "#9A8F7A",
                "Retained (decontaminated)"                  = "#1F5FA8")

# ---- packages & helpers --------------------------------------------------------
pkgs <- c("phyloseq", "ANCOMBC", "vegan", "ggplot2", "dplyr", "tidyr", "openxlsx", "ggrepel", "patchwork", "ragg")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Please install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages(for (p in pkgs) library(p, character.only = TRUE))
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
for (d in c("tables", "figures/main", "figures/supplementary")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, name, w, h, sub = "supplementary") {
  for (fm in FIG_FORMATS) {
    f <- file.path(OUT, "figures", sub, paste0(name, ".", fm))
    dev <- switch(fm, pdf = grDevices::cairo_pdf, png = ragg::agg_png, tiff = ragg::agg_tiff, NULL)
    tryCatch(if (fm == "pdf") ggsave(f, p, width = w, height = h, device = dev)
             else ggsave(f, p, width = w, height = h, dpi = FIG_DPI, device = dev),
             error = function(e) say("  ! could not write %s (open in another program?): %s", basename(f), conditionMessage(e)))
  }
  say("  figure: %s/%s", sub, name)
}
theme_pub <- function(bs = 9) theme_classic(base_size = bs) +
  theme(axis.line = element_line(linewidth = 0.35, colour = "grey20"),
        axis.ticks = element_line(linewidth = 0.35, colour = "grey20"),
        axis.text = element_text(colour = "grey15"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(face = "bold", size = bs - 0.5),
        plot.title = element_text(face = "bold", size = bs + 1.5),
        plot.title.position = "plot", plot.subtitle = element_text(colour = "grey25"),
        legend.key.size = unit(9, "pt"), legend.title = element_text(face = "bold"),
        panel.spacing.x = unit(0.5, "lines"))
theme_set(theme_pub())
counts <- function(ps) { m <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) t(m) else m }
pstar <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))


# =============================================================================
# 1. LOAD DATA AND SAMPLE TABLE
# =============================================================================
say("1. Loading %s", PS_RAW)
ps0 <- readRDS(PS_RAW)
if (taxa_are_rows(ps0)) otu_table(ps0) <- t(otu_table(ps0))
ps0 <- prune_samples(setdiff(sample_names(ps0), DROP_ALWAYS), ps0)
sn  <- sample_names(ps0)
grp_of <- function(n) ifelse(n == ARTEMIA, "Artemia (food)", ifelse(grepl("^NEG", n), "Extraction blank",
                  paste(ifelse(substr(n, 1, 1) == "A", "Aiptasia", "Waminoa"), ifelse(substr(n, 2, 2) == "F", "Fed", "Starved"))))
ord_all <- intersect(c(paste0("AF", 1:6), paste0("AS", 1:6), paste0("WF", 1:6), paste0("WS", 1:6), ARTEMIA, BLANKS), sn)
ord_ani <- setdiff(ord_all, c(ARTEMIA, BLANKS))
stopifnot(all(BLANKS %in% sn), ARTEMIA %in% sn, length(ord_ani) >= 10)
qc <- if (file.exists(QC_CSV)) read.csv(QC_CSV, stringsAsFactors = FALSE) else NULL
meta_all <- data.frame(Sample = sn, Group = grp_of(sn), row.names = sn, stringsAsFactors = FALSE)
meta_all$Organism  <- ifelse(meta_all$Group %in% names(GROUP_COL), sub(" .*", "", meta_all$Group), meta_all$Group)
meta_all$Treatment <- ifelse(meta_all$Group %in% names(GROUP_COL), sub(".* ", "", meta_all$Group), NA)
meta_all$Primer_dimer_pct <- if (!is.null(qc)) qc$pct_primer_dimer[match(sn, qc$Name)] else NA
cnt0 <- counts(ps0)[ord_all, , drop = FALSE]; ps0 <- prune_samples(ord_all, ps0)
tt0  <- as.data.frame(as(tax_table(ps0), "matrix"), stringsAsFactors = FALSE)
say("   %d libraries, %d ASVs, %s reads", nrow(cnt0), ncol(cnt0), format(sum(cnt0), big.mark = ","))


# =============================================================================
# 2. CONTAMINANT REMOVAL
# =============================================================================
say("2. Contaminant removal")
rel0 <- cnt0 / pmax(rowSums(cnt0), 1)
add_labels <- function(tt) {
  tt$label <- ifelse(!is.na(tt$Genus), tt$Genus, ifelse(!is.na(tt$Family), paste0("NA_", tt$Family), "Unassigned"))
  l <- gsub(" ", "_", tt$label)
  tt$habitat <- ifelse(l %in% MARINE, "Marine", ifelse(l %in% HUMAN_ORAL_GUT, "Human", ifelse(l %in% REAGENT, "Reagent",
                ifelse(l %in% SOIL, "Soil", "Other"))))
  tt
}
STEP_LAB <- unname(names(STATUS_COL))
# status (1 label per ASV) for a taxonomy table and a set of blanks
get_status <- function(tt, blanks) {
  tt <- add_labels(tt)
  s1  <- tt$Kingdom != "Bacteria" | is.na(tt$Kingdom) | is.na(tt$Phylum) | tt$Order %in% "Chloroplast" | tt$Family %in% "Mitochondria"
  s1[is.na(s1)] <- TRUE
  s3  <- !s1 & (tt$Genus %in% MANUAL_GENERA | tt$Family %in% MANUAL_FAMILIES)
  s3b <- !s1 & !s3 & tt$habitat == "Human"
  in_blank <- colSums(cnt0[blanks, , drop = FALSE])
  s4  <- !s1 & !s3 & !s3b & tt$habitat == "Reagent" & tt$label %in% unique(tt$label[in_blank > 0])
  bmax <- apply(rel0[blanks, , drop = FALSE], 2, max)
  amax <- apply(rel0[ord_ani, , drop = FALSE], 2, max)
  art  <- rel0[ARTEMIA, ]
  s5  <- !s1 & !s3 & !s3b & !s4 & bmax > amax & !(bmax > 0 & art > bmax)
  st <- ifelse(s1, STEP_LAB[1], ifelse(s3, STEP_LAB[2], ifelse(s3b, STEP_LAB[3], ifelse(s4, STEP_LAB[4],
        ifelse(s5, STEP_LAB[5], STEP_LAB[6])))))
  list(status = factor(st, levels = STEP_LAB), tt = tt, bmax = bmax, amax = amax, art = art, in_blank = in_blank)
}
cl <- get_status(tt0, BLANKS)
status <- cl$status; tt <- cl$tt; KEPT <- STEP_LAB[6]
keep_asv <- status == KEPT
stopifnot(sum(keep_asv) > 0)

# NEGCON2 sensitivity: what would change if both blanks (when NEGCON2 exists in the raw object) were used
cl_alt <- NULL
ps_full <- readRDS(PS_RAW); if (taxa_are_rows(ps_full)) otu_table(ps_full) <- t(otu_table(ps_full))
if (length(DROP_ALWAYS) && all(DROP_ALWAYS %in% sample_names(ps_full))) {
  cnt_save <- cnt0; rel_save <- rel0                      # get_status() reads cnt0 / rel0, so swap them briefly
  cnt0 <- counts(ps_full)[c(ord_all, DROP_ALWAYS), colnames(cnt_save)]
  rel0 <- cnt0 / pmax(rowSums(cnt0), 1)
  cl_alt <- get_status(tt0, c(BLANKS, DROP_ALWAYS))
  cnt0 <- cnt_save; rel0 <- rel_save
}
rm(ps_full)

# cleaned reads per library, minimum-depth filter
reads_raw   <- rowSums(cnt0)
reads_clean <- rowSums(cnt0[, keep_asv, drop = FALSE])
low_reads <- ord_ani[reads_clean[ord_ani] < MIN_READS]
manual_ex <- intersect(trimws(EXCLUDE_ANIMALS), ord_ani)
low <- union(low_reads, manual_ex)                      # every animal left out, for whatever reason
ani <- setdiff(ord_ani, low)
if (length(manual_ex)) say("   also excluded by choice (EXCLUDE_ANIMALS): %s", paste(manual_ex, collapse = ", "))
say("   kept %d of %d ASVs; reads retained %.1f %% (animals); removed for < %d reads: %s",
    sum(keep_asv), ncol(cnt0), 100 * sum(reads_clean[ani]) / sum(reads_raw[ani]), MIN_READS,
    if (length(low)) paste(low, collapse = ", ") else "none")

meta <- meta_all[ani, ]
meta$Organism  <- factor(meta$Organism,  levels = c("Aiptasia", "Waminoa"))
meta$Treatment <- factor(meta$Treatment, levels = c("Fed", "Starved"))
meta$Group     <- factor(meta$Group, levels = names(GROUP_COL))
meta <- meta[order(meta$Organism, meta$Treatment, rownames(meta)), ]
ani  <- rownames(meta)
ps   <- prune_taxa(colnames(cnt0)[keep_asv], prune_samples(ani, ps0))
sample_data(ps) <- sample_data(meta); ps <- prune_taxa(taxa_sums(ps) > 0, ps)
cat("\n"); print(table(meta$Organism, meta$Treatment)); cat("\n")
g_n <- table(meta$Group)

samples_tab <- data.frame(
  Sample = ord_all, Group = meta_all[ord_all, "Group"],
  Reads_after_DADA2 = as.integer(reads_raw), Reads_after_cleaning = as.integer(reads_clean),
  Retained_pct = round(100 * reads_clean / reads_raw, 1), Primer_dimer_pct_R1 = meta_all[ord_all, "Primer_dimer_pct"],
  Status = ifelse(ord_all %in% ani, "analysed", ifelse(ord_all %in% low_reads, paste0("REMOVED (< ", MIN_READS, " reads after cleaning)"),
           ifelse(ord_all %in% manual_ex, "REMOVED (excluded by choice, EXCLUDE_ANIMALS)",
                  ifelse(ord_all == ARTEMIA, "Artemia food library (reference only)", "extraction blank (contaminant reference)")))),
  row.names = NULL, check.names = FALSE)

# ---- contamination per library ----------------------------------------------------
by_status <- sapply(STEP_LAB, function(s) rowSums(cnt0[, status == s, drop = FALSE]))   # libraries x 6 statuses
contam_ps <- data.frame(Sample = ord_all, Group = meta_all[ord_all, "Group"],
                        Analysed = ifelse(ord_all %in% ani, "yes", "no"),
                        Reads_before_cleaning = rowSums(by_status), by_status[, 1:5], Reads_removed_total = rowSums(by_status[, 1:5]),
                        Contamination_pct = round(100 * rowSums(by_status[, 1:5]) / rowSums(by_status), 1),
                        Reads_retained = by_status[, 6], Retained_pct = round(100 * by_status[, 6] / rowSums(by_status), 1),
                        ASVs_before_cleaning = rowSums(cnt0 > 0), ASVs_retained = rowSums(cnt0[, keep_asv, drop = FALSE] > 0),
                        row.names = NULL, check.names = FALSE)
tot_row <- function(lab, smp) { d <- by_status[smp, , drop = FALSE]
  data.frame(Sample = lab, Group = "", Analysed = "", Reads_before_cleaning = sum(d), t(colSums(d[, 1:5, drop = FALSE])),
             Reads_removed_total = sum(d[, 1:5]), Contamination_pct = round(100 * sum(d[, 1:5]) / sum(d), 1),
             Reads_retained = sum(d[, 6]), Retained_pct = round(100 * sum(d[, 6]) / sum(d), 1),
             ASVs_before_cleaning = NA, ASVs_retained = NA, check.names = FALSE) |> setNames(names(contam_ps)) }
contam_tab <- rbind(contam_ps, tot_row("TOTAL analysed animals", ani),
                    tot_row("TOTAL Aiptasia", ani[meta$Organism == "Aiptasia"]), tot_row("TOTAL Waminoa", ani[meta$Organism == "Waminoa"]))
print(contam_tab[grepl("^TOTAL", contam_tab$Sample), c("Sample", "Reads_before_cleaning", "Contamination_pct", "Reads_retained")])

seen <- colSums(cnt0[ani, ]) > 0
step_tab <- do.call(rbind, lapply(STEP_LAB, function(s) { i <- seen & status == s
  data.frame(Filter_step = s, ASVs = sum(i), ASVs_pct = round(100 * sum(i) / sum(seen), 1),
             Reads = sum(cnt0[ani, i]), Reads_pct = round(100 * sum(cnt0[ani, i]) / sum(cnt0[ani, seen]), 1)) }))
step_tab <- rbind(data.frame(Filter_step = "ALL ASVs before cleaning (analysed animals)", ASVs = sum(seen), ASVs_pct = 100,
                             Reads = sum(cnt0[ani, seen]), Reads_pct = 100),
                  step_tab[1:5, ],
                  data.frame(Filter_step = "TOTAL removed as contaminant", ASVs = sum(step_tab$ASVs[1:5]),
                             ASVs_pct = round(sum(step_tab$ASVs_pct[1:5]), 1), Reads = sum(step_tab$Reads[1:5]),
                             Reads_pct = round(sum(step_tab$Reads_pct[1:5]), 1)),
                  data.frame(Filter_step = "LEFT AFTER REMOVAL", step_tab[6, -1]))
print(step_tab)

# every ASV with its fate (audit table)
dec <- tryCatch(decontam::isContaminant(cnt0[c(ani, BLANKS), ], neg = c(rep(FALSE, length(ani)), rep(TRUE, length(BLANKS))),
                                        method = "prevalence", threshold = 0.5), error = function(e) NULL)
asv_status <- data.frame(ASV = colnames(cnt0), tt0[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")],
                         Status = as.character(status), Habitat_list = tt$habitat,
                         Reads_in_analysed_animals = colSums(cnt0[ani, ]), Reads_in_NEGCON1 = colSums(cnt0[BLANKS, , drop = FALSE]),
                         Max_rel_pct_animal = round(100 * cl$amax, 3), Max_rel_pct_NEGCON1 = round(100 * cl$bmax, 3),
                         Rel_pct_Artemia = round(100 * cl$art, 3),
                         decontam_p = if (!is.null(dec)) round(dec$p, 4) else NA, t(cnt0[ani, ]),
                         check.names = FALSE, row.names = NULL)
asv_status <- asv_status[asv_status$Reads_in_analysed_animals > 0 | asv_status$Reads_in_NEGCON1 > 0, ]
asv_status <- asv_status[order(asv_status$Status == KEPT, -asv_status$Reads_in_analysed_animals), ]

# ---- why NEGCON1: what it adds on top of the manual blocklist --------------------
rd <- function(i) sum(cnt0[ani, i])
tot_a <- sum(cnt0[ani, ])
neg_contrib <- data.frame(
  Rule = c("1 non-bacterial / chloroplast / mitochondria", "3 + 3b manual blocklist and oral/gut list (NO blank needed)",
           "4 reagent genera that are also in NEGCON1 (NEGCON1 needed)", "5 blank-dominant ASVs (NEGCON1 needed)",
           "Retained after all rules"),
  ASVs = c(sum(seen & status == STEP_LAB[1]), sum(seen & status %in% STEP_LAB[2:3]), sum(seen & status == STEP_LAB[4]),
           sum(seen & status == STEP_LAB[5]), sum(seen & keep_asv)),
  Reads_in_analysed_animals = c(rd(status == STEP_LAB[1]), rd(status %in% STEP_LAB[2:3]), rd(status == STEP_LAB[4]),
                                rd(status == STEP_LAB[5]), rd(keep_asv)))
neg_contrib$Pct_of_animal_reads <- round(100 * neg_contrib$Reads_in_analysed_animals / tot_a, 2)
neg_contrib$Needs_NEGCON1 <- c("no", "no", "YES", "YES", "")
neg_top <- asv_status[asv_status$Reads_in_NEGCON1 > 0, ]
neg_top <- neg_top[order(-neg_top$Reads_in_NEGCON1), ][1:min(40, sum(asv_status$Reads_in_NEGCON1 > 0)),
          c("ASV", "Phylum", "Family", "Genus", "Status", "Reads_in_NEGCON1", "Max_rel_pct_NEGCON1", "Max_rel_pct_animal",
            "Rel_pct_Artemia", "Reads_in_analysed_animals")]
neg2 <- NULL
if (!is.null(cl_alt)) {
  d <- which(cl_alt$status != status & !status %in% STEP_LAB[c(1, 2, 3)])
  neg2 <- data.frame(Question = c("ASVs flagged as contaminant with NEGCON1 only", "ASVs flagged with NEGCON1 + NEGCON2",
                                  "Extra ASVs flagged if NEGCON2 were also used", "Extra reads (analysed animals) that would be removed",
                                  "... as % of the retained animal reads"),
                     Value = c(sum(status %in% STEP_LAB[4:5]), sum(cl_alt$status %in% STEP_LAB[4:5]), length(d),
                               sum(cnt0[ani, d]), round(100 * sum(cnt0[ani, d]) / sum(cnt0[ani, keep_asv]), 2)))
}
why_negcon <- data.frame(Topic = c(
  "What is NEGCON1?", "What does the manual blocklist do?", "What is wrong with a blocklist alone?",
  "What does NEGCON1 add?", "Rule 4 (reagent genera)", "Rule 5 (blank-dominant ASVs)", "Why is one blank enough here?",
  "What are the limits?", "Why NEGCON2 is not used", "How to see it in this workbook"),
  Explanation = c(
  "An extraction blank: water instead of sample, carried through DNA extraction, PCR and sequencing exactly like the animals. Every bacterial sequence in it came from kit reagents, tubes, the lab or index cross-talk, not from an animal.",
  "Removes genera that are known from the literature to live on human skin, in the mouth or gut (handling contamination) - without looking at your data.",
  "It is based on reputation only. (a) Some genera on 'reagent' lists (Pseudomonas, Flavobacterium, Sphingomonas, Bacillus ...) can also be real members of a marine animal microbiome - removing them all would delete real signal. (b) It misses contaminants that nobody put on the list.",
  "Evidence from YOUR experiment about which contaminants were actually present in YOUR kit/lab. It lets us (i) remove reagent-type genera only when they really appear in the blank (rule 4), and (ii) catch contaminant ASVs that are not on any list (rule 5).",
  "A genus on the reagent list is removed only if it is also present in NEGCON1. Reagent-type genera that are NOT in the blank are kept.",
  "An ASV is removed if its relative abundance in NEGCON1 is higher than in every single animal (a sequence that is more at home in water than in any animal is contamination). Exception: if Artemia carries it even higher, it is cross-talk from the food library and is kept.",
  "We use the blank descriptively (is a taxon present / how abundant relative to the animals), which needs no statistics. Statistical tools such as decontam need many blanks to have power; with one blank the decontam column is only reported, never used.",
  "One blank cannot show how much contamination varies between extraction batches, and low-biomass animal libraries are partly made of contaminants that were also present in the blank. Results are therefore reported as 'the decontaminated bacterial community'.",
  "NEGCON2 had too few usable reads. Sheet NEGCON2_sensitivity shows how many extra ASVs it would have flagged.",
  "Sheets: NEGCON1_contribution (what rules 4 and 5 remove beyond the blocklist), NEGCON1_top_ASVs (what is in the blank), Contamination_per_sample, figures/supplementary/NEGCON1_composition."))

# blocklist coverage: which listed names exist in this dataset
cov <- function(v, what) data.frame(List = what, Name = v, Found_in_taxonomy = v %in% c(tt$Genus, tt$Family, tt$label),
  Reads_in_analysed_animals = sapply(v, function(x) sum(cnt0[ani, which(tt$Genus %in% x | tt$Family %in% x | tt$label %in% x)])),
  row.names = NULL)
coverage <- rbind(cov(MANUAL_GENERA, "Manual blocklist - genera"), cov(MANUAL_FAMILIES, "Manual blocklist - families"),
                  cov(HUMAN_ORAL_GUT, "Human oral/gut list"), cov(REAGENT, "Reagent list"))
coverage <- coverage[coverage$Found_in_taxonomy, ]
coverage <- coverage[order(coverage$List, -coverage$Reads_in_analysed_animals), ]
# what is left: top retained genera, so blocklist gaps (e.g. renamed genera) are easy to spot
left_top <- asv_status[asv_status$Status == KEPT, ] |> group_by(Phylum, Family, Genus, Habitat_list) |>
  summarise(Reads = sum(Reads_in_analysed_animals), .groups = "drop") |> arrange(desc(Reads)) |> head(40) |> as.data.frame()
names(left_top)[names(left_top) == "Habitat_list"] <- "On_habitat_list"   # "Reagent" here = on the reagent list but NOT found in NEGCON1
left_top$Pct_of_retained_reads <- round(100 * left_top$Reads / sum(cnt0[ani, keep_asv]), 2)
saveRDS(list(status = status, ani = ani, low = low), file.path(OUT, "tables", "cleaning_state.rds"))


# =============================================================================
# 3. BAR PLOTS: top families / genera BEFORE and AFTER contaminant removal
# =============================================================================
say("3. Bar plots")
lab_rank <- function(rank) { x <- tt0[[rank]]; x[is.na(x) | x == ""] <- "Unclassified"; x }
pct_by <- function(m, lab) { a <- rowsum(t(m), lab); sweep(a, 2, colSums(a), "/") * 100 }
side_levels <- c("Reference", "Reference")  # Artemia (food) and the blank share ONE side panel (two narrow panels overlap)
ref_labels  <- function(x) ifelse(x == ARTEMIA, "Artemia (food)", ifelse(x %in% BLANKS, paste(x, "(blank)"), x))
# percent table for a rank, state = "before" (all reads) or "after" (decontaminated reads)
bar_prep <- function(rank, state) {
  lab <- lab_rank(rank)
  smp <- c(ani, ARTEMIA, if (state == "before" && SHOW_BLANK_IN_BEFORE) BLANKS)
  cols <- if (state == "after") keep_asv else rep(TRUE, ncol(cnt0))
  pct <- pct_by(cnt0[smp, cols, drop = FALSE], lab[cols])
  pa  <- pct[, ani, drop = FALSE]
  top <- setdiff(names(sort(rowMeans(pa), decreasing = TRUE)), "Unclassified")[1:TOP_N]
  list(pct = pct, top = top, smp = smp, rank = rank, state = state)
}
# top taxa + one pooled "Other / Unclassified" row (everything else, incl. reads without a name at this rank)
top_plus_other <- function(b) {
  p <- b$pct; rbind(p[b$top, , drop = FALSE], colSums(p[!rownames(p) %in% b$top, , drop = FALSE])) |>
    `rownames<-`(c(b$top, OTHER_LABEL))
}
# look of the earlier run-3 bar plots: off-white background, italic title and taxon names, no bar outlines
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"
theme_bars <- theme_classic(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.background = element_rect(fill = SURF, colour = NA),
        strip.background = element_rect(fill = "#ececE6", colour = NA), strip.text = element_text(face = "bold", size = 9, colour = INK),
        strip.clip = "off", panel.spacing = unit(8, "pt"),
        axis.line = element_line(linewidth = .3, colour = "#9a9a95"), axis.text = element_text(colour = INK2),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 8), axis.title = element_text(colour = INK2),
        legend.key.size = unit(10, "pt"), legend.text = element_text(size = 8.5, hjust = 0),
        legend.background = element_rect(fill = SURF, colour = NA),
        plot.title = element_text(face = "italic", size = 13, colour = INK), plot.subtitle = element_text(size = 9, colour = INK2))
bar_plot <- function(b, cols, title = NULL, subtitle = NULL) {
  long <- as.data.frame(top_plus_other(b)) |> tibble::rownames_to_column("Taxon") |>
    pivot_longer(-Taxon, names_to = "Sample", values_to = "pct")
  long$Group <- ifelse(long$Sample == ARTEMIA, side_levels[1], ifelse(long$Sample %in% BLANKS, side_levels[2],
                       as.character(meta[long$Sample, "Group"])))
  long$Group  <- factor(long$Group, levels = unique(c(names(GROUP_COL), side_levels)))
  long$Sample <- factor(long$Sample, levels = c(ani, ARTEMIA, BLANKS))
  long$Taxon  <- factor(long$Taxon, levels = rev(c(OTHER_LABEL, b$top)))     # grey on top of each bar
  lab_it <- function(x) parse(text = ifelse(x == OTHER_LABEL, sprintf("plain('%s')", x), sprintf("italic('%s')", x)))
  ggplot(long, aes(Sample, pct, fill = Taxon)) +
    geom_col(width = .8, colour = NA) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cols, breaks = c(b$top, OTHER_LABEL), labels = lab_it, name = b$rank) +
    scale_y_continuous(breaks = seq(0, 100, 25)) + scale_x_discrete(labels = ref_labels) +
    labs(x = NULL, y = "Relative abundance (%)", title = title, subtitle = subtitle) +
    theme_bars
}
bars <- list(); bar_tabs <- list()
for (rank in c("Family", "Genus")) {
  b_bef <- bar_prep(rank, "before"); b_aft <- bar_prep(rank, "after")
  # one colour per taxon, shared between the before and after plot
  tx <- unique(c(b_aft$top, b_bef$top))      # "after" taxa get the first (run-3) colours
  cols <- c(setNames(rep_len(TAXON_PALETTE, length(tx)), tx), setNames(OTHER_COL, OTHER_LABEL))
  rk <- tolower(sub("y$", "ies", rank))
  pb <- bar_plot(b_bef, cols, sprintf("%s composition - before removal of contaminants", rank),
                 sprintf("Top %d %s by mean relative abundance; Artemia (food) and extraction blank NEGCON1 at the right", TOP_N, rk))
  pa <- bar_plot(b_aft, cols, sprintf("%s composition - after removal of contaminants", rank),
                 sprintf("Top %d %s by mean relative abundance; Artemia (food) at the right", TOP_N, rk))
  sub <- if (rank == "Family") "main" else "supplementary"
  save_fig(pb, sprintf("Bar_%s_before_removal", rank), 11, 6, sub)
  save_fig(pa, sprintf("Bar_%s_after_removal",  rank), 11, 6, sub)
  both <- (pb + labs(title = "Before removal of contaminants", subtitle = NULL)) /
          (pa + labs(title = "After removal of contaminants", subtitle = NULL)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14), plot.background = element_rect(fill = SURF, colour = NA))
  save_fig(both, sprintf("Bar_%s_before_vs_after", rank), 11, 10.5, sub)
  bars[[rank]] <- list(before = pb, after = pa)
  mk <- function(b) { d <- top_plus_other(b)
    data.frame(Taxon = rownames(d), Mean_pct_animals = round(rowMeans(d[, ani, drop = FALSE]), 2), round(d, 2), check.names = FALSE, row.names = NULL) }
  bar_tabs[[paste0("Top", TOP_N, "_", rank, "_before")]] <- mk(b_bef)
  bar_tabs[[paste0("Top", TOP_N, "_", rank, "_after")]]  <- mk(b_aft)
}

# contamination composition per library (what was removed, by rule)
cd <- as.data.frame(100 * by_status / rowSums(by_status)) |> tibble::rownames_to_column("Sample") |>
  pivot_longer(-Sample, names_to = "Status", values_to = "pct")
cd$Sample <- factor(cd$Sample, levels = c(ani, low, ARTEMIA, BLANKS))
cd$Group  <- factor(ifelse(cd$Sample == ARTEMIA, side_levels[1], ifelse(cd$Sample %in% BLANKS, side_levels[2],
                     ifelse(cd$Sample %in% low, "Removed (< min reads)", as.character(meta_all[as.character(cd$Sample), "Group"])))),
                    levels = unique(c(names(GROUP_COL), "Removed (< min reads)", side_levels)))
cd$Status <- factor(cd$Status, levels = rev(STEP_LAB))
p_contam <- ggplot(cd, aes(Sample, pct, fill = Status)) + geom_col(width = 0.88, colour = "white", linewidth = 0.12) +
  facet_grid(~ Group, scales = "free_x", space = "free_x", labeller = label_wrap_gen(11)) +
  scale_fill_manual(values = STATUS_COL, breaks = STEP_LAB, name = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) + scale_x_discrete(labels = ref_labels) +
  labs(x = NULL, y = "Share of reads (%)", title = "What the contaminant rules removed, per library") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
save_fig(p_contam, "Contamination_per_library", 11, 5.5)

# what is in the blank, and how each of its genera was classified (why NEGCON1 matters)
nb <- data.frame(ASV = colnames(cnt0), Genus = lab_rank("Genus"), Status = status, reads = cnt0[BLANKS[1], ])
nb <- nb[nb$reads > 0, ] |> group_by(Genus, Status) |> summarise(reads = sum(reads), .groups = "drop")
topg <- (nb |> group_by(Genus) |> summarise(r = sum(reads)) |> arrange(desc(r)) |> head(15))$Genus
nb <- nb[nb$Genus %in% topg, ]; nb$pct <- 100 * nb$reads / sum(cnt0[BLANKS[1], ]); nb$Genus <- factor(nb$Genus, levels = rev(topg))
p_neg <- ggplot(nb, aes(pct, Genus, fill = Status)) + geom_col(width = 0.75) +
  scale_fill_manual(values = STATUS_COL, name = "How the ASV was classified") +
  labs(x = sprintf("%% of reads in %s", BLANKS[1]), y = NULL, title = sprintf("What is in the extraction blank (%s)", BLANKS[1]),
       subtitle = "The 15 most abundant genera and the rule that removed them") +
  theme(axis.text.y = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1))
save_fig(p_neg, "NEGCON1_composition", 8.5, 5)


# =============================================================================
# 4. ALPHA DIVERSITY
# =============================================================================
say("4. Alpha diversity")
m_alpha <- counts(ps)[ani, , drop = FALSE]; depth <- min(rowSums(m_alpha))
set.seed(1); obs <- 0; sha <- 0
for (i in 1:100) { r <- rrarefy(m_alpha, depth); obs <- obs + rowSums(r > 0) / 100; sha <- sha + diversity(r, "shannon") / 100 }
alpha <- data.frame(Sample = ani, meta[ani, c("Organism", "Treatment", "Group")], Observed_ASVs = round(obs, 1), Shannon = round(sha, 3), row.names = NULL)
alpha_tests <- do.call(rbind, lapply(c("Observed_ASVs", "Shannon"), function(v) do.call(rbind, lapply(levels(meta$Organism), function(o) {
  d <- alpha[alpha$Organism == o, ]; if (min(table(d$Treatment)) < 2) return(NULL)
  w <- suppressWarnings(wilcox.test(d[[v]] ~ d$Treatment))
  data.frame(Metric = v, Host = o, n_Fed = sum(d$Treatment == "Fed"), n_Starved = sum(d$Treatment == "Starved"),
             Median_Fed = median(d[[v]][d$Treatment == "Fed"]), Median_Starved = median(d[[v]][d$Treatment == "Starved"]),
             Test = "Wilcoxon rank-sum", p_value = signif(w$p.value, 3), Significance = pstar(w$p.value)) }))))
print(alpha_tests)
p_alpha <- alpha |> pivot_longer(c(Observed_ASVs, Shannon), names_to = "Metric") |>
  mutate(Metric = recode(Metric, Observed_ASVs = "Observed ASVs", Shannon = "Shannon index")) |>
  ggplot(aes(Group, value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey25", alpha = 0.55, linewidth = 0.35) +
  geom_jitter(width = 0.1, height = 0, size = 2.4, shape = 21, colour = "grey15", stroke = 0.4) +
  facet_wrap(~ Metric, scales = "free_y") + scale_fill_manual(values = GROUP_COL, guide = "none") +
  labs(x = NULL, y = NULL, title = sprintf("Alpha diversity (rarefied to %d reads)", depth)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_alpha, "Alpha_diversity", 6.5, 3.8)


# =============================================================================
# 5. BETA DIVERSITY: PERMANOVA + PCoA
# =============================================================================
say("5. PERMANOVA and PCoA")
make_dists <- function(m) { clr <- log(m + 0.5); clr <- clr - rowMeans(clr)
  list("Bray-Curtis" = vegdist(m / rowSums(m), "bray"), "Aitchison" = dist(clr)) }
subsets <- list("All animals" = ani, "Aiptasia" = ani[meta$Organism == "Aiptasia"], "Waminoa" = ani[meta$Organism == "Waminoa"])
perm_rows <- list(); coord_rows <- list(); pcoa <- list(); pair_rows <- list()
for (sub in names(subsets)) {
  s <- subsets[[sub]]; md <- meta[s, ]; m <- counts(ps)[s, , drop = FALSE]; m <- m[, colSums(m) > 0, drop = FALSE]
  dl <- make_dists(m)
  for (dn in names(dl)) {
    d <- dl[[dn]]; set.seed(1)
    if (sub == "All animals") { a <- adonis2(d ~ Organism * Treatment, data = md, permutations = N_PERM, by = "terms"); dg <- md$Group
    } else { a <- adonis2(d ~ Treatment, data = md, permutations = N_PERM); dg <- md$Treatment }
    set.seed(1); disp_p <- permutest(betadisper(d, dg), permutations = 999)$tab$`Pr(>F)`[1]
    a <- as.data.frame(a); a <- a[!rownames(a) %in% c("Residual", "Total"), ]
    terms <- c(Organism = "Host (Aiptasia vs Waminoa)", Treatment = "Treatment (Fed vs Starved)",
               "Organism:Treatment" = "Host x Treatment interaction", Model = "Treatment (Fed vs Starved)")
    minp <- if (sub == "All animals") NA else signif(1 / choose(nrow(md), min(table(md$Treatment))), 3)
    perm_rows[[length(perm_rows) + 1]] <- data.frame(Samples = sub, n = nrow(md), Distance = dn, Effect = unname(terms[rownames(a)]),
      R2 = round(a$R2, 3), F = round(a$F, 2), p_value = a$`Pr(>F)`, Significance = pstar(a$`Pr(>F)`),
      Result = ifelse(a$`Pr(>F)` < ALPHA, "significant", "not significant"),
      Dispersion_p = round(disp_p, 3), Smallest_possible_p = minp, row.names = NULL)
    if (sub == "All animals") {                         # pairwise between the four groups (Holm-corrected)
      pr <- do.call(rbind, lapply(combn(levels(md$Group), 2, simplify = FALSE), function(g) {
        k <- md$Group %in% g; set.seed(1); aa <- adonis2(as.dist(as.matrix(d)[k, k]) ~ Group, data = md[k, ], permutations = N_PERM)
        data.frame(Distance = dn, Group_1 = g[1], Group_2 = g[2], n_1 = sum(md$Group == g[1]), n_2 = sum(md$Group == g[2]),
                   R2 = round(aa$R2[1], 3), F = round(aa$F[1], 2), p_value = aa$`Pr(>F)`[1]) }))
      pr$p_Holm <- signif(p.adjust(pr$p_value, "holm"), 3); pr$Significance <- pstar(pr$p_Holm); pair_rows[[dn]] <- pr
    }
    pc <- cmdscale(d, k = 2, eig = TRUE); ev <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
    df <- data.frame(md[, c("Organism", "Treatment", "Group")], Sample = rownames(md), PCo1 = pc$points[, 1], PCo2 = pc$points[, 2], row.names = NULL)
    coord_rows[[length(coord_rows) + 1]] <- data.frame(Samples = sub, Distance = dn, df[, c("Sample", "Organism", "Treatment")],
      PCo1 = round(df$PCo1, 4), PCo2 = round(df$PCo2, 4), PCo1_percent = round(ev[1], 1), PCo2_percent = round(ev[2], 1))
    hull <- df |> group_by(Group) |> filter(n() >= 3) |> slice(chull(PCo1, PCo2))
    ptxt <- ifelse(a$`Pr(>F)` < 0.001, "p < 0.001", sprintf("p = %.3f", a$`Pr(>F)`))
    stat <- paste(sprintf("%s: R\u00b2 = %.2f, %s", sub("\\s*\\(.*", "", unname(terms[rownames(a)])), a$R2, ptxt), collapse = "   |   ")
    p <- ggplot(df, aes(PCo1, PCo2)) +
      geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
      geom_polygon(data = hull, aes(fill = Group, colour = Group), alpha = 0.14, linewidth = 0.3, show.legend = FALSE) +
      geom_point(aes(fill = Group, shape = Organism), size = 3.4, colour = "grey10", stroke = 0.5) +
      ggrepel::geom_text_repel(aes(label = Sample), size = 2.7, colour = "grey30", seed = 1, min.segment.length = 0.3, max.overlaps = 50) +
      scale_fill_manual(values = GROUP_COL, name = NULL, drop = TRUE) + scale_colour_manual(values = GROUP_COL, guide = "none") +
      scale_shape_manual(values = c(Aiptasia = 21, Waminoa = 22), name = NULL, drop = TRUE) +
      guides(shape = "none", fill = guide_legend(override.aes = list(shape = c(21, 21, 22, 22)[names(GROUP_COL) %in% df$Group], size = 3))) +
      labs(x = sprintf("PCo1 (%.1f %%)", ev[1]), y = sprintf("PCo2 (%.1f %%)", ev[2]),
           title = sprintf("PCoA, %s distance - %s (n = %d)", dn, sub, nrow(md)), subtitle = paste("PERMANOVA  ", stat)) +
      theme(plot.subtitle = element_text(size = 7.5), legend.position = "right")
    nmx <- sprintf("PCoA_%s_%s", gsub("[^A-Za-z]", "", dn), gsub(" ", "_", sub))
    save_fig(p, nmx, if (sub == "All animals") 7.2 else 6.4, 4.8, if (sub == "All animals") "main" else "supplementary")
    pcoa[[nmx]] <- p
  }
}
perm_tab <- do.call(rbind, perm_rows); coord_tab <- do.call(rbind, coord_rows); pair_tab <- do.call(rbind, pair_rows)
print(perm_tab[, c("Samples", "Distance", "Effect", "R2", "p_value", "Significance", "Dispersion_p")])
# sensitivity: library quality (primer-dimer %) entered first
sens_tab <- NULL
if (!anyNA(meta$Primer_dimer_pct)) {
  mm <- counts(ps)[ani, ]; mm <- mm[, colSums(mm) > 0]
  sens_tab <- do.call(rbind, lapply(names(make_dists(mm)), function(dn) { d <- make_dists(mm)[[dn]]; set.seed(1)
    a <- as.data.frame(adonis2(d ~ Primer_dimer_pct + Organism * Treatment, data = meta, permutations = N_PERM, by = "terms"))
    a <- a[!rownames(a) %in% c("Residual", "Total"), ]
    data.frame(Distance = dn, n = nrow(meta), Term = c("Primer-dimer % (entered first)", "Host", "Treatment (Fed vs Starved)", "Host x Treatment")[seq_len(nrow(a))],
               R2 = round(a$R2, 3), F = round(a$F, 2), p_value = a$`Pr(>F)`, Significance = pstar(a$`Pr(>F)`)) }))
}


# =============================================================================
# 6. ANCOM-BC2 (differential abundance)
# =============================================================================
say("6. ANCOM-BC2 (ANCOMBC %s)", as.character(packageVersion("ANCOMBC")))
# neg_lb = FALSE: with 3-6 animals per group, neg_lb = TRUE would call present taxa "absent".
glom_short <- function(ps, rank) {            # collapse to a rank, short IDs (long names crash ancombc2)
  g <- tax_glom(ps, taxrank = rank, NArm = FALSE); tax <- as.data.frame(as(tax_table(g), "matrix"), stringsAsFactors = FALSE)
  lab <- tax[[rank]]; parent <- if (rank == "Genus") tax$Family else tax$Order
  na <- is.na(lab) | lab == ""; lab[na] <- paste("unclassified", ifelse(is.na(parent[na]), "Bacteria", parent[na]))
  ids <- sprintf("T%03d", seq_len(ntaxa(g))); taxa_names(g) <- ids
  list(ps = g, key = data.frame(id = ids, Taxon = lab, Family = tax$Family, Phylum = tax$Phylum)) }
fit_ancombc2 <- function(g, fixf, grp) {     # a taxon present in one group only stops ancombc2: set it aside and refit
  drop <- character(0)
  for (i in 1:10) {
    g_try <- prune_taxa(setdiff(taxa_names(g), drop), g); set.seed(123)
    fit <- tryCatch(ancombc2(data = g_try, fix_formula = fixf, group = grp, p_adj_method = "holm", prv_cut = ANCOM_PREV,
                             lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, alpha = ALPHA, pseudo_sens = TRUE, n_cl = 1, verbose = FALSE),
                    error = function(e) e)
    if (!inherits(fit, "error")) return(list(fit = fit, dropped = drop))
    bad <- intersect(unlist(strsplit(conditionMessage(fit), "[^A-Za-z0-9_]+")), taxa_names(g_try))
    if (!length(bad)) { say("  ANCOM-BC2 error: %s", conditionMessage(fit)); return(NULL) }
    drop <- union(drop, bad) }
  NULL }
run_ancom <- function(samples, fixf, grp, comparison, rank) {
  sub <- prune_samples(samples, ps); sub <- prune_taxa(taxa_sums(sub) > 0, sub)
  gs <- glom_short(sub, rank); g <- gs$ps; key <- gs$key
  lv <- levels(droplevels(meta[samples, grp])); res <- fit_ancombc2(g, fixf, grp); if (is.null(res)) return(NULL)
  r <- res$fit$res; suf <- paste0(grp, lv[2])
  rel <- counts(g); rel <- 100 * rel / rowSums(rel); grp_of <- as.character(meta[rownames(rel), grp])
  mean_in <- function(l, id) mean(rel[grp_of == l, id]); prev_in <- function(l, id) sprintf("%d of %d", sum(rel[grp_of == l, id] > 0), sum(grp_of == l))
  k <- match(r$taxon, key$id); lfc <- r[[paste0("lfc_", suf)]]; q <- r[[paste0("q_", suf)]]
  ss <- r[[paste0("passed_ss_", suf)]]; if (is.null(ss)) ss <- NA
  out <- data.frame(Comparison = comparison, Level = rank, Taxon = key$Taxon[k], Family = key$Family[k], Phylum = key$Phylum[k],
    Reference_group = lv[1], Compared_group = lv[2],
    Mean_pct_reference = round(sapply(r$taxon, function(id) mean_in(lv[1], id)), 3), Mean_pct_compared = round(sapply(r$taxon, function(id) mean_in(lv[2], id)), 3),
    Present_in_reference = sapply(r$taxon, function(id) prev_in(lv[1], id)), Present_in_compared = sapply(r$taxon, function(id) prev_in(lv[2], id)),
    Log_fold_change = round(lfc, 3), Fold_change = round(exp(lfc), 2), Higher_in = ifelse(lfc > 0, lv[2], lv[1]),
    SE = round(r[[paste0("se_", suf)]], 3), p_value = signif(r[[paste0("p_", suf)]], 3), q_value_Holm = signif(q, 3),
    Significant = ifelse(!is.na(q) & q < ALPHA, "YES", "no"), Significance = pstar(q),
    Robust_to_pseudocount = ifelse(is.na(ss), "n/a", ifelse(ss, "yes", "no")), row.names = NULL)
  out$Verdict <- with(out, ifelse(Significant == "YES" & Robust_to_pseudocount != "no", paste("DIFFERENT: higher in", Higher_in),
                           ifelse(Significant == "YES", "significant, but not robust", "no difference")))
  out <- out[order(out$p_value), ]
  zi <- res$fit$zero_ind; zeros <- NULL
  if (!is.null(zi) && nrow(zi)) {
    zc <- grep("structural_zero", names(zi), value = TRUE); zmat <- as.matrix(zi[, zc, drop = FALSE])
    n_other <- sapply(seq_len(nrow(zi)), function(i) if (sum(zmat[i, ]) == 1) sum(rel[grp_of == lv[!zmat[i, ]], zi$taxon[i]] > 0) else 0)
    hit <- rowSums(zmat) == 1 & n_other >= 2
    if (any(hit)) { kz <- match(zi$taxon[hit], key$id)
      absent <- apply(as.matrix(zi[hit, zc, drop = FALSE]), 1, function(x) paste(sub(".*= (.*)\\)$", "\\1", zc[x]), collapse = " + "))
      zeros <- data.frame(Comparison = comparison, Level = rank, Taxon = key$Taxon[kz], Family = key$Family[kz], Absent_from = absent,
                          Present_in = sapply(zi$taxon[hit], function(id) paste(lv, sapply(lv, function(l) prev_in(l, id)), collapse = "; ")), row.names = NULL)[order(-n_other[hit]), ] } }
  aside <- if (length(res$dropped)) { kd <- match(res$dropped, key$id)
    data.frame(Comparison = comparison, Level = rank, Taxon = key$Taxon[kd], Family = key$Family[kd],
               Present_in = sapply(res$dropped, function(id) paste(lv, sapply(lv, function(l) prev_in(l, id)), collapse = "; ")),
               Reason = "only in one group (zero variance): cannot be tested", row.names = NULL) }
  say("  %-44s %-6s tested %3d | significant %d | robust %d", comparison, rank, nrow(out), sum(out$Significant == "YES"), sum(grepl("^DIFFERENT", out$Verdict)))
  list(res = out, zeros = zeros, aside = aside) }
jobs <- list(
  list(s = ani[meta$Organism == "Aiptasia"], f = "Treatment", g = "Treatment", lab = "Aiptasia: Fed vs Starved"),
  list(s = ani[meta$Organism == "Waminoa"],  f = "Treatment", g = "Treatment", lab = "Waminoa: Fed vs Starved"),
  list(s = ani, f = "Organism + Treatment", g = "Organism", lab = "Aiptasia vs Waminoa (adjusted for treatment)"))
anc <- list()
for (j in jobs) for (rk in ANCOM_RANKS)
  anc[[paste(j$lab, rk)]] <- tryCatch(run_ancom(j$s, j$f, j$g, j$lab, rk), error = function(e) { say("  ! %s %s failed: %s", j$lab, rk, conditionMessage(e)); NULL })
anc_all <- do.call(rbind, lapply(anc, `[[`, "res")); anc_zero <- do.call(rbind, lapply(anc, `[[`, "zeros")); anc_aside <- do.call(rbind, lapply(anc, `[[`, "aside"))
anc_sig <- anc_all[anc_all$Significant == "YES", ]
anc_fedstarved_sig <- anc_sig[grepl("Fed vs Starved", anc_sig$Comparison), ]
if (!nrow(anc_sig)) anc_sig <- data.frame(Result = "No taxon reached q < 0.05 (Holm) in any comparison.")
anc_summary <- anc_all |> group_by(Comparison, Level) |>
  summarise(Taxa_tested = n(), Significant_q_below_0.05 = sum(Significant == "YES"), Significant_and_robust = sum(grepl("^DIFFERENT", Verdict)),
            Smallest_p = min(p_value, na.rm = TRUE), Smallest_q = min(q_value_Holm, na.rm = TRUE), Top_taxon = Taxon[which.min(p_value)], .groups = "drop") |> as.data.frame()
print(anc_summary)

# forest plot helper: the taxa with the smallest p per comparison; filled = significant (q < 0.05)
forest <- function(level, n_top, title = NULL) {
  d <- anc_all |> filter(Level == level) |> group_by(Comparison) |> slice_min(p_value, n = n_top, with_ties = FALSE) |> ungroup() |>
    mutate(key = paste(Taxon, Comparison, sep = "___"), Comparison = factor(Comparison, levels = unique(anc_all$Comparison)),
           sig = factor(ifelse(Significant == "YES", "q < 0.05 (Holm)", "not significant"), levels = c("q < 0.05 (Holm)", "not significant")))
  ggplot(d, aes(Log_fold_change, reorder(key, Log_fold_change))) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_errorbar(aes(xmin = Log_fold_change - 1.96 * SE, xmax = Log_fold_change + 1.96 * SE, colour = sig), width = 0, orientation = "y", linewidth = 0.45) +
    geom_point(aes(fill = sig, colour = sig), size = 2.3, shape = 21, stroke = 0.4) +
    scale_colour_manual(values = c("q < 0.05 (Holm)" = "#B3262D", "not significant" = "grey45"), name = NULL, drop = FALSE) +
    scale_fill_manual(values = c("q < 0.05 (Holm)" = "#B3262D", "not significant" = "white"), name = NULL, drop = FALSE) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    facet_wrap(~ Comparison, scales = "free_y", ncol = 1, labeller = label_wrap_gen(34)) +
    labs(x = "ANCOM-BC2 log fold change (95 % CI)\n> 0: higher in Starved (Fed vs Starved) or in Waminoa", y = NULL, title = title) +
    theme(axis.text.y = element_text(face = if (level == "Genus") "italic" else "plain", size = 7), legend.position = "bottom") }
p_anc_genus  <- forest("Genus", 15, "ANCOM-BC2, genus level")
p_anc_family <- forest("Family", 15, "ANCOM-BC2, family level")
save_fig(p_anc_genus,  "ANCOMBC2_genus_top15",  7.5, 10)
save_fig(p_anc_family, "ANCOMBC2_family_top15", 7.5, 10)


# =============================================================================
# 7. MAIN FIGURE: A bar plot | B PCoA | C ANCOM-BC2   (+ each panel alone)
# =============================================================================
say("7. Panel figure A/B/C")
p_C <- forest("Genus", 8)
p_A <- bars$Family$after + labs(title = NULL, subtitle = NULL) + theme(legend.position = "right")
save_fig(p_A, "Panel_A_barplot_families", 11, 4.6, "main")
save_fig(p_C, "Panel_C_ANCOMBC2_genera",  6.4, 8, "main")
for (dn in c("Bray-Curtis", "Aitchison")) {
  p_B <- pcoa[[sprintf("PCoA_%s_All_animals", gsub("[^A-Za-z]", "", dn))]] + labs(title = NULL) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 6.5))
  save_fig(p_B, sprintf("Panel_B_PCoA_%s", gsub("[^A-Za-z]", "", dn)), 6.4, 5.2, "main")
  fig <- p_A / (p_B | p_C) + plot_layout(heights = c(1, 2.2)) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 15))
  save_fig(fig, sprintf("Figure_main_ABC_%s", gsub("[^A-Za-z]", "", dn)), 13, 12.5, "main")
}


# =============================================================================
# 8. COMPARISON TABLES (previous run, SILVA 138.1 vs 138.2)
# =============================================================================
say("8. Comparison tables")
track <- tryCatch(read.csv(file.path(dirname(PS_RAW), "read_tracking.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
cmp_track <- NULL
if (!is.null(track) && file.exists(PREV_TRACK)) {
  old <- read.csv(PREV_TRACK, stringsAsFactors = FALSE)
  cmp_track <- data.frame(Sample = track$Name, Input_pairs = track$input,
    Filtered_pct_new = round(100 * track$filtered / track$input, 1), Merged_pct_new = round(100 * track$merged / track$input, 1),
    Final_reads_new = track$final, Merged_pct_previous = round(100 * old$merged[match(track$Name, old$Name)] / old$input[match(track$Name, old$Name)], 1),
    Final_reads_previous = old$no_artefact[match(track$Name, old$Name)])
  cmp_track$Final_reads_change_pct <- round(100 * (cmp_track$Final_reads_new / cmp_track$Final_reads_previous - 1), 1)
  names(cmp_track) <- sub("_new", " (this run: trunc 240/200, SILVA 138.2)", sub("_previous", " (previous: trunc 239/198)", names(cmp_track)))
}
silva_cmp <- NULL; silva_phy <- NULL; silva_gen <- NULL
f_old <- file.path(dirname(PS_RAW), "taxa_old_db.rds")
if (file.exists(f_old)) {
  told <- as.data.frame(readRDS(f_old), stringsAsFactors = FALSE)[colnames(cnt0), , drop = FALSE]
  rds_all <- colSums(cnt0)
  rate <- function(tx, rk) c(ASVs = round(100 * mean(!is.na(tx[[rk]])), 1), Reads = round(100 * sum(rds_all[!is.na(tx[[rk]])]) / sum(rds_all), 1))
  silva_cmp <- do.call(rbind, lapply(c("Phylum", "Class", "Order", "Family", "Genus", "Species"), function(rk)
    data.frame(Rank = rk, Pct_ASVs_classified_138.1 = rate(told, rk)[1], Pct_ASVs_classified_138.2 = rate(tt0, rk)[1],
               Pct_reads_classified_138.1 = rate(told, rk)[2], Pct_reads_classified_138.2 = rate(tt0, rk)[2])))
  cl_old <- get_status(told, BLANKS)
  silva_cmp <- rbind(silva_cmp, data.frame(Rank = "Contamination removed (% of reads, analysed animals)",
    Pct_ASVs_classified_138.1 = NA, Pct_ASVs_classified_138.2 = NA,
    Pct_reads_classified_138.1 = round(100 * (1 - sum(cnt0[ani, cl_old$status == KEPT]) / sum(cnt0[ani, ])), 1),
    Pct_reads_classified_138.2 = round(100 * (1 - sum(cnt0[ani, keep_asv]) / sum(cnt0[ani, ])), 1)))
  silva_phy <- as.data.frame.matrix(rowsum(rds_all, paste(ifelse(is.na(told$Phylum), "NA", told$Phylum), "->", ifelse(is.na(tt0$Phylum), "NA", tt0$Phylum))))
  silva_phy <- data.frame(Phylum_138.1_to_138.2 = rownames(silva_phy), Reads = silva_phy[, 1], row.names = NULL)
  silva_phy$Pct_reads <- round(100 * silva_phy$Reads / sum(rds_all), 2); silva_phy <- silva_phy[order(-silva_phy$Reads), ]
  ch <- is.na(told$Genus) != is.na(tt0$Genus) | (!is.na(told$Genus) & !is.na(tt0$Genus) & told$Genus != tt0$Genus)
  silva_gen <- data.frame(ASV = colnames(cnt0)[ch], Reads = rds_all[ch], Family_138.1 = told$Family[ch], Genus_138.1 = told$Genus[ch],
                          Family_138.2 = tt0$Family[ch], Genus_138.2 = tt0$Genus[ch], row.names = NULL)
  silva_gen <- head(silva_gen[order(-silva_gen$Reads), ], 300)
}


# =============================================================================
# 9. EXCEL + SUMMARY
# =============================================================================
say("9. Writing Excel")
tt_ps <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
rank_pct <- function(rank) {                          # % of reads per taxon and animal, after cleaning
  m <- counts(ps)[ani, , drop = FALSE]; x <- tt_ps[colnames(m), rank]; x[is.na(x) | x == ""] <- "Unclassified"
  p <- pct_by(m, x); p[order(-rowMeans(p)), , drop = FALSE] }
as_sheet <- function(p, rank) setNames(data.frame(rownames(p), round(rowMeans(p), 3), round(p, 3), check.names = FALSE, row.names = NULL), c(rank, "Mean % (all)", colnames(p)))
grp_means <- function(p, rank) setNames(data.frame(rownames(p), round(sapply(levels(meta$Group), function(g) rowMeans(p[, ani[meta$Group == g], drop = FALSE])), 3),
                                                   check.names = FALSE, row.names = NULL), c(rank, paste("Mean %", levels(meta$Group))))
pc_g <- rank_pct("Genus"); pc_f <- rank_pct("Family"); pc_p <- rank_pct("Phylum")
asv_tab <- data.frame(ASV = taxa_names(ps), tt_ps, t(counts(ps))[taxa_names(ps), ani], check.names = FALSE, row.names = NULL)

sheets <- c(list(
  Samples = list(samples_tab, "Every library: reads after DADA2 and after cleaning, primer-dimer %, which were analysed / removed"),
  DADA2_read_tracking = list(track, "Reads kept at each DADA2 step per library (input, filtered, denoised, merged, non-chimeric, final)"),
  Compare_previous_run = list(cmp_track, "Merging / final reads: this run (240/200, SILVA 138.2) vs previous run (239/198, SILVA 138.1)"),
  SILVA_comparison = list(silva_cmp, "SILVA 138.1 vs 138.2: % classified per rank, and effect on the contamination estimate"),
  SILVA_phylum_renaming = list(silva_phy, "Reads by phylum name in 138.1 -> 138.2 (shows renamed phyla, e.g. Proteobacteria -> Pseudomonadota)"),
  SILVA_genus_changes = list(silva_gen, "ASVs whose genus differs between the databases (top 300 by reads)"),
  Why_NEGCON1 = list(why_negcon, "Plain-language explanation: what NEGCON1 is for and what the blocklist cannot do alone"),
  NEGCON1_contribution = list(neg_contrib, "How many reads each rule removes; rules 4 and 5 are the part that NEEDS NEGCON1"),
  NEGCON1_top_ASVs = list(neg_top, "The 40 most abundant ASVs in the blank, with the rule that removed them"),
  NEGCON2_sensitivity = list(neg2, "What would change if NEGCON2 were used as a second blank"),
  Contamination_per_sample = list(contam_tab, "Reads removed per contaminant rule, contamination %, per library + totals"),
  ASV_filter_summary = list(step_tab, "ASVs and reads removed by each rule (analysed animals)"),
  Blocklist_coverage = list(coverage, "Which blocklist names exist in this dataset and how many animal reads they carry"),
  Top_retained_genera = list(left_top, "Most abundant genera LEFT after cleaning - check for contaminants missing from the blocklist"),
  ASV_all_with_status = list(asv_status, "Every ASV before cleaning with its status (why removed / retained) and per-animal counts"),
  ASV_counts_clean = list(asv_tab, "ASV counts after cleaning, with taxonomy (analysed animals)"),
  Alpha_diversity = list(alpha, "Observed ASVs and Shannon (rarefied)"), Alpha_tests = list(alpha_tests, "Fed vs Starved, per host (Wilcoxon)"),
  PERMANOVA = list(perm_tab, "Do communities differ between groups? Behind the PCoA plots"),
  PERMANOVA_pairwise = list(pair_tab, "Pairwise PERMANOVA between the 4 groups, Holm-corrected"),
  PERMANOVA_dimer_adj = list(sens_tab, "Sensitivity: PERMANOVA with primer-dimer % entered first"),
  PCoA_coordinates = list(coord_tab, "Coordinates of every point in the PCoA plots"),
  ANCOMBC2_summary = list(anc_summary, "One line per comparison and level"), ANCOMBC2_significant = list(anc_sig, "Only taxa with q < 0.05"),
  ANCOMBC2_all_results = list(anc_all, "Every tested taxon, sorted by p"),
  ANCOMBC2_one_group = list(anc_zero, "Taxa in >= 2 samples of one group and none of the other (presence/absence, no p value)"),
  ANCOMBC2_set_aside = list(anc_aside, "Taxa that could not be tested (zero variance)")),
  setNames(lapply(names(bar_tabs), function(n) list(bar_tabs[[n]], paste("% of reads per sample;", n))), names(bar_tabs)),
  list(Genus_percent = list(as_sheet(pc_g, "Genus"), "% of reads per genus and sample (after cleaning)"),
       Genus_group_means = list(grp_means(pc_g, "Genus"), "Mean % per genus in each group"),
       Family_percent = list(as_sheet(pc_f, "Family"), "% of reads per family and sample (after cleaning)"),
       Family_group_means = list(grp_means(pc_f, "Family"), "Mean % per family in each group"),
       Phylum_percent = list(as_sheet(pc_p, "Phylum"), "% of reads per phylum and sample (after cleaning)")))
sheets <- sheets[!vapply(sheets, function(s) is.null(s[[1]]) || !nrow(as.data.frame(s[[1]])), logical(1))]

readme <- rbind(
  data.frame(Item = "Created", Explanation = format(Sys.time(), "%Y-%m-%d %H:%M")),
  data.frame(Item = "Input", Explanation = PS_RAW),
  data.frame(Item = "Analysed animals", Explanation = sprintf("%d: %s", length(ani), paste(names(g_n), g_n, sep = " n=", collapse = ", "))),
  data.frame(Item = "Removed (< min reads)", Explanation = sprintf("%s (fewer than %d reads after cleaning)", if (length(low)) paste(low, collapse = ", ") else "none", MIN_READS)),
  data.frame(Item = "Blank used", Explanation = sprintf("%s (NEGCON2 left out). See sheet Why_NEGCON1.", paste(BLANKS, collapse = ", "))),
  data.frame(Item = "Contaminant rules", Explanation = "1 non-bacterial/chloroplast/mitochondria | 3 manual blocklist | 3b oral/gut genera | 4 reagent genera also in NEGCON1 | 5 blank-dominant ASVs. Edit the blocklist section in the script to change the lists."),
  data.frame(Item = "", Explanation = ""),
  data.frame(Item = paste0("Sheet: ", names(sheets)), Explanation = vapply(sheets, `[[`, "", 2)),
  data.frame(Item = "", Explanation = ""),
  data.frame(Item = c("HOW TO READ", "p_value", "q_value_Holm", "Significance", "Log_fold_change", "Fold_change", "Robust_to_pseudocount", "R2", "Dispersion_p", "Smallest_possible_p", "Bray-Curtis", "Aitchison"),
             Explanation = c("", "Chance of a difference this large if there were none.", "p corrected for many taxa (Holm). USE THIS: q < 0.05 = significant.",
               "*** p<0.001, ** p<0.01, * p<0.05, ns = not significant.", "ANCOM-BC2 effect size (natural log); > 0 = higher in the compared group (Starved / Waminoa).",
               "exp(log fold change): 2 = twice as abundant.", "'no' = result depends on the pseudo-count: treat with caution.",
               "Share of between-sample variation explained by the grouping.", "betadisper test; if < 0.05 a PERMANOVA may reflect spread, not a shift.",
               "With few animals the smallest p a permutation test can ever give (e.g. 4 vs 3 animals: 1/35 = 0.029).",
               "Distance on relative abundances.", "Compositional distance (CLR counts); same logic as ANCOM-BC2.")))
wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom", wrapText = TRUE)
add_sheet <- function(nm, df) { nm <- substr(nm, 1, 31); addWorksheet(wb, nm); writeData(wb, nm, as.data.frame(df), headerStyle = hs)
  freezePane(wb, nm, firstRow = TRUE); setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = if (ncol(df) > 40) 12 else "auto") }
add_sheet("README", readme); setColWidths(wb, "README", cols = 1:2, widths = c(30, 130))
for (nm in names(sheets)) add_sheet(nm, sheets[[nm]][[1]])
setColWidths(wb, "Why_NEGCON1", cols = 1:2, widths = c(34, 140))
addStyle(wb, "Why_NEGCON1", createStyle(wrapText = TRUE, valign = "top"), rows = 2:(nrow(why_negcon) + 1), cols = 1:2, gridExpand = TRUE)
sig_style <- createStyle(fgFill = "#E2EFDA")
for (nm in c("ANCOMBC2_all_results", "PERMANOVA", "PERMANOVA_pairwise", "Alpha_tests")) if (nm %in% names(sheets)) {
  d <- sheets[[nm]][[1]]; pcol <- intersect(c("q_value_Holm", "p_Holm", "p_value"), names(d))[1]
  rows <- which(!is.na(d[[pcol]]) & d[[pcol]] < ALPHA) + 1
  if (length(rows)) addStyle(wb, substr(nm, 1, 31), sig_style, rows = rows, cols = seq_len(ncol(d)), gridExpand = TRUE) }
xlsx <- file.path(OUT, "tables", "16S_all_tables.xlsx")
ok <- tryCatch({ saveWorkbook(wb, xlsx, overwrite = TRUE); TRUE }, error = function(e) FALSE)
if (!ok) { xlsx <- file.path(OUT, "tables", sprintf("16S_all_tables_%s.xlsx", format(Sys.time(), "%H%M%S"))); saveWorkbook(wb, xlsx, overwrite = TRUE)
  say("  the workbook was open in Excel: saved as %s", basename(xlsx)) }

# plain-text summary
hit_f <- if (nrow(anc_fedstarved_sig)) paste(anc_fedstarved_sig$Taxon, collapse = ", ") else "none"
pt <- perm_tab[perm_tab$Effect %in% c("Treatment (Fed vs Starved)"), ]
writeLines(c(
  sprintf("16S RESULTS SUMMARY  (%s)", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("Animals analysed: %d (%s). Removed for < %d reads: %s. Excluded by choice: %s.", length(ani), paste(names(g_n), g_n, sep = "=", collapse = ", "), MIN_READS, if (length(low_reads)) paste(low_reads, collapse = ", ") else "none", if (length(manual_ex)) paste(manual_ex, collapse = ", ") else "none"),
  sprintf("Contamination removed (analysed animals): %.1f %% of reads (Aiptasia %.1f %%, Waminoa %.1f %%).",
          contam_tab$Contamination_pct[contam_tab$Sample == "TOTAL analysed animals"], contam_tab$Contamination_pct[contam_tab$Sample == "TOTAL Aiptasia"], contam_tab$Contamination_pct[contam_tab$Sample == "TOTAL Waminoa"]),
  sprintf("Rules needing NEGCON1 (4 + 5) removed %.1f %% of animal reads; blocklist (3 + 3b) removed %.1f %%.", sum(neg_contrib$Pct_of_animal_reads[3:4]), neg_contrib$Pct_of_animal_reads[2]),
  "", "PERMANOVA Fed vs Starved:", sprintf("  %-10s %-12s R2 = %.3f  p = %.4f  %s  (smallest possible p: %s)", pt$Samples, pt$Distance, pt$R2, pt$p_value, pt$Significance, pt$Smallest_possible_p),
  "", sprintf("ANCOM-BC2 Fed vs Starved, taxa with q < 0.05: %s", hit_f),
  sprintf("ANCOM-BC2 Aiptasia vs Waminoa, taxa with q < 0.05: %d", sum(anc_all$Significant == "YES" & grepl("^Aiptasia vs", anc_all$Comparison))),
  "", paste("Excel:", xlsx)), file.path(OUT, "RESULTS_SUMMARY.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "R_session_info.txt"))
cat("\n"); cat(readLines(file.path(OUT, "RESULTS_SUMMARY.txt")), sep = "\n")
say("DONE. Tables: %s | Figures: %s", xlsx, file.path(OUT, "figures"))
