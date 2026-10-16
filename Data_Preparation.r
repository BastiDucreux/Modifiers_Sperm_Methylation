#!/usr/bin/env Rscript
# EPICv2 methylation preprocessing pipeline (guide script)
#
# Raw data (idat files) and metadata are available on GEO: GSE334252.
# This script is provided as a guide to the preprocessing steps and is not
# meant to run as-is: paths must be adapted.
#
# Expected content of path_to_data:
#   - GSE334252_RAW/        idat files (GEO RAW archive, extracted)
#   - metadata.csv          GEO sample table exported to CSV; columns used
#                           here: idat_file_Grn, Age, BMI
#   - EPICv2_ReactiveLoci.txt   cross-reactive / reactive probes list
# The last step builds the R environment used for downstream analyses.
#
# Sample QC: exclude samples with >1% of CpGs failing detection p-value.
# Probe filter: NA-tolerant detP filter (a single NA does not drop the CpG).

set.seed(100)
library(minfi)
library(minfiData)
library(wateRmelon)
library(maxprobes)
library(dplyr)
library(tidyverse)
library(IlluminaHumanMethylation27kmanifest)
library(IlluminaHumanMethylation27kanno.ilmn12.hg19)
library(MASS)
library(readxl)
library(ChAMP)
library(ggplot2)
library(limma)
library(vroom)
library(openxlsx)
library(stringr)
library(FactoMineR)
library(PCAtools)
library(DMRcate)
library(missMethyl)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

date <- format(Sys.Date(), "%Y%m%d")

DATA_DIR   <- "path_to_data"
OUT_DIR    <- file.path(DATA_DIR, "output")
QC_DIR     <- file.path(OUT_DIR, "QC")
XLOCI_FILE <- file.path(DATA_DIR, "EPICv2_ReactiveLoci.txt")
META_CSV   <- file.path(DATA_DIR, "metadata.csv")
IDAT_DIR   <- file.path(DATA_DIR, "GSE334252_RAW")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(QC_DIR,  showWarnings = FALSE)

# ===========================================================
# STEP 1: Metadata
# ===========================================================
cat("==========================================\n")
cat("=== STEP 1: Loading metadata ===\n")
cat("==========================================\n")

targets <- read.csv(META_CSV, header = TRUE)

# Slide and array position are encoded in the idat name: <Slide>_<Array>_Grn.idat
targets$ID       <- sub("_Grn\\.idat$", "", targets$idat_file_Grn)

# Locate the idat files; GEO may prefix names with the GSM accession and gzip them
grn_files <- list.files(IDAT_DIR, pattern = "_Grn\\.idat(\\.gz)?$",
                        full.names = TRUE, recursive = TRUE)
grn_ids   <- sub("_Grn\\.idat(\\.gz)?$", "", sub("^GSM[0-9]+_", "", basename(grn_files)))
targets$Basename <- sub("_Grn\\.idat(\\.gz)?$", "", grn_files[match(targets$ID, grn_ids)])
stopifnot(!anyNA(targets$Basename))

targets$Slide    <- sub("_.*$", "", targets$ID)
targets$Array    <- sub("^[^_]*_", "", targets$ID)
cat("targets dimensions:", nrow(targets), "x", ncol(targets), "\n")

# ===========================================================
# STEP 2: Reading idat files
# ===========================================================
cat("==========================================\n")
cat("=== STEP 2: Reading idat files ===\n")
cat("==========================================\n")

rgSet <- read.metharray.exp(targets = targets)
rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

colnames(rgSet) <- targets$ID
cat("rgSet dimensions:", dim(rgSet), "\n")

write.table(targets, file.path(OUT_DIR, paste0("metadata_", date, ".txt")),
            sep = " ", quote = FALSE, row.names = FALSE)

# ===========================================================
# STEP 3: QC on raw data
# ===========================================================
cat("==========================================\n")
cat("=== STEP 3: QC ===\n")
cat("==========================================\n")

MSet <- preprocessRaw(rgSet)
qc   <- getQC(MSet)

tiff(file.path(QC_DIR, "QC.tif"), res = 300, width = 1400, height = 1200)
plotQC(qc)
dev.off()
cat("QC plot OK\n")

# 10.5 is the median intensity cutoff used by minfi's plotQC
meds        <- (qc$mMed + qc$uMed) / 2
bad_samples <- rownames(qc)[which(meds < 10.5)]
cat("Low-quality samples (median < 10.5):",
    ifelse(length(bad_samples) == 0, "none", paste(bad_samples, collapse = ", ")), "\n")

# ===========================================================
# STEP 3b: Sample QC - exclude if >1% of CpGs fail detP
# ===========================================================
cat("==========================================\n")
cat("=== STEP 3b: Sample QC (1% detP threshold) ===\n")
cat("==========================================\n")

detP_qc <- detectionP(rgSet)
pct_fail_per_sample <- 100 * colMeans(detP_qc >= 0.01, na.rm = TRUE)

bad_samples_detP <- colnames(rgSet)[pct_fail_per_sample > 1]
cat("Samples excluded (>1% CpGs failing detP):", length(bad_samples_detP), "\n")
if (length(bad_samples_detP) > 0) {
  for (s in bad_samples_detP) {
    cat("  -", s, ":", round(pct_fail_per_sample[s], 2), "% failing\n")
  }
}

keep_samples_detP <- !(colnames(rgSet) %in% bad_samples_detP)
rgSet   <- rgSet[, keep_samples_detP]
targets <- targets[keep_samples_detP, ]

cat("Number of samples after detP QC:", ncol(rgSet), "\n")

# ===========================================================
# STEP 4: Noob + BMIQ normalization
# ===========================================================
cat("==========================================\n")
cat("=== STEP 4: Noob + BMIQ normalization ===\n")
cat("==========================================\n")

mSet     <- preprocessNoob(rgSet)
BetaBMIQ <- wateRmelon::BMIQ(mSet)
cat("BetaBMIQ dimensions:", dim(BetaBMIQ), "\n")

# ===========================================================
# STEP 5: detP probe filtering
# ===========================================================
cat("==========================================\n")
cat("=== STEP 5: detP filtering ===\n")
cat("==========================================\n")

# Keep a CpG if every non-NA detP value is < 0.01
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01, na.rm = TRUE) == rowSums(!is.na(detP))
cat("CpGs kept after detP filter:", sum(keep), "/", nrow(detP), "\n")
cat("CpGs with at least one NA in detP:", sum(rowSums(is.na(detP)) > 0), "\n")
mSetQC <- mSet[keep, ]

# ===========================================================
# STEP 6: Filter SNPs, cross-reactive and non-cg probes
# ===========================================================
cat("==========================================\n")
cat("=== STEP 6: Filtering SNPs, cross-reactive and non-cg probes ===\n")
cat("==========================================\n")

t           <- mapToGenome(mSetQC)
mSetQCnosnp <- dropLociWithSnps(t)

xloci     <- read.table(XLOCI_FILE)
xloci     <- xloci[, 1]
mSetFinal <- mSetQCnosnp[!(rownames(mSetQCnosnp) %in% xloci), ]

mSetFinal <- mSetFinal[grepl("^cg", rownames(mSetFinal)), ]
cat("CpGs after SNP, cross-reactive and non-cg filtering:", nrow(mSetFinal), "\n")

# ===========================================================
# STEP 7: Beta and M-values
# ===========================================================
cat("==========================================\n")
cat("=== STEP 7: Computing beta and M-values ===\n")
cat("==========================================\n")

beta <- BetaBMIQ[rownames(mSetFinal), ]
Mval <- log2(beta / (1 - beta))

beta <- data.frame(beta)
Mval <- data.frame(Mval)
beta$CpGs <- rownames(beta)
Mval$CpGs <- rownames(Mval)

# ===========================================================
# STEP 8: Export
# ===========================================================
cat("==========================================\n")
cat("=== STEP 8: Exporting files ===\n")
cat("==========================================\n")

vroom_write(beta, file.path(OUT_DIR, paste0("beta_", date, ".txt")))
vroom_write(Mval, file.path(OUT_DIR, paste0("Mval_", date, ".txt")))

# ===========================================================
# STEP 9: Analysis environment
# ===========================================================
cat("==========================================\n")
cat("=== STEP 9: Building analysis environment ===\n")
cat("==========================================\n")

beta <- vroom(file.path(OUT_DIR, paste0("beta_", date, ".txt")))
Mval <- vroom(file.path(OUT_DIR, paste0("Mval_", date, ".txt")))

# targets was already filtered by the detP sample QC (step 3b)
targets_full <- targets

# data.frame() prefixes sample names starting with a digit with "X"
targets_full$ID <- paste0("X", targets_full$ID)
rownames(targets_full) <- targets_full$ID

beta <- data.frame(beta)
Mval <- data.frame(Mval)
rownames(beta) <- beta$CpGs
rownames(Mval) <- Mval$CpGs

# Same sample order as the metadata
beta <- beta[, targets_full$ID]
Mval <- Mval[, targets_full$ID]

# Technical covariates (batch correction)
Slide <- as.factor(targets_full$Slide)
Array <- as.factor(targets_full$Array)

# Physical covariates
Age <- targets_full$Age
BMI <- targets_full$BMI

save.image(file.path(OUT_DIR, "env.RData"))
cat("Environment saved:", file.path(OUT_DIR, "env.RData"), "\n")

cat("==========================================\n")
cat("=== PIPELINE COMPLETE ===\n")
cat("Samples excluded by detP QC:", length(bad_samples_detP), "\n")
cat("beta     :", file.path(OUT_DIR, paste0("beta_", date, ".txt")), "\n")
cat("Mval     :", file.path(OUT_DIR, paste0("Mval_", date, ".txt")), "\n")
cat("metadata :", file.path(OUT_DIR, paste0("metadata_", date, ".txt")), "\n")
cat("QC       :", QC_DIR, "\n")
cat("==========================================\n")