#!/usr/bin/env Rscript
# Residential environment and current occupational exposures vs DNA
# methylation (guide script).
#
# Occupational exposures: participants exposed to a given agent (professionally
# exposed AND exposed to this agent) are compared with participants never
# exposed (not professionally exposed AND no exposure to this agent).
# Number of significant CpGs (adj. p < 0.05), with and without age adjustment.
#
# Inputs: beta_<date>.txt and Mval_<date>.txt written by
# methylation_preprocessing.R, and metadata.csv.
# Raw data and metadata: GEO GSE334252.

library(vroom)
library(limma)
library(stringr)
library(ggplot2)
library(dplyr)

DATA_DIR    <- "path_to_data"
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")
OUT_DIR     <- file.path(DATA_DIR, "exposure")
dir.create(OUT_DIR, showWarnings = FALSE)

# ===========================================================
# 1. Load data and metadata
# ===========================================================
# Most recent <prefix>_*.txt file written by the preprocessing script
read_methylation <- function(prefix) {
  f <- tail(sort(list.files(PREPROC_DIR, pattern = paste0("^", prefix, "_.*\\.txt$"),
                            full.names = TRUE)), 1)
  m <- as.data.frame(vroom(f, show_col_types = FALSE))
  rownames(m) <- m$CpGs
  m[, colnames(m) != "CpGs"]
}

Mval <- read_methylation("Mval")
beta <- read_methylation("beta")

targets_full <- read.csv(META_CSV, header = TRUE, na.strings = c("", "NA"))

# Slide and array position are encoded in the idat name: <Slide>_<Array>_Grn.idat
targets_full$ID    <- sub("_Grn\\.idat$", "", targets_full$idat_file_Grn)
targets_full$Slide <- sub("_.*$", "", targets_full$ID)
targets_full$Array <- sub("^[^_]*_", "", targets_full$ID)

# Sample names in the data files start with "X" (added by data.frame() in the
# preprocessing). Align metadata and beta on the Mval columns.
mval_ids     <- sub("^X", "", colnames(Mval))
targets_full <- targets_full[match(mval_ids, targets_full$ID), ]
stopifnot(!anyNA(targets_full$ID))
beta <- beta[, colnames(Mval)]

Age   <- targets_full$Age
Array <- as.factor(targets_full$Array)
Slide <- as.factor(targets_full$Slide)

# ===========================================================
# 2. Exposure variables
# ===========================================================
# Rural is the reference (control) group
residential_environment <- factor(targets_full$Residential_environment,
                                  levels = c("Rural", "Urban"))
professional_exposure   <- targets_full$ExposureActual_professional

# Detailed exposures are left empty when the participant is not exposed to the
# agent: empty or missing = "No"
na_to_no <- function(x) ifelse(is.na(x) | x == "", "No", x)

specific_exposures <- list(
  "Metals"    = na_to_no(targets_full$ExposureActual_metal),
  "Chemicals" = na_to_no(targets_full$ExposureActual_chemicals),
  "Heat"      = na_to_no(targets_full$ExposureActual_heat),
  "Smoke"     = na_to_no(targets_full$ExposureActual_smoke)
)

# ===========================================================
# 3. DMPs
# ===========================================================
# Residential environment: coefficient 2 of the design
run_residential <- function(exposure) {
  designs <- list(
    unadjusted   = model.matrix(~ exposure + Array + Slide),
    age_adjusted = model.matrix(~ exposure + Age + Array + Slide)
  )
  
  lapply(designs, function(design) {
    fit2 <- eBayes(lmFit(Mval, design))
    topTable(fit2, number = Inf, coef = 2)
  })
}

# Occupational exposure to a given agent: "Yes_Yes" vs "No_No"
# (professional exposure _ exposure to the agent)
run_occupational <- function(specific) {
  exposure_group <- paste0(professional_exposure, "_", specific)
  
  missing_groups <- setdiff(c("Yes_Yes", "No_No"), exposure_group)
  if (length(missing_groups) > 0) {
    stop("Group(s) not found: ", paste(missing_groups, collapse = ", "),
         ". Groups present: ", paste(sort(unique(exposure_group)), collapse = ", "))
  }
  
  designs <- list(
    unadjusted   = model.matrix(~ 0 + exposure_group + Array + Slide),
    age_adjusted = model.matrix(~ 0 + exposure_group + Age + Array + Slide)
  )
  
  lapply(designs, function(design) {
    fit   <- lmFit(Mval, design)
    contr <- makeContrasts("exposure_groupYes_Yes - exposure_groupNo_No",
                           levels = colnames(design))
    fit2  <- eBayes(contrasts.fit(fit, contr))
    topTable(fit2, number = Inf, coef = 1)
  })
}

results <- c(
  list("Residential environment" = run_residential(residential_environment)),
  lapply(specific_exposures, run_occupational)
)

n_sig <- function(r) {
  c(DMPs              = sum(r$unadjusted$adj.P.Val < 0.05),
    DMPs_age_adjusted = sum(r$age_adjusted$adj.P.Val < 0.05))
}

summary_table <- data.frame(Exposure = names(results),
                            t(sapply(results, n_sig)),
                            row.names = NULL)

cat("\nNumber of significant CpGs (adj. p < 0.05):\n")
print(summary_table)

