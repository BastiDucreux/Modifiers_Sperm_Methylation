#!/usr/bin/env Rscript
# BMI and DNA methylation: number of significant CpGs (adj. p < 0.05) for
# several BMI codings, with and without age adjustment (guide script).
#
# Inputs: Mval_<date>.txt written by methylation_preprocessing.R and
# metadata.csv. Raw data and metadata: GEO GSE334252.

library(vroom)
library(limma)

DATA_DIR    <- "path_to_data"
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")

# ===========================================================
# 1. Load M-values and metadata
# ===========================================================
mval_file <- tail(sort(list.files(PREPROC_DIR, pattern = "^Mval_.*\\.txt$",
                                  full.names = TRUE)), 1)
Mval <- vroom(mval_file, show_col_types = FALSE)
cpgs <- Mval$CpGs
Mval <- as.data.frame(Mval[, colnames(Mval) != "CpGs"])
rownames(Mval) <- cpgs

targets_full <- read.csv(META_CSV, header = TRUE)

# Slide and array position are encoded in the idat name: <Slide>_<Array>_Grn.idat
targets_full$ID    <- sub("_Grn\\.idat$", "", targets_full$idat_file_Grn)
targets_full$Slide <- sub("_.*$", "", targets_full$ID)
targets_full$Array <- sub("^[^_]*_", "", targets_full$ID)

# Sample names in the Mval file start with "X" (added by data.frame() in the
# preprocessing). Align metadata on the Mval columns.
mval_ids     <- sub("^X", "", colnames(Mval))
targets_full <- targets_full[match(mval_ids, targets_full$ID), ]
stopifnot(!anyNA(targets_full$ID))

Age   <- targets_full$Age
Array <- as.factor(targets_full$Array)
Slide <- as.factor(targets_full$Slide)
BMI   <- targets_full$BMI

# ===========================================================
# 2. Helper: limma fit and number of significant CpGs
# ===========================================================
# Without contrasts, tests coefficient 2 of the design (continuous BMI).
# With contrasts, tests each contrast.
run_dmp <- function(label, design, contrasts = NULL) {
  fit <- lmFit(Mval, design)
  
  if (is.null(contrasts)) {
    fit2  <- eBayes(fit)
    coefs <- 2
    names(coefs) <- "BMI"
  } else {
    contr <- makeContrasts(contrasts = contrasts, levels = colnames(design))
    fit2  <- eBayes(contrasts.fit(fit, contr))
    coefs <- seq_along(contrasts)
    names(coefs) <- contrasts
  }
  
  res <- lapply(coefs, function(k) topTable(fit2, number = Inf, coef = k))
  
  cat("\n", label, "\n", sep = "")
  print(sapply(res, function(x) sum(x$adj.P.Val < 0.05)))
  
  invisible(res)
}

# ===========================================================
# 3. Four BMI categories (<18.5, 18.5-25, 25-30, >30; reference 18.5-25)
# ===========================================================
BMI_cat <- ifelse(BMI < 18.5, "under18.5",
                  ifelse(BMI <= 25, "18.5_25",
                         ifelse(BMI <= 30, "25_30", "over30")))

contr_cat <- c("BMI_cat25_30 - BMI_cat18.5_25",
               "BMI_catover30 - BMI_cat18.5_25")

dmp_cat4     <- run_dmp("4 categories",
                        model.matrix(~ 0 + BMI_cat + Array + Slide), contr_cat)
dmp_cat4_age <- run_dmp("4 categories, adjusted for age",
                        model.matrix(~ 0 + BMI_cat + Age + Array + Slide), contr_cat)

# ===========================================================
# 4. Two BMI categories (normal < 25 vs overweight >= 25)
# ===========================================================
BMI_cat <- ifelse(BMI < 25, "normal", "overweight")

contr_bin <- "BMI_catoverweight - BMI_catnormal"

dmp_cat2     <- run_dmp("2 categories",
                        model.matrix(~ 0 + BMI_cat + Array + Slide), contr_bin)
dmp_cat2_age <- run_dmp("2 categories, adjusted for age",
                        model.matrix(~ 0 + BMI_cat + Age + Array + Slide), contr_bin)

# ===========================================================
# 5. Continuous BMI
# ===========================================================
dmp_cont     <- run_dmp("Continuous BMI",
                        model.matrix(~ BMI + Array + Slide))
dmp_cont_age <- run_dmp("Continuous BMI, adjusted for age",
                        model.matrix(~ BMI + Age + Array + Slide))

# ===========================================================
# 6. Four BMI categories with a 20 cutoff (<20, 20-25, 25-30, >30; reference 20-25)
# ===========================================================
BMI_cat <- ifelse(BMI < 20, "under20",
                  ifelse(BMI <= 25, "20_25",
                         ifelse(BMI <= 30, "25_30", "over30")))

contr_cat20 <- c("BMI_cat25_30 - BMI_cat20_25",
                 "BMI_catover30 - BMI_cat20_25")

dmp_cat20     <- run_dmp("4 categories (20 cutoff)",
                         model.matrix(~ 0 + BMI_cat + Array + Slide), contr_cat20)
dmp_cat20_age <- run_dmp("4 categories (20 cutoff), adjusted for age",
                         model.matrix(~ 0 + BMI_cat + Age + Array + Slide), contr_cat20)