#!/usr/bin/env Rscript
# SVD of M-values and association of the principal components with covariates
# (guide script).
#
# Requires the output of methylation_preprocessing.R (Mval_<date>.txt) and
# metadata.csv. Raw data and metadata: GEO GSE334252.

library(vroom)
library(ggplot2)
library(reshape2)

DATA_DIR    <- "path_to_data"
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")
OUT_DIR     <- file.path(DATA_DIR, "SVD")
dir.create(OUT_DIR, showWarnings = FALSE)

# ===========================================================
# 1. Load preprocessed data
# ===========================================================
cat("Loading preprocessed data...\n")

# Most recent Mval file written by the preprocessing script
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
rownames(targets_full) <- colnames(Mval)

# ===========================================================
# 2. Covariates
# ===========================================================
pd <- data.frame(
  `Array batch`             = as.factor(targets_full$Slide),
  `Array position`          = as.factor(targets_full$Array),
  `Age (years)`             = as.numeric(targets_full$Age),
  `BMI (kg/m2)`             = as.numeric(targets_full$BMI),
  `Residential environment` = as.factor(targets_full$Residential_environment),
  `Professional exposure`   = as.factor(targets_full$ExpositionActual_or_Previous_professional),
  `Smoking (active)`        = as.factor(targets_full$Smoking),
  `Smoking (pack-year)`     = as.numeric(ifelse(is.na(targets_full$Pack_year), 0,
                                                targets_full$Pack_year)),
  `HbCO (%)`                = as.numeric(targets_full$HbCO),
  `DNA Fragmentation`       = as.numeric(targets_full$Fragmentation),
  `Sperm parameters`        = as.factor(targets_full$Sperm_param_A_astheno_T_terato_O_oligo_NI_not_identified),
  `Early pregnancy`         = as.factor(targets_full$Early_pregnancy),
  row.names = rownames(targets_full), check.names = FALSE
)

# Set factor levels with fewer than 2 samples to NA
for (col in colnames(pd)) {
  if (is.factor(pd[[col]])) {
    tab <- table(pd[[col]])
    pd[[col]][pd[[col]] %in% names(tab[tab < 2])] <- NA
    pd[[col]] <- droplevels(pd[[col]])
  }
}

# ===========================================================
# 3. SVD
# ===========================================================
cat("\nComputing SVD...\n")

# Center each CpG (no scaling)
Mval_scaled <- t(scale(t(Mval), center = TRUE, scale = FALSE))

svd_res <- svd(Mval_scaled)

n_components <- min(10, ncol(svd_res$v))
pcs <- svd_res$v[, 1:n_components]
colnames(pcs) <- paste0("PC", 1:n_components)
rownames(pcs) <- colnames(Mval)

variance_explained <- (svd_res$d^2 / sum(svd_res$d^2)) * 100

cat("Variance explained by the first PCs:\n")
print(round(variance_explained[1:n_components], 2))

# ===========================================================
# 4. PC ~ covariate association tests
# ===========================================================
cat("\nTesting PC ~ covariate associations...\n")

results <- data.frame()

for (i in 1:n_components) {
  pc_name <- paste0("PC", i)
  
  for (cov in colnames(pd)) {
    valid_idx <- !is.na(pd[[cov]])
    pc_vals   <- pcs[valid_idx, i]
    cov_vals  <- pd[[cov]][valid_idx]
    
    if (length(unique(cov_vals)) < 2) next
    
    if (is.numeric(cov_vals)) {
      # Linear regression for continuous covariates
      test_res <- summary(lm(pc_vals ~ cov_vals))
      pval <- test_res$coefficients[2, 4]
      rsq  <- test_res$r.squared
    } else {
      # ANOVA for categorical covariates
      test_res <- summary(aov(pc_vals ~ cov_vals))
      pval     <- test_res[[1]][["Pr(>F)"]][1]
      ss_total <- sum((pc_vals - mean(pc_vals))^2)
      ss_resid <- sum(test_res[[1]][["Sum Sq"]][-1])
      rsq      <- 1 - (ss_resid / ss_total)
    }
    
    results <- rbind(results, data.frame(
      PC        = pc_name,
      Covariate = cov,
      Pvalue    = pval,
      R2        = rsq
    ))
  }
}

results$FDR <- p.adjust(results$Pvalue, method = "BH")

# ===========================================================
# 5. Heatmap
# ===========================================================
cat("\nCreating plot...\n")

results$Covariate <- gsub("\\.", " ", results$Covariate)

heatmap_data <- dcast(results, Covariate ~ PC, value.var = "Pvalue")
rownames(heatmap_data) <- heatmap_data$Covariate
heatmap_data$Covariate <- NULL

# Order covariates by the first PC where they reach p < 0.05
# (ties broken by p-value; never significant -> last)
covariate_first_sig <- sapply(rownames(heatmap_data), function(cov) {
  pvals     <- as.numeric(heatmap_data[cov, ])
  first_sig <- which(pvals < 0.05)[1]
  if (is.na(first_sig)) {
    return(999)
  }
  return(first_sig + pvals[first_sig] / 1000)
})

covariate_order <- names(sort(covariate_first_sig))

heatmap_long <- melt(as.matrix(heatmap_data))
colnames(heatmap_long) <- c("Covariate", "PC", "Pvalue")
heatmap_long$Covariate <- as.character(heatmap_long$Covariate)

heatmap_long$PC <- factor(heatmap_long$PC,
                          levels = paste0("PC", 1:n_components))
heatmap_long$Covariate <- factor(heatmap_long$Covariate,
                                 levels = rev(covariate_order))

heatmap_long$Pvalue_cat <- cut(heatmap_long$Pvalue,
                               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                               labels = c("p < 0.001", "0.001 \u2264 p < 0.01",
                                          "0.01 \u2264 p < 0.05", "p \u2265 0.05"),
                               right = FALSE)

# x-axis labels with the percentage of variance explained
pc_labels <- paste0("PC", 1:n_components, "\n(",
                    round(variance_explained[1:n_components], 1), "%)")
names(pc_labels) <- paste0("PC", 1:n_components)

color_palette <- c("p < 0.001"           = "blue",
                   "0.001 \u2264 p < 0.01" = "purple",
                   "0.01 \u2264 p < 0.05"  = "red",
                   "p \u2265 0.05"         = "#E8E8E8")

SVD_plot <- ggplot(heatmap_long, aes(x = PC, y = Covariate, fill = Pvalue_cat)) +
  geom_tile(color = "grey50", linewidth = 0.5) +
  scale_fill_manual(
    values = color_palette,
    name = "Significance",
    drop = FALSE
  ) +
  scale_x_discrete(labels = pc_labels) +
  geom_hline(yintercept = seq(0.5, length(unique(heatmap_long$Covariate)) + 0.5, 1),
             color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = seq(0.5, n_components + 0.5, 1),
             color = "grey80", linewidth = 0.3) +
  labs(
    x = "Principal Component",
    y = "Covariate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank()
  )

# ===========================================================
# 6. Global figure theme and export
# ===========================================================
theme_global <- theme(
  text          = element_text(size = 10),
  axis.title    = element_text(size = 10, face = "bold"),
  axis.text     = element_text(size = 9),
  legend.text   = element_text(size = 9),
  legend.title  = element_text(size = 10, face = "bold"),
  plot.title    = element_text(size = 11, face = "bold"),
  plot.subtitle = element_text(size = 10)
)

FigSVD <- SVD_plot + theme_global
print(FigSVD)

saveRDS(FigSVD, file.path(OUT_DIR, "FigSVD.rds"))

cat("\n==========================================\n")
cat("SVD COMPLETE\n")
cat("Top significant associations:\n")
print(head(results[order(results$Pvalue), ], 20))
cat("==========================================\n")