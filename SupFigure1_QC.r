#!/usr/bin/env Rscript
# Reproduces Supplementary Figure 1 (guide script):
#   A. Sample QC: median methylated vs unmethylated intensity
#   B. Beta-value density of each sample
#   C. Trimodal CpGs (candidate genotype-driven CpGs), top 500 by clustering score
#
# Inputs (raw data and metadata: GEO GSE334252):
#   - GSE334252_RAW/ (idat files) and metadata.csv, for panel A
#   - beta_<date>.txt written by methylation_preprocessing.R, for panels B and C

library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(vroom)
library(matrixStats)
library(mclust)
library(ggplot2)
library(patchwork)

DATA_DIR    <- "path_to_data"
IDAT_DIR    <- file.path(DATA_DIR, "GSE334252_RAW")
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")
OUT_DIR     <- file.path(DATA_DIR, "FigureS1")
dir.create(OUT_DIR, showWarnings = FALSE)

# ===========================================================
# A. Sample QC (raw intensities)
# ===========================================================
cat("Reading idat files...\n")

targets <- read.csv(META_CSV, header = TRUE, na.strings = c("", "NA"))

# Slide and array position are encoded in the idat name: <Slide>_<Array>_Grn.idat
targets$ID <- sub("_Grn\\.idat$", "", targets$idat_file_Grn)

# Locate the idat files; GEO may prefix names with the GSM accession and gzip them
grn_files <- list.files(IDAT_DIR, pattern = "_Grn\\.idat(\\.gz)?$",
                        full.names = TRUE, recursive = TRUE)
grn_ids   <- sub("_Grn\\.idat(\\.gz)?$", "", sub("^GSM[0-9]+_", "", basename(grn_files)))
targets$Basename <- sub("_Grn\\.idat(\\.gz)?$", "", grn_files[match(targets$ID, grn_ids)])
stopifnot(!anyNA(targets$Basename))

rgSet <- read.metharray.exp(targets = targets)
rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

MSet <- preprocessRaw(rgSet)
qc   <- getQC(MSet)

qc_data <- data.frame(
  mMed = qc$mMed,
  uMed = qc$uMed,
  sample_index = 1:nrow(qc)
)

# Cutoff line: mean of the two median intensities = 10.5, i.e. mMed + uMed = 21
x_seq  <- seq(8, 13, length.out = 100)
y_line <- 21 - x_seq

FigQC <- ggplot(qc_data, aes(x = mMed, y = uMed)) +
  # Red area: below the line
  geom_ribbon(data = data.frame(x = x_seq, ymin = 8, ymax = y_line),
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#FFB3B3", alpha = 0.3, inherit.aes = FALSE) +
  # Green area: above the line
  geom_ribbon(data = data.frame(x = x_seq, ymin = y_line, ymax = 13),
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#B3FFB3", alpha = 0.3, inherit.aes = FALSE) +
  geom_abline(intercept = 21, slope = -1,
              linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(shape = 1, color = "black", size = 2) +
  annotate("text", x = 10.5, y = 9.5,
           label = "Bad intensity samples",
           color = "#CC0000", size = 3.5, fontface = "italic") +
  annotate("text", x = 10.5, y = 11.5,
           label = "Good intensity samples",
           color = "#007700", size = 3.5, fontface = "italic") +
  scale_x_continuous(breaks = c(9, 11, 13), limits = c(8, 13)) +
  scale_y_continuous(breaks = c(9, 11, 13), limits = c(8, 13)) +
  labs(x = "Meth median intensity (log2)",
       y = "Unmeth median intensity (log2)") +
  theme_classic()

rm(rgSet, MSet)

# ===========================================================
# Beta values (preprocessed)
# ===========================================================
# Most recent beta_*.txt file written by the preprocessing script
beta_file <- tail(sort(list.files(PREPROC_DIR, pattern = "^beta_.*\\.txt$",
                                  full.names = TRUE)), 1)
beta_df <- as.data.frame(vroom(beta_file, show_col_types = FALSE))
rownames(beta_df) <- beta_df$CpGs
beta_mat <- as.matrix(beta_df[, colnames(beta_df) != "CpGs"])
rm(beta_df)

# ===========================================================
# B. Beta-value density of each sample
# ===========================================================
density_df <- do.call(rbind, lapply(seq_len(ncol(beta_mat)), function(i) {
  d <- density(beta_mat[, i], na.rm = TRUE, n = 512, from = -0.05, to = 1.05)
  data.frame(x = d$x, y = d$y, sample = colnames(beta_mat)[i])
}))

FigDensity <- ggplot(density_df, aes(x = x, y = y, group = sample)) +
  geom_line(alpha = 0.5) +
  xlim(-0.05, 1.05) +
  xlab("Beta-value") +
  ylab("Density") +
  theme_classic()

# ===========================================================
# C. Trimodal CpGs
# ===========================================================
cat("Searching trimodal CpGs (this step is long)...\n")

# Most variable CpGs (top 5% of the variance)
vars     <- rowVars(beta_mat)
names(vars) <- rownames(beta_mat)
high_var <- names(vars[vars > quantile(vars, 0.95)])

# Number of mixture components (1 to 3) for each CpG
mclust_res <- sapply(high_var, function(cpg) {
  x <- na.omit(as.numeric(beta_mat[cpg, ]))
  tryCatch({
    fit <- Mclust(x, G = 1:3, verbose = FALSE)
    fit$G
  }, error = function(e) NA)
})

# CpGs with 3 components are trimodal
trimodal_cpgs <- names(mclust_res[!is.na(mclust_res) & mclust_res == 3])
cat("Trimodal CpGs detected:", length(trimodal_cpgs), "\n")

# Clustering score: 1 - (within-cluster variance / total variance)
mclust_score <- sapply(trimodal_cpgs, function(cpg) {
  x <- na.omit(as.numeric(beta_mat[cpg, ]))
  tryCatch({
    fit            <- Mclust(x, G = 3, verbose = FALSE)
    cluster_assign <- apply(fit$z, 1, which.max)
    var_intra      <- mean(tapply(x, cluster_assign, var, na.rm = TRUE), na.rm = TRUE)
    var_total      <- var(x, na.rm = TRUE)
    1 - (var_intra / var_total)
  }, error = function(e) NA)
})

# Best-separated clusters first
trimodal_cpgs_top <- head(names(sort(mclust_score, decreasing = TRUE)), 1000)

# Top 500 CpGs for the figure
top_cpgs <- head(trimodal_cpgs_top, 500)

data_plot <- do.call(rbind, lapply(seq_along(top_cpgs), function(i) {
  x <- as.numeric(beta_mat[top_cpgs[i], ])
  data.frame(CpG_index = i, CpG_name = top_cpgs[i], beta = x[!is.na(x)])
}))

FigMeQTL <- ggplot(data_plot, aes(x = CpG_index, y = beta)) +
  geom_point(alpha = 0.5, size = 1, color = "black") +
  labs(x = "CpG (ranked by trimodal score)",
       y = "Methylation level (Beta value)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ===========================================================
# Figure assembly
# ===========================================================
# Same style for all panels
force_style <- function(p) {
  p + theme_classic() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black")
    )
}

pA <- force_style(FigQC)
pB <- force_style(FigDensity)
pC <- force_style(FigMeQTL)

# | = side by side, / = stacked
combined_final <- (pA | pB) / (pC | plot_spacer()) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

print(combined_final)

saveRDS(combined_final, file.path(OUT_DIR, "SupplementaryFigure1.rds"))

# Figure size (inches): adjust to the journal requirements
ggsave(file.path(OUT_DIR, "SupplementaryFigure1.pdf"), combined_final,
       width = 12, height = 10)