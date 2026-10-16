#!/usr/bin/env Rscript
# Interaction between age and DNA fragmentation on DNA methylation (guide
# script): DMPs for the Age x Fragmentation interaction and plots of the top 6
# CpGs.
#
# Inputs: beta_<date>.txt and Mval_<date>.txt written by
# methylation_preprocessing.R, and metadata.csv.
# Raw data and metadata: GEO GSE334252.

library(vroom)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(ggplot2)
library(ggpubr)

DATA_DIR    <- "path_to_data"
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")
OUT_DIR     <- file.path(DATA_DIR, "fragmentation")
dir.create(OUT_DIR, showWarnings = FALSE)

# Fragmentation cutoff defining the two groups shown in the plots
FRAG_CUTOFF <- 0.2

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

model_data <- data.frame(
  Age           = as.numeric(targets_full$Age),
  Fragmentation = as.numeric(targets_full$Fragmentation),
  Array         = factor(targets_full$Array),
  Slide         = factor(targets_full$Slide)
)

# Samples with missing age or fragmentation are excluded from the model
keep <- complete.cases(model_data)

# ===========================================================
# 2. Age x Fragmentation interaction model
# ===========================================================
formula_inter <- "~ Age * Fragmentation + Array + Slide"
cat("Model:", formula_inter, "\n")

design <- model.matrix(as.formula(formula_inter), data = model_data[keep, ])

inter_coef <- "Age:Fragmentation"

# M-values
fit2 <- eBayes(lmFit(Mval[, keep], design))
DMPs_inter <- topTable(fit2, coef = inter_coef, number = Inf)

# Beta values
fit2beta <- eBayes(lmFit(beta[, keep], design))
DMPbeta_inter <- topTable(fit2beta, coef = inter_coef, number = Inf)
colnames(DMPbeta_inter) <- paste0("beta_", colnames(DMPbeta_inter))

Full_Inter_Frag <- merge(DMPs_inter, DMPbeta_inter, by = 0)
Full_Inter_Frag$beta_logFC <- 100 * Full_Inter_Frag$beta_logFC
Full_Inter_Frag <- Full_Inter_Frag[order(Full_Inter_Frag$adj.P.Val), ]
colnames(Full_Inter_Frag)[1] <- "CpG"

cat("Significant CpGs (adj. p < 0.05) for the Age:Fragmentation interaction:",
    sum(Full_Inter_Frag$adj.P.Val < 0.05), "\n")

# ===========================================================
# 3. Top 6 CpGs: methylation vs age in low and high fragmentation
# ===========================================================
annEPIC <- data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))

top6_inter <- Full_Inter_Frag$CpG[1:6]

theme_global <- theme(
  text          = element_text(size = 10),
  axis.title    = element_text(size = 10, face = "bold"),
  axis.title.x  = element_text(size = 10, face = "bold"),
  axis.title.y  = element_text(size = 10, face = "bold"),
  axis.text     = element_text(size = 9),
  axis.text.x   = element_text(size = 9),
  axis.text.y   = element_text(size = 9),
  legend.text   = element_text(size = 9),
  legend.title  = element_text(size = 10, face = "bold"),
  plot.title    = element_text(size = 11, face = "bold"),
  plot.subtitle = element_text(size = 10),
  strip.text    = element_text(size = 9),
  plot.margin   = margin(t = 10, r = 15, b = 10, l = 15)
)

plot_list <- list()

for (i in 1:6) {
  cpg_id <- top6_inter[i]
  
  gene_name <- annEPIC[cpg_id, "UCSC_RefGene_Name"]
  if (is.na(gene_name) || gene_name == "") gene_name <- "Intergenic"
  
  df_plot <- data.frame(
    Age = model_data$Age,
    Beta = as.numeric(beta[cpg_id, ]),
    Frag_Group = factor(ifelse(model_data$Fragmentation > FRAG_CUTOFF, "High Frag", "Low Frag"),
                        levels = c("Low Frag", "High Frag"))
  )
  df_plot <- df_plot[complete.cases(df_plot), ]
  
  plot_list[[i]] <- ggplot(df_plot, aes(x = Age, y = Beta, color = Frag_Group)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    scale_color_manual(values = c("High Frag" = "#E41A1C", "Low Frag" = "#377EB8"),
                       labels = c("Low Frag" = "Low Fragmentation", "High Frag" = "High Fragmentation")) +
    labs(
      title = cpg_id,
      subtitle = paste(gene_name, "\nP =", formatC(Full_Inter_Frag$P.Value[i], format = "e", digits = 2)),
      x = "Age", y = "Beta value",
      color = "DNA Fragmentation"
    ) +
    theme_pubr() +
    theme_global
}

# One legend shared by the 6 panels
final_plot <- ggarrange(plotlist = plot_list,
                        ncol = 3, nrow = 2,
                        common.legend = TRUE,
                        legend = "right")

print(final_plot)

saveRDS(final_plot, file.path(OUT_DIR, "FigAgeFragmentation_interaction.rds"))

# Figure size (inches): adjust to the journal requirements
ggsave(file.path(OUT_DIR, "FigAgeFragmentation_interaction.pdf"), final_plot,
       width = 12, height = 7)