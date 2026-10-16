#!/usr/bin/env Rscript
# Reproduces Figure 2 (age-associated DNA methylation) from the preprocessed
# methylation data (guide script).
#
# Inputs (raw data and metadata: GEO GSE334252):
#   - beta_<date>.txt and Mval_<date>.txt, written by methylation_preprocessing.R
#   - metadata.csv (GEO sample table; columns used: idat_file_Grn, Age)
#
# Steps: limma DMPs -> power simulation -> volcano -> DMRcate DMRs -> GO on DMRs
#        -> example CpGs -> promoter DMRs -> PCA -> DMR annotation and enrichment
#        -> Figure 2 assembly.

set.seed(100)

library(vroom)
library(minfi)
library(limma)
library(DMRcate)
library(missMethyl)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(org.Hs.eg.db)
library(GenomicRanges)
library(annotatr)
library(ggtext)
library(egg)
library(grid)
library(cowplot)
library(openxlsx)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)

DATA_DIR    <- "path_to_data"
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")
OUT_DIR     <- file.path(DATA_DIR, "Figure2")
dir.create(OUT_DIR, showWarnings = FALSE)

# ===========================================================
# 1. Load preprocessed data and metadata
# ===========================================================
cat("Loading data...\n")

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

targets_full <- read.csv(META_CSV, header = TRUE)

# Slide and array position are encoded in the idat name: <Slide>_<Array>_Grn.idat
targets_full$ID    <- sub("_Grn\\.idat$", "", targets_full$idat_file_Grn)
targets_full$Slide <- sub("_.*$", "", targets_full$ID)
targets_full$Array <- sub("^[^_]*_", "", targets_full$ID)

# Sample names in the data files start with "X" (added by data.frame() in the
# preprocessing). Align metadata and beta on the Mval columns.
mval_ids     <- sub("^X", "", colnames(Mval))
targets_full <- targets_full[match(mval_ids, targets_full$ID), ]
stopifnot(!anyNA(targets_full$ID))
rownames(targets_full) <- colnames(Mval)
beta <- beta[, colnames(Mval)]

Age   <- targets_full$Age
Array <- as.factor(targets_full$Array)
Slide <- as.factor(targets_full$Slide)

# ===========================================================
# 2. Differentially methylated positions (DMPs) with age
# ===========================================================
cat("Fitting DMP models...\n")

design <- model.matrix(~ Age + Array + Slide)

fit  <- lmFit(Mval, design)
fit2 <- eBayes(fit)
DMPs <- topTable(fit2, num = Inf, coef = 2)

fitbeta  <- lmFit(beta, design)
fit2beta <- eBayes(fitbeta)
DMPbeta  <- topTable(fit2beta, num = Inf, coef = 2)
colnames(DMPbeta) <- paste0("beta_", colnames(DMPbeta))

FullDMP <- merge(DMPs, DMPbeta, by = 0)
# Delta beta per year, in %
FullDMP$beta_logFC <- 100 * FullDMP$beta_logFC
FullDMP <- FullDMP[order(FullDMP$adj.P.Val), ]

# ===========================================================
# 3. Power simulation (n = 93, genome-wide threshold)
# ===========================================================
cat("\nRunning power simulation...\n")

N_FIXED <- 93

# Residual standard deviations of the limma model fitted on beta values
empirical_sigma <- fitbeta$sigma
n_tests         <- nrow(beta)
alpha_gw        <- 0.05 / n_tests
real_ages       <- targets_full$Age[!is.na(targets_full$Age)]

cat("Number of CpGs tested:", n_tests, "\n")
cat("Genome-wide threshold (Bonferroni):", format(alpha_gw, scientific = TRUE), "\n")

slope_grid <- seq(0, 0.5, by = 0.01)

# For each slope (%/year), simulate n_perm CpGs: residual SD drawn from the
# empirical distribution, ages resampled from the observed ages, then test the
# age effect.
power_slope_fixedN <- function(n, slope_vec_pct, sigma_pool, age_pool,
                               n_perm = 100000, alpha_gw) {
  
  power_gw <- numeric(length(slope_vec_pct))
  
  for (r in seq_along(slope_vec_pct)) {
    slope <- slope_vec_pct[r] / 100
    
    pvals <- numeric(n_perm)
    for (p in seq_len(n_perm)) {
      sigma_i <- sample(sigma_pool, 1)
      age_i   <- sample(age_pool, n, replace = TRUE)
      y <- slope * age_i + rnorm(n, mean = 0, sd = sigma_i)
      mod <- lm(y ~ age_i)
      pvals[p] <- summary(mod)$coefficients["age_i", 4]
    }
    
    power_gw[r] <- mean(pvals < alpha_gw) * 100
  }
  
  data.frame(slope_pct = slope_vec_pct, power_gw = power_gw)
}

res_slope <- power_slope_fixedN(
  n = N_FIXED, slope_vec_pct = slope_grid, sigma_pool = empirical_sigma,
  age_pool = real_ages, n_perm = 1000, alpha_gw = alpha_gw
)

# Slope at which power crosses 80% (linear interpolation)
find_crossing_80 <- function(res) {
  idx <- which(res$power_gw >= 80)[1]
  if (is.na(idx) || idx == 1) return(NA)
  
  x1 <- res$slope_pct[idx - 1]; y1 <- res$power_gw[idx - 1]
  x2 <- res$slope_pct[idx];     y2 <- res$power_gw[idx]
  
  x1 + (80 - y1) * (x2 - x1) / (y2 - y1)
}

slope_80 <- find_crossing_80(res_slope)
cat("Slope at which genome-wide power reaches 80%:", round(slope_80, 4), "%/year\n")

p_power <- ggplot(res_slope, aes(x = slope_pct, y = power_gw)) +
  geom_line(color = "firebrick", linewidth = 1) +
  geom_point(size = 1.3, color = "firebrick") +
  geom_hline(yintercept = 80, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = slope_80, linetype = "dashed", color = "grey40") +
  annotate("point", x = slope_80, y = 80, size = 3, color = "black") +
  annotate("text", x = slope_80, y = 72,
           label = paste0(round(slope_80, 3), " %/year"),
           hjust = -0.1, size = 4, fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  labs(x = "Methylation slope (%/year)", y = "Power (%)",
       title = paste0("Genome-wide power vs age slope (n = ", N_FIXED, ")"),
       subtitle = paste0("Threshold: p < ", format(alpha_gw, scientific = TRUE))) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_power)

# Minimum effect size (delta beta per year, in %) used to call DMPs and DMRs
MIN_EFFECT_PCT <- 0.3
cat("Effect-size threshold used for DMPs:", MIN_EFFECT_PCT, "%/year\n")

# ===========================================================
# 4. Volcano plot
# ===========================================================
Fig1 <- FullDMP %>%
  mutate(significance = ifelse(adj.P.Val < 0.05 & abs(beta_logFC) > MIN_EFFECT_PCT,
                               "Significant", "Not significant")) %>%
  mutate(direction = ifelse(beta_logFC < 0, "Hypomethylation", "Hypermethylation")) %>%
  mutate(`Differential methylation` = ifelse(significance == "Not significant",
                                             "Not significant", direction))
Fig1$`Differential methylation` <- as.factor(Fig1$`Differential methylation`)

n_hypo  <- sum(Fig1$`Differential methylation` == "Hypomethylation")
n_hyper <- sum(Fig1$`Differential methylation` == "Hypermethylation")

FigVolcano <- ggplot(Fig1, aes(beta_logFC, -log10(P.Value))) +
  geom_point(aes(col = `Differential methylation`), size = 1) +
  scale_color_manual(values = c("Not significant" = "grey60",
                                "Hypomethylation" = "blue",
                                "Hypermethylation" = "red")) +
  labs(x = "Delta beta per year (%)", y = "-log10(p-value)") +
  annotate("text", x = min(Fig1$beta_logFC) * 0.98, y = max(-log10(Fig1$P.Value)) * 1.05,
           label = paste0("Hypomethylated CpGs: ", format(n_hypo, big.mark = ",")),
           hjust = 0, size = 4, color = "blue", fontface = "bold") +
  annotate("text", x = max(Fig1$beta_logFC) * 0.98, y = max(-log10(Fig1$P.Value)) * 1.05,
           label = paste0("Hypermethylated CpGs: ", format(n_hyper, big.mark = ",")),
           hjust = 1, size = 4, color = "red", fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, margin = margin(5, 0, 0, 0)),
    axis.text.y = element_text(size = 10, margin = margin(0, 0, 0, 0)),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

# ===========================================================
# 5. Differentially methylated regions (DMRs)
# ===========================================================
cat("\nCalling DMRs...\n")

annEPIC <- data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))

myannotation <- cpg.annotate("array", object = as.matrix(Mval), what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type = "differential",
                             design = design, coef = 2, fdr = 0.05)

dmrcoutput <- dmrcate(myannotation, lambda = 1000, C = 2, min.cpgs = 2)
DMR_ranges <- extractRanges(dmrcoutput, genome = "hg38")

data <- merge(FullDMP, annEPIC, by.y = 0, by.x = "Row.names", all.x = TRUE)

# Mean delta beta over each DMR
Annotated_ALL_GRanges <- GRanges(as.character(data$chr),
                                 IRanges(start = c(data$pos), end = c(data$pos + 1)),
                                 BetaEffect = data$beta_logFC,
                                 ind.fdr = data$adj.P.Val,
                                 is.sig = data$adj.P.Val < 0.05 & abs(data$beta_logFC) > MIN_EFFECT_PCT,
                                 probe = data$Row.names)
A <- findOverlaps(Annotated_ALL_GRanges, DMR_ranges)
DMR_beta <- merge(data.frame(A), data.frame(Annotated_ALL_GRanges),
                  by.x = "queryHits", by.y = 0)

DMR_beta_mean_per_region <- DMR_beta %>%
  group_by(subjectHits) %>%
  summarise(MeanBeta = mean(BetaEffect), MaxBeta = max(BetaEffect),
            MinBeta = min(BetaEffect),
            nbCpGneg = sum(sign(BetaEffect < 0)),
            nbCpGpos = sum(sign(BetaEffect > 0))) %>%
  mutate(ratio = ifelse(nbCpGneg > nbCpGpos, nbCpGneg / (nbCpGpos + nbCpGneg),
                        nbCpGpos / (nbCpGpos + nbCpGneg))) %>%
  mutate(sign = ifelse(ratio >= 0.75 & nbCpGneg > nbCpGpos, -1,
                       ifelse(ratio >= 0.75 & nbCpGneg < nbCpGpos, 1, 0))) %>%
  mutate(consistency = (sign == sign(MeanBeta)))

DMR_final <- merge(DMR_beta_mean_per_region, DMR_ranges,
                   by.x = "subjectHits", by.y = 0)

vroom_write(FullDMP %>% filter(adj.P.Val < 0.05 & abs(beta_logFC) > MIN_EFFECT_PCT),
            file.path(OUT_DIR, "DMP_Age.tbl"))
vroom_write(DMR_final, file.path(OUT_DIR, "DMR_Age.tbl"))

New_DMR_ranges <- GRanges(
  as.character(DMR_final$seqnames),
  IRanges(start = c(DMR_final$start), end = c(DMR_final$end)),
  MeanBeta = DMR_final$MeanBeta,
  MaxBeta = DMR_final$MaxBeta,
  MinBeta = DMR_final$MinBeta,
  nbCpGneg = DMR_final$nbCpGneg,
  nbCpGpos = DMR_final$nbCpGpos,
  sign = DMR_final$sign,
  no.cpgs = DMR_final$no.cpgs,
  min_smoothed_fdr = DMR_final$min_smoothed_fdr,
  overlapping.genes = DMR_final$overlapping.genes
)
GenomeInfoDb::genome(New_DMR_ranges) <- "hg38"

# ===========================================================
# 6. GO enrichment of DMRs (terms clustered by shared genes)
# ===========================================================
cat("\nGO enrichment of DMRs...\n")

all.cpg <- FullDMP$Row.names

gst.region <- goregion(DMR_ranges, all.cpg = all.cpg,
                       collection = "GO", array.type = "EPIC_V2", plot.bias = FALSE)
gst.region <- gst.region %>% filter(FDR < 0.05)
gst.region$GOID <- rownames(gst.region)

go2eg <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = gst.region$GOID,
                               keytype = "GOALL",
                               columns = c("ENTREZID", "GOALL")) %>%
  distinct(GOALL, ENTREZID, .keep_all = TRUE)

gene_sets <- split(go2eg$ENTREZID, go2eg$GOALL)
gene_sets <- gene_sets[gst.region$GOID]
gene_sets <- gene_sets[!sapply(gene_sets, is.null)]

# Jaccard similarity between GO terms, clustered at >= 20% shared genes
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

term_ids <- names(gene_sets)
n <- length(term_ids)
sim_matrix <- matrix(0, n, n, dimnames = list(term_ids, term_ids))

for (i in seq_len(n)) {
  for (j in i:n) {
    s <- jaccard(gene_sets[[i]], gene_sets[[j]])
    sim_matrix[i, j] <- s
    sim_matrix[j, i] <- s
  }
}

dist_mat   <- as.dist(1 - sim_matrix)
hclust_res <- hclust(dist_mat, method = "complete")
clusters   <- cutree(hclust_res, h = 1 - 0.20)

gst.region <- gst.region %>%
  mutate(Cluster = clusters[GOID]) %>%
  filter(!is.na(Cluster)) %>%
  arrange(Cluster, FDR)

cluster_sizes <- gst.region %>%
  group_by(ONTOLOGY, Cluster) %>%
  summarise(n_terms = n(), .groups = "drop")

# One representative term per cluster (best FDR), top 10 per ontology
gst_plot <- gst.region %>%
  group_by(ONTOLOGY, Cluster) %>%
  slice_min(order_by = FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(cluster_sizes, by = c("ONTOLOGY", "Cluster")) %>%
  mutate(TERM_label = paste0(TERM, " (", n_terms, ifelse(n_terms > 1, " terms)", " term)"))) %>%
  group_by(ONTOLOGY) %>%
  top_n(10, -FDR) %>%
  ungroup()

FigGO_DMR <- ggplot(gst_plot, aes(x = FDR, y = reorder(TERM_label, -FDR),
                                  color = FDR, size = DE)) +
  xlim(0, 0.05) +
  facet_wrap(~ ONTOLOGY, ncol = 1, scales = "free") +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  ylab("") +
  xlab("") +
  ggtitle("GO enrichment analysis (clusters, \u226520% shared genes)")

# ===========================================================
# 7. Example CpGs (one hyper- and one hypomethylated with age)
# ===========================================================
plot_cpg <- function(cpg, line_color) {
  fdr_val <- FullDMP$adj.P.Val[FullDMP$Row.names == cpg]
  
  df <- data.frame(
    Beta = as.numeric(beta[cpg, ]),
    Age  = targets_full$Age
  )
  
  ggplot(df, aes(x = Age, y = Beta)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = line_color) +
    labs(title = cpg,
         subtitle = paste("FDR =", format(fdr_val, scientific = TRUE, digits = 2)),
         x = "Age (years)", y = "Beta value (corrected)") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
}

FigCpG_hyper <- plot_cpg("cg12434747_TC21", "red")
FigCpG_hypo  <- plot_cpg("cg21962768_BC21", "blue")

# ===========================================================
# 8. DMRs overlapping promoters (two most significant DMRs with >= 10 CpGs)
# ===========================================================
cat("\nAnnotating promoter DMRs...\n")

annotations_promoter <- build_annotations(genome = "hg38",
                                          annotations = "hg38_genes_promoters")

dm_annotated <- annotate_regions(
  regions = New_DMR_ranges,
  annotations = annotations_promoter,
  ignore.strand = TRUE,
  quiet = FALSE
)

# Remove duplicate annotations per region
df_dm_annotated <- data.frame(dm_annotated) %>%
  group_by(seqnames, start, end, annot.type) %>%
  slice_head() %>%
  ungroup()

DMR_promoters <- df_dm_annotated %>%
  filter(!is.na(annot.symbol)) %>%
  filter(annot.type == "hg38_genes_promoters")

DMR_promoters_strong <- DMR_promoters %>% filter(no.cpgs >= 10)
top2_dmrs <- head(DMR_promoters_strong %>% arrange(min_smoothed_fdr), 2)

cpgs_granges <- GRanges(
  seqnames = data$chr,
  ranges = IRanges(start = data$pos, end = data$pos + 1),
  probe = data$Row.names,
  beta_logFC = data$beta_logFC,
  pval = data$P.Value,
  adj_pval = data$adj.P.Val
)

# Physical panel width per CpG (cm); adjust to taste
bar_slot_cm <- 0.6

FigDMR <- list()

for (i in 1:2) {
  dmr <- top2_dmrs[i, ]
  
  cat("\n=== DMR", i, "===\n")
  cat("Chromosome:", dmr$seqnames, "\n")
  cat("Position:", dmr$start, "-", dmr$end, "\n")
  cat("FDR:", format(dmr$min_smoothed_fdr, scientific = TRUE), "\n")
  
  dmr_range <- GRanges(
    seqnames = as.character(dmr$seqnames),
    ranges = IRanges(start = dmr$start, end = dmr$end)
  )
  
  overlaps <- findOverlaps(cpgs_granges, dmr_range)
  cpgs_in_dmr <- data$Row.names[queryHits(overlaps)]
  cat("Number of CpGs:", length(cpgs_in_dmr), "\n")
  
  plot_data <- data[data$Row.names %in% cpgs_in_dmr, ]
  plot_data <- plot_data %>% arrange(pos)
  plot_data$CpG_ordered <- factor(plot_data$Row.names, levels = plot_data$Row.names)
  
  n_cpgs <- nrow(plot_data)
  
  gene_name <- ifelse(!is.null(dmr$overlapping.genes) && !is.na(dmr$overlapping.genes),
                      as.character(dmr$overlapping.genes),
                      "Intergenic")
  
  # If several genes, keep the first one
  if (grepl(",", gene_name)) {
    gene_name <- strsplit(gene_name, ", ")[[1]][1]
  }
  cat("Gene:", gene_name, "\n")
  
  # Gene name in italics (markdown), except "Intergenic"
  gene_label <- ifelse(gene_name == "Intergenic", gene_name, paste0("*", gene_name, "*"))
  
  p <- ggplot(plot_data, aes(x = CpG_ordered, y = beta_logFC)) +
    geom_col(aes(fill = beta_logFC), color = "black", width = 0.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, name = "Delta beta\nper year (%)") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    labs(title = paste0(gene_label, " - ", dmr$seqnames, ":",
                        format(dmr$start, big.mark = ","), "-",
                        format(dmr$end, big.mark = ",")),
         subtitle = paste("FDR =", format(dmr$min_smoothed_fdr, scientific = TRUE, digits = 2),
                          "| n CpGs =", n_cpgs),
         x = "CpG sites (ordered by genomic position)",
         y = "Delta beta per year (%)") +
    theme_classic() +
    theme(
      text          = element_text(size = 10),
      axis.title.x  = element_text(size = 10, face = "bold"),
      axis.title.y  = element_text(size = 10, face = "bold"),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      axis.text.y   = element_text(size = 9),
      legend.text   = element_text(size = 9),
      legend.title  = element_text(size = 10, face = "bold"),
      plot.title    = ggtext::element_markdown(size = 11, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right",
      plot.margin   = margin(t = 10, r = 15, b = 10, l = 15)
    )
  
  # Panel width proportional to the number of CpGs (same bar width in both DMRs)
  FigDMR[[i]] <- set_panel_size(p,
                                width  = unit(n_cpgs * bar_slot_cm, "cm"),
                                height = unit(6, "cm"))
}

FigDMR1 <- FigDMR[[1]]
FigDMR2 <- FigDMR[[2]]

# ===========================================================
# 9. PCA (PC1 vs PC4, coloured by age)
# ===========================================================
cat("\nComputing PCA...\n")

pca_res <- prcomp(t(Mval), center = TRUE, scale. = FALSE)

n_components <- min(10, ncol(pca_res$x))
pcs <- pca_res$x[, 1:n_components]
rownames(pcs) <- colnames(Mval)

variance_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

cat("Variance explained by the first PCs:\n")
print(round(variance_explained[1:n_components], 2))

pca_df <- data.frame(
  PC1 = pcs[, "PC1"],
  PC4 = pcs[, "PC4"],
  Age = as.numeric(targets_full$Age)
)

mod_pc4 <- lm(PC4 ~ Age, data = pca_df)
r2_pc4  <- summary(mod_pc4)$r.squared
p_pc4   <- summary(mod_pc4)$coefficients[2, 4]
cat("PC4 ~ Age: R2 =", round(r2_pc4, 3), "| p =", format.pval(p_pc4, digits = 3), "\n")

FigPCA <- ggplot(pca_df, aes(x = PC1, y = PC4, color = Age)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC4 (", round(variance_explained[4], 1), "%)"),
    color = "Age"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11)
  )

# ===========================================================
# 10. DMR genomic annotation (pie chart) and CpG enrichment in DMRs
# ===========================================================
cat("\nAnnotating DMRs by genomic region...\n")

format_region <- function(x) {
  x <- gsub("hg38_genes_", "", x)
  x <- ifelse(x == "5UTRs", "5'UTRs", x)
  x <- ifelse(x == "5'UTRs", x, str_to_sentence(x))
  x
}

annots_list <- c("hg38_genes_promoters", "hg38_genes_exons", "hg38_genes_introns",
                 "hg38_genes_5UTRs", "hg38_genes_3UTRs", "hg38_genes_intergenic")
annotations_all <- build_annotations(genome = "hg38", annotations = annots_list)

dm_annotated_all <- annotate_regions(regions = New_DMR_ranges,
                                     annotations = annotations_all,
                                     ignore.strand = TRUE)
df_annot <- data.frame(dm_annotated_all)

original_ids <- data.frame(DMR_id = paste(GenomeInfoDb::seqnames(New_DMR_ranges),
                                          start(New_DMR_ranges),
                                          end(New_DMR_ranges), sep = "_"))

# One region per DMR (first annotation returned)
annotated_ids <- df_annot %>%
  mutate(DMR_id = paste(seqnames, start, end, sep = "_")) %>%
  group_by(DMR_id) %>%
  dplyr::slice(1) %>%
  mutate(Region = format_region(annot.type)) %>%
  dplyr::select(DMR_id, Region)

pie_data <- original_ids %>%
  left_join(annotated_ids, by = "DMR_id") %>%
  mutate(Region = ifelse(is.na(Region), "Not annotated", Region)) %>%
  group_by(Region) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round(count / sum(count) * 100, 1),
         legend_label = paste0(Region, " (n=", count, ", ", percentage, "%)"))

FigPie <- ggplot(pie_data, aes(x = "", y = count, fill = legend_label)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Genomic Regions") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

# Background: all tested CpGs. EPICv2 probe IDs can carry a replicate suffix
# (e.g. _TC21), removed here to match the annotation.
cpg_ids <- gsub("_TC21", "", rownames(DMPs))

ann <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
tested_cpgs <- ann[cpg_ids, ]
tested_cpgs <- tested_cpgs[!is.na(tested_cpgs$chr), ]

all_cpgs_ranges <- GRanges(
  seqnames = tested_cpgs$chr,
  ranges = IRanges(start = tested_cpgs$pos, end = tested_cpgs$pos),
  strand = tested_cpgs$strand,
  Name = rownames(tested_cpgs)
)
GenomeInfoDb::genome(all_cpgs_ranges) <- "hg38"

all_cpgs_annotated <- annotate_regions(
  regions = all_cpgs_ranges,
  annotations = annotations_all,
  ignore.strand = TRUE
)
df_all_cpgs <- data.frame(all_cpgs_annotated)

background_counts <- df_all_cpgs %>%
  mutate(Region = format_region(annot.type)) %>%
  group_by(Name) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  group_by(Region) %>%
  summarise(background_count = n(), .groups = "drop")

total_tested   <- length(all_cpgs_ranges)
annotated_cpgs <- length(unique(df_all_cpgs$Name))
not_annotated  <- total_tested - annotated_cpgs

if (not_annotated > 0) {
  background_counts <- bind_rows(
    background_counts,
    data.frame(Region = "Not annotated", background_count = not_annotated)
  )
}

# CpGs located in DMRs
overlaps <- findOverlaps(all_cpgs_ranges, New_DMR_ranges)
cpgs_in_dmrs <- all_cpgs_ranges[queryHits(overlaps)]

dmr_cpgs_annotated <- annotate_regions(
  regions = cpgs_in_dmrs,
  annotations = annotations_all,
  ignore.strand = TRUE
)
df_dmr_cpgs <- data.frame(dmr_cpgs_annotated)

dmr_counts <- df_dmr_cpgs %>%
  mutate(Region = format_region(annot.type)) %>%
  group_by(Name) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  group_by(Region) %>%
  summarise(dmr_count = n(), .groups = "drop")

total_dmr_cpgs     <- length(cpgs_in_dmrs)
annotated_dmr_cpgs <- length(unique(df_dmr_cpgs$Name))
dmr_not_annotated  <- total_dmr_cpgs - annotated_dmr_cpgs

if (dmr_not_annotated > 0) {
  dmr_counts <- bind_rows(
    dmr_counts,
    data.frame(Region = "Not annotated", dmr_count = dmr_not_annotated)
  )
}

enrichment_table <- background_counts %>%
  full_join(dmr_counts, by = "Region") %>%
  replace_na(list(background_count = 0, dmr_count = 0)) %>%
  mutate(
    background_total = sum(background_count),
    dmr_total = sum(dmr_count),
    background_pct = (background_count / background_total) * 100,
    dmr_pct = (dmr_count / dmr_total) * 100,
    fold_enrichment = dmr_pct / background_pct,
    p_value = NA_real_
  )

# Fisher's exact test for each region
for (i in 1:nrow(enrichment_table)) {
  in_region_in_dmr   <- enrichment_table$dmr_count[i]
  in_region_not_dmr  <- enrichment_table$background_count[i] - in_region_in_dmr
  not_region_in_dmr  <- enrichment_table$dmr_total[i] - in_region_in_dmr
  not_region_not_dmr <- enrichment_table$background_total[i] -
    enrichment_table$background_count[i] - not_region_in_dmr
  
  contingency <- matrix(c(in_region_in_dmr, not_region_in_dmr,
                          in_region_not_dmr, not_region_not_dmr),
                        nrow = 2, byrow = TRUE)
  
  enrichment_table$p_value[i] <- fisher.test(contingency)$p.value
}

enrichment_table <- enrichment_table %>%
  mutate(
    adj_p_value = p.adjust(p_value, method = "BH"),
    significance = case_when(
      adj_p_value < 0.001 ~ "***",
      adj_p_value < 0.01 ~ "**",
      adj_p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(desc(fold_enrichment))

print("=== ENRICHMENT OF CpGs IN DMRs ===")
print(enrichment_table %>%
        dplyr::select(Region, dmr_count, background_count, fold_enrichment,
                      adj_p_value, significance))

openxlsx::write.xlsx(enrichment_table, file.path(OUT_DIR, "CpG_enrichment_in_DMRs.xlsx"))

# Distribution of all tested CpGs vs CpGs in DMRs
comparison_data <- enrichment_table %>%
  dplyr::select(Region, background_pct, dmr_pct) %>%
  pivot_longer(cols = c(background_pct, dmr_pct),
               names_to = "Type",
               values_to = "Percentage") %>%
  mutate(Type = ifelse(Type == "background_pct",
                       "All tested CpGs",
                       "CpGs in DMRs"))

# One significance label per region, above the highest bar
sig_labels <- comparison_data %>%
  group_by(Region) %>%
  summarise(Percentage = max(Percentage), .groups = "drop") %>%
  left_join(enrichment_table %>% dplyr::select(Region, significance), by = "Region")

FigEnrichment <- ggplot(comparison_data,
                        aes(x = Region, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(data = sig_labels, aes(x = Region, y = Percentage, label = significance),
            inherit.aes = FALSE, hjust = -0.3, size = 4) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Genomic regions",
    y = "Percentage (%)",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# ===========================================================
# 11. Figure 2 assembly
# ===========================================================
cat("\nAssembling Figure 2...\n")

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

FigCpG_hyper  <- FigCpG_hyper  + theme_global
FigCpG_hypo   <- FigCpG_hypo   + theme_global
FigPCA        <- FigPCA        + theme_global
FigEnrichment <- FigEnrichment + theme_global

FigGO_DMR <- FigGO_DMR + theme_global +
  theme(
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.justification = "left",
    legend.location      = "plot",
    plot.title           = element_blank()
  )

FigVolcano <- FigVolcano +
  theme_global +
  theme(
    legend.position   = "right",
    legend.direction  = "vertical",
    legend.box        = "vertical",
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    legend.margin     = margin(t = 0, r = 0, b = 0, l = 5),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 5)
  ) +
  guides(color = guide_legend(title = "", ncol = 1, override.aes = list(size = 3)))

CpG_combined <- plot_grid(
  FigCpG_hyper, FigCpG_hypo,
  ncol = 1, nrow = 2,
  align = "v", axis = "lr"
)

Pie_Enrichment_combined <- plot_grid(
  FigPie, FigEnrichment,
  ncol = 1, nrow = 2,
  align = "v", axis = "lr"
)

# Row 1: PCA / volcano
row1 <- plot_grid(
  FigPCA, FigVolcano,
  ncol = 2,
  rel_widths = c(1, 1.5),
  align = "h", axis = "tb"
)

# Row 2: example CpGs (stacked) / promoter DMR 1 / promoter DMR 2
row2 <- plot_grid(
  CpG_combined, FigDMR1, FigDMR2,
  ncol = 3,
  rel_widths = c(1, 1, 1),
  align = "hv", axis = "tblr"
)

# Row 3: GO of DMRs / (pie chart + enrichment stacked)
row3 <- plot_grid(
  FigGO_DMR, Pie_Enrichment_combined,
  ncol = 2,
  rel_widths = c(1.3, 1),
  align = "hv", axis = "tblr"
)

combined_final <- plot_grid(
  row1, row2, row3,
  ncol = 1,
  rel_heights = c(1, 1.1, 1.1)
)

print(combined_final)

# Figure size (inches): adjust to the journal requirements
ggsave(file.path(OUT_DIR, "Figure2.pdf"), combined_final, width = 16, height = 14)

cat("\nFigure 2 written to:", file.path(OUT_DIR, "Figure2.pdf"), "\n")