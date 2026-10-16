#!/usr/bin/env Rscript
# Sperm parameters and DNA methylation (guide script): sperm parameter class
# and DNA fragmentation, with and without age adjustment; DMPs and DMRs
# associated with DNA fragmentation.
#
# Inputs: beta_<date>.txt and Mval_<date>.txt written by
# methylation_preprocessing.R, and metadata.csv.
# Raw data and metadata: GEO GSE334252.

library(vroom)
library(limma)
library(minfi)
library(DMRcate)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(GenomicRanges)
library(annotatr)
library(gridExtra)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)

DATA_DIR    <- "path_to_data"
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")

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
rownames(targets_full) <- colnames(Mval)
beta <- beta[, colnames(Mval)]

Age   <- targets_full$Age
Array <- as.factor(targets_full$Array)
Slide <- as.factor(targets_full$Slide)

Fragmentation <- targets_full$Fragmentation

# Sperm parameter classes. Adapt the labels to the GEO table if needed.
param_levels <- c("OAT_OT_OA", "T_isolated", "NI", "AT_A")
param_raw <- targets_full$Sperm_param_A_astheno_T_terato_O_oligo_NI_not_identified
unexpected <- setdiff(unique(param_raw[!is.na(param_raw)]), param_levels)
if (length(unexpected) > 0) {
  stop("Unexpected sperm parameter values: ", paste(unexpected, collapse = ", "))
}
Param <- factor(param_raw, levels = param_levels)

# Samples with missing values are excluded from the corresponding models
keep_param <- !is.na(Param)
keep_frag  <- !is.na(Fragmentation)

# ===========================================================
# 2. Helper: limma fit and number of significant CpGs
# ===========================================================
# Without contrasts, tests coefficient 2 of the design. With contrasts, tests
# each contrast. `samples` selects the Mval columns matching the design rows.
# Returns the topTable of each test.
run_dmp <- function(label, design, contrasts = NULL,
                    samples = rep(TRUE, ncol(Mval))) {
  fit <- lmFit(Mval[, samples], design)
  
  if (is.null(contrasts)) {
    fit2  <- eBayes(fit)
    coefs <- 2
    names(coefs) <- colnames(design)[2]
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
# 3. Sperm parameter classes (each class vs the mean of the three others)
# ===========================================================
contr_param <- unname(sapply(param_levels, function(l) {
  others <- setdiff(param_levels, l)
  paste0("Param", l, " - (", paste0("Param", others, collapse = " + "), ")/3")
}))

run_dmp("Sperm parameter classes",
        model.matrix(~ 0 + Param + Array + Slide), contr_param,
        samples = keep_param)
run_dmp("Sperm parameter classes, adjusted for age",
        model.matrix(~ 0 + Param + Age + Array + Slide), contr_param,
        samples = keep_param)

# ===========================================================
# 4. DNA fragmentation
# ===========================================================
design_frag     <- model.matrix(~ Fragmentation + Array + Slide)
design_frag_age <- model.matrix(~ Fragmentation + Age + Array + Slide)

dmp_frag <- run_dmp("Fragmentation", design_frag, samples = keep_frag)[[1]]
run_dmp("Fragmentation, adjusted for age", design_frag_age, samples = keep_frag)

# Sensitivity analysis without the extreme individual (fragmentation = 0.82)
keep_noext <- !is.na(Fragmentation) & Fragmentation != 0.82

design_noext <- model.matrix(~ Fragmentation + Array + Slide,
                             data = targets_full[keep_noext, ])
dmp_sub <- run_dmp("Fragmentation, without extreme individual",
                   design_noext, samples = keep_noext)[[1]]

# Significant CpGs (all individuals) and their results without the extreme
sig <- dmp_frag[dmp_frag$adj.P.Val < 0.05, ]
m   <- match(rownames(sig), rownames(dmp_sub))

table_frag <- data.frame(
  CpG_full              = rownames(sig),
  CpG_ID                = str_extract(rownames(sig), "^cg[0-9]+"),
  Suffix                = str_extract(rownames(sig), "(?<=_)[A-Za-z0-9]+$"),
  Slope_all             = sig$logFC,
  Slope_without_extreme = dmp_sub$logFC[m],
  FDR_all               = sig$adj.P.Val,
  FDR_without_extreme   = dmp_sub$adj.P.Val[m],
  Direction             = ifelse(sig$logFC < 0, "Hypomethylated", "Hypermethylated"),
  stringsAsFactors = FALSE
)

# Gene annotation (EPICv2); fall back on the probe ID without suffix
annEPIC <- data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))

clean_genes <- function(x) {
  sapply(x, function(g) {
    if (is.na(g) || g == "") return(NA_character_)
    paste(unique(strsplit(g, ";")[[1]]), collapse = ";")
  }, USE.NAMES = FALSE)
}

idx <- match(table_frag$CpG_full, rownames(annEPIC))
idx[is.na(idx)] <- match(table_frag$CpG_ID[is.na(idx)],
                         sub("_.*$", "", rownames(annEPIC)))
table_frag$Gene <- clean_genes(annEPIC$UCSC_RefGene_Name[idx])

# Robust: still significant without the extreme, same direction
table_frag$Robust <- table_frag$FDR_without_extreme < 0.05 &
  sign(table_frag$Slope_all) == sign(table_frag$Slope_without_extreme)

table_frag$Slope_all             <- round(table_frag$Slope_all, 4)
table_frag$Slope_without_extreme <- round(table_frag$Slope_without_extreme, 4)
table_frag$FDR_all               <- formatC(table_frag$FDR_all, format = "e", digits = 2)
table_frag$FDR_without_extreme   <- formatC(table_frag$FDR_without_extreme, format = "e", digits = 2)

print(table_frag)
cat("CpGs still significant without the extreme:", sum(table_frag$Robust), "/", nrow(table_frag), "\n")

# Grid of the top CpGs vs fragmentation
plot_top_cpgs <- function(dmp, beta_mat, frag, n, ncol, point_color,
                          line_color, fill_color, digits) {
  plots <- lapply(rownames(dmp)[1:n], function(cpg) {
    subtitle_text <- paste0("Slope = ", round(dmp[cpg, "logFC"], 4),
                            "\nAdj.P = ", round(dmp[cpg, "adj.P.Val"], digits))
    
    df_plot <- data.frame(Beta = as.numeric(beta_mat[cpg, ]), Fragmentation = frag)
    
    ggplot(df_plot, aes(x = Fragmentation, y = Beta)) +
      geom_point(size = 2, alpha = 0.6, color = point_color) +
      geom_smooth(method = "lm", se = TRUE, color = line_color, fill = fill_color, alpha = 0.3) +
      labs(
        title = str_wrap(cpg, width = 15),
        subtitle = subtitle_text,
        x = "Fragmentation (%)",
        y = "Beta value"
      ) +
      theme_classic() +
      theme(
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")
      )
  })
  grid.arrange(grobs = plots, ncol = ncol)
}

# Top 9 CpGs without the extreme individual
plot_top_cpgs(dmp_sub, beta[, keep_noext], Fragmentation[keep_noext], n = 9, ncol = 3,
              point_color = "#2c3e50", line_color = "#e74c3c", fill_color = "#fab1a0",
              digits = 4)

# Top 10 CpGs, all individuals
plot_top_cpgs(dmp_frag, beta[, keep_frag], Fragmentation[keep_frag], n = 10, ncol = 5,
              point_color = "black", line_color = "red", fill_color = "pink",
              digits = 3)

# ===========================================================
# 5. DMPs on beta values (delta beta per fragmentation unit, in %)
# ===========================================================
make_full_dmp <- function(design, samples = rep(TRUE, ncol(Mval))) {
  dmp_m <- topTable(eBayes(lmFit(Mval[, samples], design)), number = Inf, coef = 2)
  dmp_b <- topTable(eBayes(lmFit(beta[, samples], design)), number = Inf, coef = 2)
  cat("Significant CpGs (M-values / beta values):",
      sum(dmp_m$adj.P.Val < 0.05), "/", sum(dmp_b$adj.P.Val < 0.05), "\n")
  colnames(dmp_b) <- paste0("beta_", colnames(dmp_b))
  
  full <- merge(dmp_m, dmp_b, by = 0)
  full$beta_logFC <- 100 * full$beta_logFC
  full[order(full$adj.P.Val), ]
}

cat("\nFragmentation, adjusted for age:\n")
FullDMP_age <- make_full_dmp(design_frag_age, keep_frag)
cat("Hypomethylated:", sum(FullDMP_age$adj.P.Val < 0.05 & FullDMP_age$beta_logFC < 0),
    "| Hypermethylated:", sum(FullDMP_age$adj.P.Val < 0.05 & FullDMP_age$beta_logFC > 0), "\n")

cat("\nFragmentation:\n")
FullDMP <- make_full_dmp(design_frag, keep_frag)
cat("Hypomethylated:", sum(FullDMP$adj.P.Val < 0.05 & FullDMP$beta_logFC < 0),
    "| Hypermethylated:", sum(FullDMP$adj.P.Val < 0.05 & FullDMP$beta_logFC > 0), "\n")

# ===========================================================
# 6. Differentially methylated regions (DMRs) with fragmentation
# ===========================================================
cat("\nCalling DMRs...\n")

myannotation <- cpg.annotate("array", object = as.matrix(Mval[, keep_frag]), what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type = "differential",
                             design = design_frag, coef = 2, fdr = 0.05)

dmrcoutput <- dmrcate(myannotation, lambda = 1000, C = 2, min.cpgs = 2)
DMR_ranges <- extractRanges(dmrcoutput, genome = "hg38")

data <- merge(FullDMP, annEPIC, by.y = 0, by.x = "Row.names", all.x = TRUE)

# Mean delta beta over each DMR
Annotated_ALL_GRanges <- GRanges(as.character(data$chr),
                                 IRanges(start = c(data$pos), end = c(data$pos + 1)),
                                 BetaEffect = data$beta_logFC,
                                 ind.fdr = data$adj.P.Val,
                                 is.sig = data$adj.P.Val < 0.05,
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

print(DMR_final[1:10, c("MeanBeta", "MaxBeta", "MinBeta", "nbCpGneg", "nbCpGpos",
                        "seqnames", "overlapping.genes")])

# ===========================================================
# 7. DMRs overlapping promoters (top 3, at least 3 CpGs)
# ===========================================================
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

annotations <- build_annotations(genome = "hg38", annotations = "hg38_genes_promoters")

dm_annotated <- annotate_regions(
  regions = New_DMR_ranges,
  annotations = annotations,
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

DMR_promoters_strong <- DMR_promoters %>% filter(no.cpgs >= 3)
top_dmrs <- head(DMR_promoters_strong %>% arrange(min_smoothed_fdr), 3)

cpgs_granges <- GRanges(
  seqnames = data$chr,
  ranges = IRanges(start = data$pos, end = data$pos + 1),
  probe = data$Row.names,
  beta_logFC = data$beta_logFC,
  pval = data$P.Value,
  adj_pval = data$adj.P.Val
)

FigDMR <- list()

for (i in 1:3) {
  dmr <- top_dmrs[i, ]
  
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
  
  gene_name <- ifelse(!is.null(dmr$overlapping.genes) && !is.na(dmr$overlapping.genes),
                      as.character(dmr$overlapping.genes),
                      "Intergenic")
  
  # If several genes, keep the first one
  if (grepl(",", gene_name)) {
    gene_name <- strsplit(gene_name, ", ")[[1]][1]
  }
  cat("Gene:", gene_name, "\n")
  
  p <- ggplot(plot_data, aes(x = CpG_ordered, y = beta_logFC)) +
    geom_col(aes(fill = beta_logFC), color = "black", width = 0.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, name = "Delta beta\nper fragmentation unit (%)") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    labs(title = paste0(gene_name, " - ", dmr$seqnames, ":",
                        format(dmr$start, big.mark = ","), "-",
                        format(dmr$end, big.mark = ",")),
         subtitle = paste("FDR =", format(dmr$min_smoothed_fdr, scientific = TRUE, digits = 2),
                          "| n CpGs =", length(cpgs_in_dmr)),
         x = "CpG sites (ordered by genomic position)",
         y = "Delta beta per fragmentation unit (%)") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right"
    )
  
  print(p)
  FigDMR[[i]] <- p
}

# ===========================================================
# 8. Individual CpGs of the top DMRs vs fragmentation
# ===========================================================
Scatter_DMR <- list()

for (i in 1:nrow(top_dmrs)) {
  dmr <- top_dmrs[i, ]
  
  dmr_range <- GRanges(
    seqnames = as.character(dmr$seqnames),
    ranges = IRanges(start = dmr$start, end = dmr$end)
  )
  
  overlaps <- findOverlaps(cpgs_granges, dmr_range)
  cpgs_in_dmr <- data$Row.names[queryHits(overlaps)]
  
  beta_subset <- as.data.frame(beta[cpgs_in_dmr, , drop = FALSE])
  beta_subset$CpG <- rownames(beta_subset)
  
  plot_data_long <- beta_subset %>%
    pivot_longer(cols = -CpG, names_to = "Sample", values_to = "Beta")
  
  pheno_df <- data.frame(Sample = colnames(Mval), Fragmentation = as.numeric(Fragmentation))
  plot_data_long <- left_join(plot_data_long, pheno_df, by = "Sample")
  
  gene_name <- if (!is.null(dmr$overlapping.genes) && !is.na(dmr$overlapping.genes)) {
    strsplit(as.character(dmr$overlapping.genes), ", ")[[1]][1]
  } else { "Intergenic" }
  
  p <- ggplot(plot_data_long, aes(x = Fragmentation, y = Beta)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE, fill = "mistyrose") +
    facet_wrap(~CpG, scales = "free_y") +
    labs(
      title = paste0("Individual CpG Behavior - Gene: ", gene_name),
      subtitle = paste0("DMR: ", dmr$seqnames, ":", dmr$start, "-", dmr$end,
                        " | FDR = ", format(dmr$min_smoothed_fdr, scientific = TRUE, digits = 2)),
      x = "Fragmentation",
      y = "Beta values (Methylation %)"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p)
  Scatter_DMR[[i]] <- p
}