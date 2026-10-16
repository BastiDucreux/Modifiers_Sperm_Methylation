#!/usr/bin/env Rscript
# DMRs associated with pack-years (guide script): minimal code that calls the
# regions, writes the DMR table and draws the figure of the top promoter DMRs.
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
library(egg)
library(ggtext)
library(grid)
library(gridExtra)
library(dplyr)
library(ggplot2)

DATA_DIR    <- "path_to_data"
PREPROC_DIR <- file.path(DATA_DIR, "output")
META_CSV    <- file.path(DATA_DIR, "metadata.csv")
OUT_DIR     <- file.path(DATA_DIR, "smoking")
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

Array <- as.factor(targets_full$Array)
Slide <- as.factor(targets_full$Slide)

# Missing pack-years = non-smoker
Pack_year <- ifelse(is.na(targets_full$Pack_year), 0, targets_full$Pack_year)

# ===========================================================
# 2. DMPs with pack-years (M-values and beta values)
# ===========================================================
design <- model.matrix(~ Pack_year + Array + Slide)

DMPs    <- topTable(eBayes(lmFit(Mval, design)), number = Inf, coef = 2)
DMPbeta <- topTable(eBayes(lmFit(beta, design)), number = Inf, coef = 2)
colnames(DMPbeta) <- paste0("beta_", colnames(DMPbeta))

FullDMP <- merge(DMPs, DMPbeta, by = 0)
# Delta beta per pack-year, in %
FullDMP$beta_logFC <- 100 * FullDMP$beta_logFC

# CpG positions
annEPIC <- data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
data <- merge(FullDMP, annEPIC[, c("chr", "pos")], by.x = "Row.names", by.y = 0, all.x = TRUE)

# ===========================================================
# 3. DMRs
# ===========================================================
cat("Calling DMRs...\n")

myannotation <- cpg.annotate("array", object = as.matrix(Mval), what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type = "differential",
                             design = design, coef = 2, fdr = 0.05)

dmrcoutput <- dmrcate(myannotation, lambda = 1000, C = 2, min.cpgs = 2)
DMR_ranges <- extractRanges(dmrcoutput, genome = "hg38")

# ===========================================================
# 4. Delta beta over each DMR
# ===========================================================
Annotated_ALL_GRanges <- GRanges(as.character(data$chr),
                                 IRanges(start = c(data$pos), end = c(data$pos + 1)),
                                 BetaEffect = data$beta_logFC,
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
  # Direction (-1 / +1) if at least 75% of the CpGs go the same way, otherwise 0
  mutate(sign = ifelse(ratio >= 0.75 & nbCpGneg > nbCpGpos, -1,
                       ifelse(ratio >= 0.75 & nbCpGneg < nbCpGpos, 1, 0))) %>%
  mutate(consistency = (sign == sign(MeanBeta)))

DMR_final <- merge(DMR_beta_mean_per_region, DMR_ranges,
                   by.x = "subjectHits", by.y = 0)
DMR_final <- DMR_final[order(DMR_final$min_smoothed_fdr), ]

cat("Number of DMRs:", nrow(DMR_final), "\n")
print(head(DMR_final[, c("seqnames", "start", "end", "no.cpgs", "min_smoothed_fdr",
                         "MeanBeta", "overlapping.genes")], 10))

vroom_write(DMR_final, file.path(OUT_DIR, "DMR_pack_year.tbl"))

# ===========================================================
# 5. Figure: top 3 DMRs overlapping promoters (at least 3 CpGs)
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

top_dmrs <- df_dm_annotated %>%
  filter(!is.na(annot.symbol), annot.type == "hg38_genes_promoters", no.cpgs >= 3) %>%
  arrange(min_smoothed_fdr) %>%
  head(3)

cpgs_granges <- GRanges(
  seqnames = data$chr,
  ranges = IRanges(start = data$pos, end = data$pos + 1),
  probe = data$Row.names,
  beta_logFC = data$beta_logFC
)

# Physical panel width per CpG (cm); adjust to taste
bar_slot_cm <- 0.6

FigDMR <- list()

for (i in seq_len(nrow(top_dmrs))) {
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
  
  plot_data <- data[data$Row.names %in% cpgs_in_dmr, ] %>% arrange(pos)
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
                         midpoint = 0, name = "Delta beta\nper pack-year (%)") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    labs(title = paste0(gene_label, " - ", dmr$seqnames, ":",
                        format(dmr$start, big.mark = ","), "-",
                        format(dmr$end, big.mark = ",")),
         subtitle = paste("FDR =", format(dmr$min_smoothed_fdr, scientific = TRUE, digits = 2)),
         x = "CpG sites (ordered by genomic position)",
         y = "Delta beta per pack-year (%)") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      plot.title = ggtext::element_markdown(size = 12, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right"
    )
  
  # Panel width proportional to the number of CpGs (same bar width in all DMRs)
  FigDMR[[i]] <- set_panel_size(p,
                                width  = unit(n_cpgs * bar_slot_cm, "cm"),
                                height = unit(6, "cm"))
}

for (i in seq_along(FigDMR)) {
  saveRDS(FigDMR[[i]], file.path(OUT_DIR, paste0("FigDMR", i, ".rds")))
}

combined_dmrs <- arrangeGrob(grobs = FigDMR, ncol = 3)

grid.newpage()
grid.draw(combined_dmrs)

saveRDS(combined_dmrs, file.path(OUT_DIR, "FigTobacco_DMRs.rds"))

# Figure size (inches): adjust to the journal requirements
pdf(file.path(OUT_DIR, "FigTobacco_DMRs.pdf"), width = 18, height = 4)
grid.draw(combined_dmrs)
dev.off()