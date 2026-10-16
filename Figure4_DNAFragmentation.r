#!/usr/bin/env Rscript
# Minimal code to reproduce the DMR figure for DNA fragmentation: the three most
# significant promoter DMRs (at least 3 CpGs), shown as delta beta per CpG
# (guide script).
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
OUT_DIR     <- file.path(DATA_DIR, "fragmentation")
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

Fragmentation <- targets_full$Fragmentation
Array         <- as.factor(targets_full$Array)
Slide         <- as.factor(targets_full$Slide)
keep          <- !is.na(Fragmentation)

# ===========================================================
# 2. DMPs with DNA fragmentation (M-values and beta values)
# ===========================================================
design <- model.matrix(~ Fragmentation + Array + Slide)

DMPs    <- topTable(eBayes(lmFit(Mval[, keep], design)), number = Inf, coef = 2)
DMPbeta <- topTable(eBayes(lmFit(beta[, keep], design)), number = Inf, coef = 2)
colnames(DMPbeta) <- paste0("beta_", colnames(DMPbeta))

FullDMP <- merge(DMPs, DMPbeta, by = 0)
# Delta beta per unit of fragmentation, in %
FullDMP$beta_logFC <- 100 * FullDMP$beta_logFC

# CpG positions
annEPIC <- data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
data <- merge(FullDMP, annEPIC[, c("chr", "pos")], by.x = "Row.names", by.y = 0, all.x = TRUE)

# ===========================================================
# 3. DMRs
# ===========================================================
myannotation <- cpg.annotate("array", object = as.matrix(Mval[, keep]), what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type = "differential",
                             design = design, coef = 2, fdr = 0.05)

dmrcoutput <- dmrcate(myannotation, lambda = 1000, C = 2, min.cpgs = 2)
DMR_ranges <- extractRanges(dmrcoutput, genome = "hg38")

# ===========================================================
# 4. Top 3 promoter DMRs (at least 3 CpGs)
# ===========================================================
annotations <- build_annotations(genome = "hg38", annotations = "hg38_genes_promoters")

dm_annotated <- annotate_regions(regions = DMR_ranges, annotations = annotations,
                                 ignore.strand = TRUE, quiet = FALSE)

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

# ===========================================================
# 5. Figure
# ===========================================================
# Physical panel width per CpG (cm); adjust to taste
bar_slot_cm <- 0.6

list_dmr_plots <- list()

for (i in 1:nrow(top_dmrs)) {
  dmr <- top_dmrs[i, ]
  
  dmr_range <- GRanges(seqnames = as.character(dmr$seqnames),
                       ranges = IRanges(start = dmr$start, end = dmr$end))
  
  overlaps <- findOverlaps(cpgs_granges, dmr_range)
  cpgs_in_dmr <- data$Row.names[queryHits(overlaps)]
  
  plot_data <- data[data$Row.names %in% cpgs_in_dmr, ] %>% arrange(pos)
  plot_data$CpG_ordered <- factor(plot_data$Row.names, levels = plot_data$Row.names)
  n_cpgs <- nrow(plot_data)
  
  gene_name <- ifelse(!is.null(dmr$overlapping.genes) && !is.na(dmr$overlapping.genes),
                      as.character(dmr$overlapping.genes), "Intergenic")
  # If several genes, keep the first one
  if (grepl(",", gene_name)) {
    gene_name <- strsplit(gene_name, ", ")[[1]][1]
  }
  # Gene name in italics (markdown), except "Intergenic"
  gene_label <- ifelse(gene_name == "Intergenic", gene_name, paste0("*", gene_name, "*"))
  
  p <- ggplot(plot_data, aes(x = CpG_ordered, y = beta_logFC)) +
    geom_col(aes(fill = beta_logFC), color = "black", width = 0.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Delta beta\nper fragmentation unit (%)") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    labs(title = paste0(gene_label, " - ", dmr$seqnames, ":",
                        format(dmr$start, big.mark = ","), "-",
                        format(dmr$end, big.mark = ",")),
         subtitle = paste("FDR =", format(dmr$min_smoothed_fdr, scientific = TRUE, digits = 2)),
         x = "CpG sites", y = "Delta beta (%)") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = ggtext::element_markdown(size = 11, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right"
    )
  
  # Panel width proportional to the number of CpGs (same bar width in all DMRs)
  list_dmr_plots[[i]] <- set_panel_size(p, width = unit(n_cpgs * bar_slot_cm, "cm"),
                                        height = unit(6, "cm"))
}

combined_dmrs <- arrangeGrob(grobs = list_dmr_plots, ncol = 3)

grid.newpage()
grid.draw(combined_dmrs)

