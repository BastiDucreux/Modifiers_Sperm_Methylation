#!/usr/bin/env Rscript
# Reproduces Figure 3 (meta-analysis of age-associated DNA methylation) from the
# meta-analysis session (guide script).
#
# Input: session_20260826.rds, which must contain:
#   all_dmp_filt   per-dataset DMP results (CpG, dataset, beta_logFC, SE_beta)
#   meta_res       meta-analysis results per CpG (CpG, meta_deltabeta, meta_pval, meta_FDR)
#   data_annotated meta-analysis results with genomic positions (chr, start, end)
#   DMR_ranges     DMRs (GRanges)
#   DMR_final      DMR table (seqnames, start, end, no.cpgs, min_smoothed_fdr, MeanDeltaBeta)
# The "This dataset" cohort is deposited on GEO: GSE334252.

library(metafor)
library(GenomicRanges)
library(IRanges)
library(missMethyl)
library(annotatr)
library(cowplot)
library(dplyr)
library(ggplot2)

DATA_DIR     <- "path_to_data"
SESSION_FILE <- file.path(DATA_DIR, "session_20260826.rds")
OUT_DIR      <- file.path(DATA_DIR, "Figure3")
dir.create(OUT_DIR, showWarnings = FALSE)

session <- readRDS(SESSION_FILE)

all_dmp_filt   <- session$all_dmp_filt
meta_res       <- session$meta_res
data_annotated <- session$data_annotated
DMR_ranges     <- session$DMR_ranges
DMR_final      <- session$DMR_final
rm(session)

# Public GEO series keep their accession; the cohort of this study is "This dataset"
all_dmp_filt$dataset <- ifelse(grepl("^GSE", all_dmp_filt$dataset),
                               all_dmp_filt$dataset, "This dataset")

# ===========================================================
# 1. Volcano plot (y axis compressed above -log10(p) = 40)
# ===========================================================
volcano_data <- meta_res %>%
  mutate(
    significance = ifelse(meta_FDR < 0.05, "Significant", "Not significant"),
    direction = ifelse(meta_deltabeta < 0, "Hypomethylation", "Hypermethylation"),
    `Differential methylation` = ifelse(
      significance == "Not significant",
      "Not significant",
      ifelse(direction == "Hypomethylation", "Hypomethylation", "Hypermethylation")
    ),
    log10p = -log10(meta_pval)
  )

break_point <- 40
compression_factor <- 5

compress_scale <- function(x, bp = 40, cf = 5) {
  ifelse(x <= bp, x, bp + (x - bp) / cf)
}

volcano_data <- volcano_data %>%
  mutate(log10p_transformed = compress_scale(log10p, bp = break_point, cf = compression_factor))

n_hypo  <- sum(volcano_data$`Differential methylation` == "Hypomethylation", na.rm = TRUE)
n_hyper <- sum(volcano_data$`Differential methylation` == "Hypermethylation", na.rm = TRUE)
max_transformed_y <- max(volcano_data$log10p_transformed, na.rm = TRUE)

# Axis labels show real values at the transformed positions
y_breaks_all <- sort(unique(c(seq(0, break_point, by = 10), seq(50, 200, by = 50))))
y_breaks_transformed <- compress_scale(y_breaks_all, bp = break_point, cf = compression_factor)

FigVolcano_meta <- ggplot(volcano_data, aes(x = meta_deltabeta, y = log10p_transformed)) +
  geom_point(aes(col = `Differential methylation`), size = 1, alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Not significant" = "black",
      "Hypomethylation" = "blue",
      "Hypermethylation" = "red"
    )
  ) +
  scale_y_continuous(
    breaks = y_breaks_transformed,
    labels = y_breaks_all,
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  annotate(
    "text",
    x = min(volcano_data$meta_deltabeta, na.rm = TRUE),
    y = max_transformed_y * 1.05,
    label = paste0("Hypomethylated CpGs: ", n_hypo),
    hjust = 0, vjust = 0, size = 4.5, color = "blue", fontface = "bold"
  ) +
  annotate(
    "text",
    x = max(volcano_data$meta_deltabeta, na.rm = TRUE),
    y = max_transformed_y * 1.05,
    label = paste0("Hypermethylated CpGs: ", n_hyper),
    hjust = 1, vjust = 0, size = 4.5, color = "red", fontface = "bold"
  ) +
  labs(
    x = "Delta beta per year (%)",
    y = "-log10(p-value)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ===========================================================
# 2. Forest plots (most significant hyper- and hypomethylated CpG)
# ===========================================================
top_up <- meta_res %>%
  dplyr::filter(meta_FDR < 0.05, meta_deltabeta > 0) %>%
  dplyr::arrange(meta_pval) %>%
  head(1) %>%
  pull(CpG)

top_down <- meta_res %>%
  dplyr::filter(meta_FDR < 0.05, meta_deltabeta < 0) %>%
  dplyr::arrange(meta_pval) %>%
  head(1) %>%
  pull(CpG)

make_forest_plot <- function(cpg, status) {
  cpg_data <- all_dmp_filt[all_dmp_filt$CpG == cpg, ]
  
  # Random-effects meta-analysis (DerSimonian-Laird)
  res_plot <- rma(
    yi = beta_logFC,
    sei = SE_beta,
    data = cpg_data,
    slab = dataset,
    method = "DL"
  )
  
  col_main <- ifelse(status == "HYPER", "#CC0000", "#0066CC")
  
  forest_data <- data.frame(
    study = cpg_data$dataset,
    estimate = cpg_data$beta_logFC,
    lower = cpg_data$beta_logFC - 1.96 * cpg_data$SE_beta,
    upper = cpg_data$beta_logFC + 1.96 * cpg_data$SE_beta,
    weight = 1 / cpg_data$SE_beta^2
  )
  
  forest_data <- rbind(
    forest_data,
    data.frame(
      study = sprintf("RE Model (I\u00b2 = %.1f%%)", res_plot$I2),
      estimate = res_plot$beta[1],
      lower = res_plot$ci.lb,
      upper = res_plot$ci.ub,
      weight = NA
    )
  )
  
  # Reverse order so that the pooled estimate is at the bottom
  forest_data$study <- factor(forest_data$study, levels = rev(forest_data$study))
  forest_data$is_summary <- c(rep(FALSE, nrow(forest_data) - 1), TRUE)
  
  ggplot(forest_data, aes(x = estimate, y = study)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper, color = is_summary),
      height = 0.3,
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = is_summary, size = is_summary, shape = is_summary)
    ) +
    scale_color_manual(
      values = c("FALSE" = col_main, "TRUE" = col_main),
      guide = "none"
    ) +
    scale_size_manual(
      values = c("FALSE" = 3, "TRUE" = 4),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c("FALSE" = 15, "TRUE" = 18),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(-2.5, 2.5),
      breaks = seq(-2.5, 2.5, by = 1.25)
    ) +
    labs(
      title = cpg,
      x = "Delta beta per year (%)",
      y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12, color = col_main, hjust = 0.5),
      axis.title.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )
}

FigForest_hyper_meta <- make_forest_plot(top_up, "HYPER")
FigForest_hypo_meta  <- make_forest_plot(top_down, "HYPO")

# ===========================================================
# 3. GO enrichment of the DMRs
# ===========================================================
all.cpg <- meta_res$CpG

gst_dmr <- goregion(
  DMR_ranges,
  all.cpg = all.cpg,
  collection = "GO",
  array.type = c("EPIC"),
  plot.bias = FALSE
)
gst_dmr <- gst_dmr %>% dplyr::filter(FDR < 0.05)

FigGO_DMR_meta <- gst_dmr %>%
  group_by(ONTOLOGY) %>%
  top_n(10, -FDR) %>%
  ggplot(aes(x = FDR, y = reorder(TERM, -FDR), color = FDR, size = DE)) +
  xlim(0, 0.05) +
  facet_wrap(~ ONTOLOGY, ncol = 1, scales = "free") +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  ylab("") +
  xlab("") +
  ggtitle("GO enrichment analysis")

# ===========================================================
# 4. Top promoter DMRs (most significant hyper- and hypomethylated)
# ===========================================================

# Per-dataset delta beta of each CpG (lines in the DMR plots)
datasets_list <- all_dmp_filt %>%
  dplyr::group_by(dataset, CpG) %>%
  dplyr::summarise(beta_logFC = mean(beta_logFC, na.rm = TRUE), .groups = "drop") %>%
  dplyr::rename(Row.names = CpG) %>%
  split(.$dataset)

Annotated_ALL_GRanges <- GRanges(
  seqnames = as.character(data_annotated$chr),
  ranges = IRanges(start = data_annotated$start, end = data_annotated$end),
  DeltaBeta = data_annotated$meta_deltabeta,
  FDR = data_annotated$meta_FDR,
  probe = data_annotated$CpG
)

New_DMR_ranges <- GRanges(
  as.character(DMR_final$seqnames),
  IRanges(start = DMR_final$start, end = DMR_final$end),
  no.cpgs = DMR_final$no.cpgs,
  min_smoothed_fdr = DMR_final$min_smoothed_fdr,
  MeanDeltaBeta = DMR_final$MeanDeltaBeta
)
GenomeInfoDb::genome(New_DMR_ranges) <- "hg19"

annotations <- build_annotations(genome = "hg19", annotations = "hg19_genes_promoters")
dm_annotated <- annotate_regions(regions = New_DMR_ranges, annotations = annotations,
                                 ignore.strand = TRUE, quiet = TRUE)

# DMRs in promoters with at least 10 CpGs
df_dm_annotated <- data.frame(dm_annotated) %>%
  dplyr::group_by(seqnames, start, end, annot.type) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(annot.symbol), annot.type == "hg19_genes_promoters", no.cpgs >= 10)

top_hyper_prom <- df_dm_annotated %>%
  dplyr::filter(MeanDeltaBeta > 0) %>% dplyr::arrange(min_smoothed_fdr) %>% head(1)
top_hypo_prom <- df_dm_annotated %>%
  dplyr::filter(MeanDeltaBeta < 0) %>% dplyr::arrange(min_smoothed_fdr) %>% head(1)

cpgs_in_dmr_row <- function(dmr_row) {
  dmr_range <- GRanges(seqnames = as.character(dmr_row$seqnames),
                       ranges = IRanges(start = dmr_row$start, end = dmr_row$end))
  Annotated_ALL_GRanges$probe[queryHits(findOverlaps(Annotated_ALL_GRanges, dmr_range))]
}

# Common colour scale for both DMR plots
all_bar_values <- data_annotated$meta_deltabeta[
  data_annotated$CpG %in% c(cpgs_in_dmr_row(top_hyper_prom), cpgs_in_dmr_row(top_hypo_prom))
]
global_min <- min(all_bar_values, na.rm = TRUE)
global_max <- max(all_bar_values, na.rm = TRUE)

# Bars: meta-analysis delta beta; lines: delta beta in each individual dataset
plot_dmr_promoter <- function(dmr_row, plot_title_prefix, global_min, global_max) {
  
  cpgs_in_dmr <- cpgs_in_dmr_row(dmr_row)
  
  meta_data <- data_annotated[data_annotated$CpG %in% cpgs_in_dmr, ] %>% arrange(start)
  meta_data$CpG_ordered <- factor(meta_data$CpG, levels = unique(meta_data$CpG))
  
  plot_data_lines <- do.call(rbind, lapply(datasets_list, function(ds) {
    ds_sub <- ds[ds$Row.names %in% cpgs_in_dmr, ]
    if (nrow(ds_sub) > 0) {
      data.frame(CpG = ds_sub$Row.names, beta_logFC = ds_sub$beta_logFC, dataset = ds$dataset[1])
    } else { NULL }
  }))
  
  if (!is.null(plot_data_lines) && nrow(plot_data_lines) > 0) {
    plot_data_lines$CpG_ordered <- factor(plot_data_lines$CpG, levels = levels(meta_data$CpG_ordered))
  }
  
  p <- ggplot() +
    geom_col(data = meta_data, aes(x = CpG_ordered, y = meta_deltabeta, fill = meta_deltabeta),
             color = "black", width = 0.7, alpha = 0.8) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(global_min, global_max),
      name = "Meta Delta Beta (%)"
    )
  
  if (!is.null(plot_data_lines) && nrow(plot_data_lines) > 0) {
    p <- p +
      geom_line(data = plot_data_lines, aes(x = CpG_ordered, y = beta_logFC, group = dataset, color = dataset), linewidth = 1) +
      geom_point(data = plot_data_lines, aes(x = CpG_ordered, y = beta_logFC, color = dataset), size = 2)
  }
  
  p +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    scale_color_brewer(palette = "Set1", name = "Individual Datasets") +
    labs(
      title = paste(plot_title_prefix, ":", dmr_row$annot.symbol),
      subtitle = paste0(dmr_row$seqnames, ":", format(dmr_row$start, big.mark = ","), "-", format(dmr_row$end, big.mark = ","),
                        " | Meta FDR = ", format(dmr_row$min_smoothed_fdr, scientific = TRUE, digits = 2)),
      x = "CpG sites (genomic position)",
      y = "Delta beta / Beta logFC"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

FigDMR_hyper_meta <- plot_dmr_promoter(top_hyper_prom, "Hypermethylated", global_min, global_max)
FigDMR_hypo_meta  <- plot_dmr_promoter(top_hypo_prom, "Hypomethylated", global_min, global_max)

# ===========================================================
# 5. Figure 3 assembly
# ===========================================================
theme_global <- theme(
  text          = element_text(size = 10),
  axis.title    = element_text(size = 10, face = "bold"),
  axis.text     = element_text(size = 9),
  legend.text   = element_text(size = 9),
  legend.title  = element_text(size = 10, face = "bold"),
  plot.title    = element_text(size = 11, face = "bold"),
  plot.subtitle = element_text(size = 10),
  strip.text    = element_text(size = 9),
  plot.margin   = margin(t = 10, r = 15, b = 10, l = 15)
)

FigVolcano_meta      <- FigVolcano_meta      + theme_global
FigForest_hyper_meta <- FigForest_hyper_meta + theme_global
FigForest_hypo_meta  <- FigForest_hypo_meta  + theme_global
FigGO_DMR_meta       <- FigGO_DMR_meta       + theme_global

Forest_combined <- plot_grid(
  FigForest_hyper_meta,
  FigForest_hypo_meta,
  ncol = 1,
  nrow = 2,
  labels = c("B", "C"),
  label_size = 16
)

# One legend shared by the two DMR plots
shared_legend <- get_legend(
  FigDMR_hypo_meta + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_hyper <- FigDMR_hyper_meta + theme(legend.position = "none")
p_hypo  <- FigDMR_hypo_meta  + theme(legend.position = "none")

p_grid <- plot_grid(
  p_hyper,
  p_hypo,
  ncol = 1,
  nrow = 2,
  label_size = 16,
  align = "v"
)

DMR_combined <- plot_grid(
  p_grid,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

DMR_combined <- DMR_combined +
  theme(
    plot.margin = margin(5, 30, 5, 5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.8, "lines")
  )

row1 <- plot_grid(
  FigVolcano_meta,
  Forest_combined,
  labels = c("A", ""),
  label_size = 20,
  ncol = 2,
  rel_widths = c(1.2, 0.8)
)

row2 <- plot_grid(
  DMR_combined,
  FigGO_DMR_meta,
  labels = c("", "F"),
  label_size = 20,
  ncol = 2,
  rel_widths = c(1, 1.1)
)

combined_final_meta <- plot_grid(
  row1,
  row2,
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 1)
)

print(combined_final_meta)

saveRDS(combined_final_meta, file.path(OUT_DIR, "FigCombined_meta.rds"))

# Figure size (inches): adjust to the journal requirements
ggsave(file.path(OUT_DIR, "Figure3.pdf"), combined_final_meta, width = 16, height = 12)

cat("Figure 3 written to:", file.path(OUT_DIR, "Figure3.pdf"), "\n")