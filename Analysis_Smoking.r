#!/usr/bin/env Rscript
# Smoking and DNA methylation (guide script): smoking status, cigarettes per
# day, duration, pack-years and HbCO, with and without age adjustment; DMPs and
# DMRs associated with pack-years.
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
library(patchwork)
library(stringr)
library(tidyr)
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
rownames(targets_full) <- colnames(Mval)
beta <- beta[, colnames(Mval)]

Age   <- targets_full$Age
Array <- as.factor(targets_full$Array)
Slide <- as.factor(targets_full$Slide)

# Smoking status, reference first. Adapt the labels to the GEO table if needed.
smoking_levels <- c("Non_smoker", "Smoker", "Former_smoker")
unexpected <- setdiff(unique(targets_full$Smoking), smoking_levels)
if (length(unexpected) > 0) {
  stop("Unexpected Smoking values: ", paste(unexpected, collapse = ", "))
}
Smoking <- factor(targets_full$Smoking, levels = smoking_levels)

# Missing consumption = no consumption
zero_if_na <- function(x) ifelse(is.na(x), 0, x)
Cigarettes_per_day <- zero_if_na(targets_full$Daily_cigarette_consumption)
Duration           <- zero_if_na(targets_full$Duration_of_consumption)
Pack_year          <- zero_if_na(targets_full$Pack_year)
HbCO               <- targets_full$HbCO

# ===========================================================
# 2. Helper: limma fit and number of significant CpGs
# ===========================================================
# Without contrasts, tests coefficient 2 of the design. With contrasts, tests
# each contrast. Returns the topTable of each test.
run_dmp <- function(label, design, contrasts = NULL) {
  fit <- lmFit(Mval, design)
  
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
# 3. Smoking status (reference: non-smokers)
# ===========================================================
contr_smoking <- c("SmokingSmoker - SmokingNon_smoker",
                   "SmokingFormer_smoker - SmokingNon_smoker")

run_dmp("Smoking status",
        model.matrix(~ 0 + Smoking + Array + Slide), contr_smoking)
run_dmp("Smoking status, adjusted for age",
        model.matrix(~ 0 + Smoking + Age + Array + Slide), contr_smoking)

run_dmp("Current + former smokers vs non-smokers, adjusted for age",
        model.matrix(~ 0 + Smoking + Age + Array + Slide),
        "(SmokingSmoker + SmokingFormer_smoker)/2 - SmokingNon_smoker")

# ===========================================================
# 4. Cigarettes per day
# ===========================================================
run_dmp("Cigarettes per day (continuous)",
        model.matrix(~ Cigarettes_per_day + Array + Slide))
run_dmp("Cigarettes per day (continuous), adjusted for age",
        model.matrix(~ Cigarettes_per_day + Age + Array + Slide))

# >= 20 per day vs 0 (every observed value below 20 is its own group)
cigarettes_cat <- ifelse(Cigarettes_per_day >= 20, "over20", Cigarettes_per_day)
contr_cigarettes <- "cigarettes_catover20 - cigarettes_cat0"

run_dmp(">= 20 cigarettes per day vs 0",
        model.matrix(~ 0 + cigarettes_cat + Array + Slide), contr_cigarettes)
run_dmp(">= 20 cigarettes per day vs 0, adjusted for age",
        model.matrix(~ 0 + cigarettes_cat + Age + Array + Slide), contr_cigarettes)

# ===========================================================
# 5. Duration of consumption
# ===========================================================
run_dmp("Duration of consumption",
        model.matrix(~ Duration + Array + Slide))
run_dmp("Duration of consumption, adjusted for age",
        model.matrix(~ Duration + Age + Array + Slide))

# ===========================================================
# 6. Pack-years
# ===========================================================
design_pack     <- model.matrix(~ Pack_year + Array + Slide)
design_pack_age <- model.matrix(~ Pack_year + Age + Array + Slide)

dmp_pack <- run_dmp("Pack-years", design_pack)[[1]]
run_dmp("Pack-years, adjusted for age", design_pack_age)

# Sensitivity analysis without extreme individuals (>= 60 pack-years)
keep          <- Pack_year < 60
Mval_sub      <- Mval[, keep]
Pack_year_sub <- Pack_year[keep]
Array_sub     <- Array[keep]
Slide_sub     <- Slide[keep]

design_sub <- model.matrix(~ Pack_year_sub + Array_sub + Slide_sub)
dmp_sub    <- topTable(eBayes(lmFit(Mval_sub, design_sub)), number = Inf, coef = 2)
cat("\nDMPs without extreme individuals:", sum(dmp_sub$adj.P.Val < 0.05), "\n")

# Significant CpGs (all individuals) and their results without extremes
sig <- dmp_pack[dmp_pack$adj.P.Val < 0.05, ]
m   <- match(rownames(sig), rownames(dmp_sub))

table_cpgs <- data.frame(
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

idx <- match(table_cpgs$CpG_full, rownames(annEPIC))
idx[is.na(idx)] <- match(table_cpgs$CpG_ID[is.na(idx)],
                         sub("_.*$", "", rownames(annEPIC)))
table_cpgs$Gene <- clean_genes(annEPIC$UCSC_RefGene_Name[idx])

# Robust: still significant without extremes, same direction
table_cpgs$Robust <- table_cpgs$FDR_without_extreme < 0.05 &
  sign(table_cpgs$Slope_all) == sign(table_cpgs$Slope_without_extreme)

table_cpgs$Slope_all             <- round(table_cpgs$Slope_all, 4)
table_cpgs$Slope_without_extreme <- round(table_cpgs$Slope_without_extreme, 4)
table_cpgs$FDR_all               <- formatC(table_cpgs$FDR_all, format = "e", digits = 2)
table_cpgs$FDR_without_extreme   <- formatC(table_cpgs$FDR_without_extreme, format = "e", digits = 2)

print(table_cpgs)
cat("CpGs still significant without extremes:", sum(table_cpgs$Robust), "/", nrow(table_cpgs), "\n")

# ===========================================================
# 7. Top 10 CpGs vs pack-years (all individuals and without extremes)
# ===========================================================
top10_cpgs <- rownames(dmp_pack)[1:10]

plot_list <- list()
for (cpg in top10_cpgs) {
  
  # Slope and FDR: all individuals, then (without extremes) in brackets
  logFC_val <- dmp_pack[cpg, "logFC"]
  adj_p     <- dmp_pack[cpg, "adj.P.Val"]
  
  logFC_val_sub <- dmp_sub[cpg, "logFC"]
  adj_p_sub     <- dmp_sub[cpg, "adj.P.Val"]
  
  # Blue = hypomethylation, red = hypermethylation
  line_color <- ifelse(logFC_val < 0, "blue", "red")
  fill_color <- ifelse(logFC_val < 0, "lightblue", "pink")
  
  df_plot <- data.frame(
    Beta = as.numeric(beta[cpg, ]),
    PackYears = Pack_year,
    Extreme = Pack_year >= 60
  )
  
  subtitle_text <- paste0("Slope = ", round(logFC_val, 4), " (", round(logFC_val_sub, 4), ")",
                          "\nFDR = ", formatC(adj_p, format = "e", digits = 2),
                          " (", formatC(adj_p_sub, format = "e", digits = 2), ")")
  
  p <- ggplot(df_plot, aes(x = PackYears, y = Beta)) +
    geom_point(aes(color = Extreme), size = 2, alpha = 0.6) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "grey60"), guide = "none") +
    geom_smooth(method = "lm", se = TRUE,
                color = line_color,
                fill = fill_color,
                alpha = 0.3) +
    geom_smooth(data = subset(df_plot, !Extreme), aes(x = PackYears, y = Beta),
                method = "lm", se = TRUE,
                color = "darkgreen",
                fill = "lightgreen",
                alpha = 0.2, inherit.aes = FALSE) +
    labs(
      title = str_wrap(cpg, width = 15),
      subtitle = subtitle_text,
      x = "Pack-years",
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
  
  plot_list[[cpg]] <- p
}

FigTobacco_DMP <- wrap_plots(plot_list, ncol = 5)
saveRDS(FigTobacco_DMP, file.path(OUT_DIR, "FigTobacco_DMP.rds"))

# ===========================================================
# 8. DMPs on beta values (delta beta per pack-year, in %)
# ===========================================================
make_full_dmp <- function(design) {
  dmp_m <- topTable(eBayes(lmFit(Mval, design)), number = Inf, coef = 2)
  dmp_b <- topTable(eBayes(lmFit(beta, design)), number = Inf, coef = 2)
  cat("Significant CpGs (M-values / beta values):",
      sum(dmp_m$adj.P.Val < 0.05), "/", sum(dmp_b$adj.P.Val < 0.05), "\n")
  colnames(dmp_b) <- paste0("beta_", colnames(dmp_b))
  
  full <- merge(dmp_m, dmp_b, by = 0)
  full$beta_logFC <- 100 * full$beta_logFC
  full[order(full$adj.P.Val), ]
}

cat("\nPack-years:\n")
FullDMP <- make_full_dmp(design_pack)
cat("Hypomethylated:", sum(FullDMP$adj.P.Val < 0.05 & FullDMP$beta_logFC < 0),
    "| Hypermethylated:", sum(FullDMP$adj.P.Val < 0.05 & FullDMP$beta_logFC > 0), "\n")

cat("\nPack-years, adjusted for age:\n")
FullDMP_age <- make_full_dmp(design_pack_age)
cat("Hypomethylated:", sum(FullDMP_age$adj.P.Val < 0.05 & FullDMP_age$beta_logFC < 0),
    "| Hypermethylated:", sum(FullDMP_age$adj.P.Val < 0.05 & FullDMP_age$beta_logFC > 0), "\n")

# ===========================================================
# 9. Differentially methylated regions (DMRs) with pack-years
# ===========================================================
cat("\nCalling DMRs...\n")

myannotation <- cpg.annotate("array", object = as.matrix(Mval), what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type = "differential",
                             design = design_pack, coef = 2, fdr = 0.05)

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
# 10. DMRs overlapping promoters (top 3, at least 3 CpGs)
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

# Physical panel width per CpG (cm); adjust to taste
bar_slot_cm <- 0.6

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

for (i in 1:3) {
  saveRDS(FigDMR[[i]], file.path(OUT_DIR, paste0("FigDMR", i, ".rds")))
}

combined_dmrs <- arrangeGrob(FigDMR[[1]], FigDMR[[2]], FigDMR[[3]], ncol = 3)
grid.newpage()
grid.draw(combined_dmrs)
saveRDS(combined_dmrs, file.path(OUT_DIR, "FigTobacco_DMRs.rds"))

# ===========================================================
# 11. Individual CpGs of the top DMRs vs pack-years
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
  
  pheno_df <- data.frame(Sample = colnames(Mval), PackYears = Pack_year)
  plot_data_long <- left_join(plot_data_long, pheno_df, by = "Sample")
  
  gene_name <- if (!is.null(dmr$overlapping.genes) && !is.na(dmr$overlapping.genes)) {
    strsplit(as.character(dmr$overlapping.genes), ", ")[[1]][1]
  } else { "Intergenic" }
  
  p <- ggplot(plot_data_long, aes(x = PackYears, y = Beta)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE, fill = "mistyrose") +
    facet_wrap(~CpG, scales = "free_y") +
    labs(
      title = paste0("Individual CpG Behavior - Gene: ", gene_name),
      subtitle = paste0("DMR: ", dmr$seqnames, ":", dmr$start, "-", dmr$end,
                        " | FDR = ", format(dmr$min_smoothed_fdr, scientific = TRUE, digits = 2)),
      x = "Pack-years",
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

# ===========================================================
# 12. Annotated table of pack-year DMPs
# ===========================================================
Full_DMP_Annotated <- merge(
  dmp_pack,
  annEPIC[, c("UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island", "chr", "pos")],
  by = 0
)
colnames(Full_DMP_Annotated)[1] <- "CpG"

sig_dmps <- Full_DMP_Annotated[Full_DMP_Annotated$adj.P.Val < 0.05, ]
sig_dmps <- sig_dmps[order(sig_dmps$adj.P.Val), ]
sig_dmps$UCSC_RefGene_Name <- sapply(strsplit(as.character(sig_dmps$UCSC_RefGene_Name), ";"), `[`, 1)

cat("\nNumber of DMPs:", nrow(sig_dmps), "\n")
print(head(sig_dmps))

# ===========================================================
# 13. HbCO
# ===========================================================
run_dmp("HbCO (continuous)",
        model.matrix(~ HbCO + Array + Slide))
run_dmp("HbCO (continuous), adjusted for age",
        model.matrix(~ HbCO + Age + Array + Slide))

# HbCO < 2 vs >= 2
HbCO_bin <- factor(ifelse(HbCO < 2, "<2", ">=2"))

run_dmp("HbCO >= 2 vs < 2",
        model.matrix(~ HbCO_bin + Array + Slide))
run_dmp("HbCO >= 2 vs < 2, adjusted for age",
        model.matrix(~ HbCO_bin + Age + Array + Slide))