#!/usr/bin/env Rscript
# Reproduces Supplementary Figure 2 (guide script): global methylation of the
# repetitive elements LINE-1 and Alu vs age, BMI, smoking, pack-years, HbCO,
# residential environment, professional exposure, DNA fragmentation, sperm
# parameters (OATS) and early pregnancy.
#
# Statistical tests are run on M-values; plots show beta values.
#
# Inputs (raw data and metadata: GEO GSE334252):
#   - metadata.csv (GEO sample table)
#   - RE.tbl: repetitive element methylation per sample (columns Basename,
#     beta_Line1.Mean, Mval_Line1.Mean, beta_Alu.Mean, Mval_Alu.Mean)

set.seed(100)

library(ggplot2)
library(ggpubr)
library(patchwork)

DATA_DIR <- "path_to_data"
META_CSV <- file.path(DATA_DIR, "metadata.csv")
RE_FILE  <- file.path(DATA_DIR, "RE.tbl")
OUT_DIR  <- file.path(DATA_DIR, "FigureS2")
dir.create(OUT_DIR, showWarnings = FALSE)

# ===========================================================
# 1. Load data
# ===========================================================
RE <- read.table(RE_FILE, header = TRUE)

targets_full <- read.csv(META_CSV, header = TRUE, na.strings = c("", "NA"))

# Sample names in RE.tbl start with "X"; <Slide>_<Array> comes from the idat name
targets_full$ID <- paste0("X", sub("_Grn\\.idat$", "", targets_full$idat_file_Grn))

NEW_RE <- merge(targets_full, RE, by.x = "ID", by.y = "Basename")

# ===========================================================
# 2. Variables
# ===========================================================
pd <- data.frame(
  Age = as.numeric(NEW_RE$Age),
  BMI = as.numeric(NEW_RE$BMI),
  # Missing pack-years = non-smoker
  Pack_year = as.numeric(ifelse(is.na(NEW_RE$Pack_year), 0, NEW_RE$Pack_year)),
  Smoking = factor(NEW_RE$Smoking,
                   levels = c("Smoker", "Former_smoker", "Non_smoker"),
                   labels = c("Active\nsmoking", "Former\nsmoking", "No\nsmoking")),
  HbCO_concentration = as.numeric(NEW_RE$HbCO),
  Fragmentation = as.numeric(NEW_RE$Fragmentation),
  OATS = factor(ifelse(is.na(NEW_RE$Sperm_param_A_astheno_T_terato_O_oligo_NI_not_identified),
                       NA,
                       NEW_RE$Sperm_param_A_astheno_T_terato_O_oligo_NI_not_identified == "OAT_OT_OA"),
                levels = c(TRUE, FALSE),
                labels = c("Yes", "No")),
  Pregnancy_outcome = factor(NEW_RE$Early_pregnancy, levels = c("Yes", "No")),
  Professional_exposure = factor(NEW_RE$ExpositionActual_or_Previous_professional,
                                 levels = c("Yes", "No")),
  Residential_environment = factor(NEW_RE$Residential_environment,
                                   levels = c("Urban", "Rural")),
  
  beta_Line1_Mean = NEW_RE$beta_Line1.Mean,
  Mval_Line1_Mean = NEW_RE$Mval_Line1.Mean,
  beta_Alu_Mean = NEW_RE$beta_Alu.Mean,
  Mval_Alu_Mean = NEW_RE$Mval_Alu.Mean,
  row.names = NEW_RE$ID
)

# ===========================================================
# 3. Plot settings
# ===========================================================
repeats <- list(
  LINE1 = list(name = "LINE-1", beta = "beta_Line1_Mean", mval = "Mval_Line1_Mean", color = "#7A0177"),
  Alu   = list(name = "Alu",    beta = "beta_Alu_Mean",   mval = "Mval_Alu_Mean",   color = "#CCB044")
)

# Common y limits for all panels
y_min_val <- min(c(pd$beta_Line1_Mean, pd$beta_Alu_Mean), na.rm = TRUE)
y_max_val <- max(c(pd$beta_Line1_Mean, pd$beta_Alu_Mean), na.rm = TRUE)
y_margin  <- (y_max_val - y_min_val) * 0.05
y_limits  <- c(y_min_val - y_margin, y_max_val + y_margin)

# Common theme (aligns titles and tags)
theme_align <- theme_classic() +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text    = element_text(size = 9),
    plot.title   = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 6)),
    plot.tag     = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0, 1.05)
  )

fmt_p <- function(p) if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))

# P-value of the association between a variable and a repeat (M-values)
#   cor     : Pearson correlation (continuous variables)
#   kruskal : Kruskal-Wallis (3 groups)
#   wilcox  : Wilcoxon (2 groups)
test_p <- function(var, test, mval) {
  if (test == "cor") {
    cor.test(pd[[var]], pd[[mval]])$p.value
  } else {
    f <- as.formula(paste(mval, "~", var))
    if (test == "kruskal") {
      kruskal.test(f, data = pd)$p.value
    } else {
      wilcox.test(f, data = pd, exact = FALSE)$p.value
    }
  }
}

# One panel: scatter plot (continuous variable) or boxplot (categorical)
make_panel <- function(spec, rep, tag = NULL) {
  p_lab <- fmt_p(test_p(spec$var, spec$test, rep$mval))
  title <- paste0(rep$name, " (", p_lab, ")")
  
  if (spec$test == "cor") {
    p <- ggplot(pd, aes(x = .data[[spec$var]], y = .data[[rep$beta]])) +
      geom_point(size = 2, alpha = 0.7, color = rep$color) +
      geom_smooth(method = "lm", se = TRUE, color = rep$color, fill = rep$color, alpha = 0.2) +
      scale_y_continuous(limits = y_limits) +
      labs(title = title, x = spec$xlab, y = "Beta value", tag = tag)
  } else {
    p <- ggboxplot(pd, x = spec$var, y = rep$beta,
                   ylab = "Beta value", xlab = spec$xlab,
                   add = "jitter", color = rep$color) +
      scale_y_continuous(limits = y_limits) +
      labs(title = title, tag = tag)
  }
  
  p + theme_align
}

# Covariates in the figure order (one row each); panel letters A-E then F-J
left_specs <- list(
  list(var = "Age",                 test = "cor",     xlab = "Age (years)"),
  list(var = "BMI",                 test = "cor",     xlab = "BMI (kg/m\u00b2)"),
  list(var = "Smoking",             test = "kruskal", xlab = "Smoking status"),
  list(var = "Pack_year",           test = "cor",     xlab = "Pack-years"),
  list(var = "HbCO_concentration",  test = "cor",     xlab = "HbCO concentration (%)")
)

right_specs <- list(
  list(var = "Residential_environment", test = "wilcox", xlab = "Residential environment"),
  list(var = "Professional_exposure",   test = "wilcox", xlab = "Professional exposure"),
  list(var = "Fragmentation",           test = "cor",    xlab = "DNA fragmentation (%)"),
  list(var = "OATS",                    test = "wilcox", xlab = "OAT"),
  list(var = "Pregnancy_outcome",       test = "wilcox", xlab = "Pregnancy outcome")
)

# ===========================================================
# 4. Panels and assembly
# ===========================================================
cat("Generating plots for LINE-1 and Alu...\n")

# Each row: left covariate (LINE-1, Alu) and right covariate (LINE-1, Alu);
# the panel letter is shown on the LINE-1 panel only
rows <- lapply(1:5, function(i) {
  l <- left_specs[[i]]
  r <- right_specs[[i]]
  
  make_panel(l, repeats$LINE1, tag = LETTERS[i]) |
    make_panel(l, repeats$Alu) |
    make_panel(r, repeats$LINE1, tag = LETTERS[i + 5]) |
    make_panel(r, repeats$Alu)
})

combined_plot <- Reduce(`/`, rows) &
  theme(plot.margin = margin(t = 12, r = 5, b = 5, l = 5))

print(combined_plot)

saveRDS(combined_plot, file.path(OUT_DIR, "SupplementaryFigure2.rds"))

# Figure size (inches): adjust to the journal requirements
ggsave(file.path(OUT_DIR, "SupplementaryFigure2.pdf"), combined_plot,
       width = 16, height = 18)

cat("Analysis complete!\n")