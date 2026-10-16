#!/usr/bin/env Rscript
# Reproduces Supplementary Figure 3 (guide script): pyrosequencing methylation
# of H19, SNRPN and LINE-1 vs age (OATS excluded), BMI, smoking status and OATS.
#
# Statistical tests are run on M-values; the displayed effect is computed on
# methylation percentages (beta for regressions, delta beta for group
# comparisons).
#
# Input: Pyrosequencing_data.xlsx, written by make_pyrosequencing_data.R.

library(openxlsx)
library(ggplot2)
library(ggpubr)
library(patchwork)

DATA_DIR <- "path_to_data"
OUT_DIR  <- file.path(DATA_DIR, "FigureS3")
dir.create(OUT_DIR, showWarnings = FALSE)

# ===========================================================
# 1. Load data
# ===========================================================
pyro <- read.xlsx(file.path(DATA_DIR, "Pyrosequencing_data.xlsx"), sheet = "Data")

pyro$OATS          <- factor(pyro$OATS, levels = c("No", "Yes"))
pyro$Smoking_group <- factor(pyro$Smoking_group, levels = c("Never/Former", "Current"))

# ===========================================================
# 2. Markers, colours, labels
# ===========================================================
meth_vars <- c("H19_CTCF3B_seq1_mean", "SNRPN_mean", "LINE1_mean")

meth_labels <- list(
  expression(bolditalic("H19")),
  expression(bolditalic("SNRPN")),
  expression(bold("LINE-1"))
)

gene_colors <- c(
  "H19_CTCF3B_seq1_mean" = "#009E73",
  "SNRPN_mean"           = "#CC79A7",
  "LINE1_mean"           = "#0072B2"
)

# Same y-axis range for all panels of a gene
gene_ylims <- list(
  "H19_CTCF3B_seq1_mean" = c(95, 100),
  "SNRPN_mean"           = c(0, 10),
  "LINE1_mean"           = c(45, 75)
)

# ===========================================================
# 3. M-values (columns "_M"), used for all tests
#    beta = % / 100 ; M = log2(beta / (1 - beta))
#    eps avoids infinite values if a point is at 0 or 100 %
# ===========================================================
to_M <- function(pct, eps = 1e-4) {
  b <- pct / 100
  b <- pmin(pmax(b, eps), 1 - eps)
  log2(b / (1 - b))
}

for (v in meth_vars) {
  pyro[[paste0(v, "_M")]] <- to_M(as.numeric(pyro[[v]]))
}

# Delta beta = mean(beta, group 2) - mean(beta, group 1), beta = % / 100
delta_beta <- function(data, var, group_col) {
  beta_vals <- as.numeric(data[[var]]) / 100
  diff(tapply(beta_vals, data[[group_col]], mean, na.rm = TRUE))
}

theme_panel <- function() {
  theme_pubr() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 8),
      plot.tag   = element_text(face = "bold", size = 14)
    )
}

# ===========================================================
# 4. Panel functions
# ===========================================================
# Continuous variable: linear regression. p-value on M-values; displayed beta
# is the slope on methylation percentages (change per unit of the variable).
scatter_panel <- function(data, x_var, x_lab, i) {
  var <- meth_vars[i]
  col_gene <- gene_colors[[var]]
  
  mod    <- lm(as.formula(paste0(var, "_M ~ ", x_var)), data = data)
  pval   <- summary(mod)$coefficients[2, 4]
  effect <- summary(lm(as.formula(paste0(var, " ~ ", x_var)), data = data))$coefficients[2, 1]
  
  ggplot(data, aes(x = .data[[x_var]], y = .data[[var]])) +
    geom_point(alpha = 0.6, color = col_gene, size = 2) +
    geom_smooth(method = "lm", color = col_gene, fill = col_gene, alpha = 0.15, linewidth = 1) +
    theme_panel() +
    coord_cartesian(ylim = gene_ylims[[var]]) +
    labs(
      title    = meth_labels[[i]],
      subtitle = paste0("\u03b2 = ", round(effect, 3), ", p = ", round(pval, 4)),
      x = x_lab,
      y = "DNA methylation (%)"
    )
}

# Two groups: Wilcoxon test on M-values, delta beta on methylation percentages
box_panel <- function(data, group_var, x_lab, i) {
  var <- meth_vars[i]
  col_gene <- gene_colors[[var]]
  
  wt   <- wilcox.test(as.formula(paste0(var, "_M ~ ", group_var)),
                      data = data, exact = FALSE)
  pval <- wt$p.value
  db   <- delta_beta(data, var, group_var)
  
  ggplot(data, aes(x = .data[[group_var]], y = .data[[var]])) +
    geom_boxplot(color = col_gene, fill = col_gene, alpha = 0.2, outlier.shape = NA) +
    geom_jitter(color = col_gene, width = 0.15, alpha = 0.5, size = 1.5) +
    theme_panel() +
    coord_cartesian(ylim = gene_ylims[[var]]) +
    labs(
      title    = meth_labels[[i]],
      subtitle = paste0("\u0394\u03b2 = ", round(db, 4), ", p = ", round(pval, 4)),
      x = x_lab,
      y = "DNA methylation (%)"
    )
}

# ===========================================================
# 5. Panels
# ===========================================================
# A. Age (OATS excluded)
df_age   <- pyro[pyro$OATS != "Yes" | is.na(pyro$OATS), ]
plots_age <- lapply(seq_along(meth_vars), function(i) scatter_panel(df_age, "Age", "Age (years)", i))

# B. BMI
df_bmi   <- pyro[!is.na(pyro$BMI), ]
plots_bmi <- lapply(seq_along(meth_vars), function(i) scatter_panel(df_bmi, "BMI", "BMI", i))

# C. Smoking (current vs never/former)
df_smoking    <- pyro[!is.na(pyro$Smoking_group), ]
plots_smoking <- lapply(seq_along(meth_vars), function(i) box_panel(df_smoking, "Smoking_group", "Smoking Status", i))

# D. OATS (yes vs no)
plots_oats <- lapply(seq_along(meth_vars), function(i) box_panel(pyro, "OATS", "OAT", i))

# ===========================================================
# 6. Assembly
# ===========================================================
plots_age[[1]]     <- plots_age[[1]]     + labs(tag = "A")
plots_bmi[[1]]     <- plots_bmi[[1]]     + labs(tag = "B")
plots_smoking[[1]] <- plots_smoking[[1]] + labs(tag = "C")
plots_oats[[1]]    <- plots_oats[[1]]    + labs(tag = "D")

panel_A <- wrap_plots(plots_age,     nrow = 1)
panel_B <- wrap_plots(plots_bmi,     nrow = 1)
panel_C <- wrap_plots(plots_smoking, nrow = 1)
panel_D <- wrap_plots(plots_oats,    nrow = 1)

final_figure <- panel_A / panel_B / panel_C / panel_D +
  plot_layout(heights = c(1, 1, 1, 1))

print(final_figure)

saveRDS(final_figure, file.path(OUT_DIR, "SupplementaryFigure3.rds"))

# Figure size (inches): adjust to the journal requirements
ggsave(file.path(OUT_DIR, "SupplementaryFigure3.pdf"), final_figure,
       width = 10, height = 12)