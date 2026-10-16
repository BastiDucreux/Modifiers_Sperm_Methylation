#!/usr/bin/env Rscript
# Epigenetic age acceleration plots: generation of panels A to J
# showing the relationship between AgeAccelerationResidual and various covariates.
#
# Input: EpigeneticAge_results.rds

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)

DATA_DIR    <- "path_to_data"
RDS_FILE    <- file.path(DATA_DIR, "EpigeneticAge_results.rds")
OUT_DIR     <- file.path(DATA_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE)

cat("Loading RDS results", sep = "\n")

# ===========================================================
# 1. Load data and metadata
# ===========================================================
results <- readRDS(RDS_FILE)
pd <- results$pd

# Correction of Pregnancy_outcome: since there are only 2 possible categories,
# test only the pattern without accent ("Pas de grossesse") and everything else
# becomes "Yes". This avoids Unicode encoding issues (NFC/NFD) on the accent.
pd$Pregnancy_outcome <- factor(
  ifelse(grepl("^Pas de grossesse", results$Table_cor$Issue_tentative), "No", "Yes"),
  levels = c("Yes", "No")
)

# ── Global theme ──────────────────────────────────────────
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
  plot.title    = element_blank(),
  plot.subtitle = element_text(size = 10),
  strip.text    = element_text(size = 9),
  plot.margin   = margin(t = 8, r = 15, b = 5, l = 15)
)

# Theme for panel tags (A, B, C...)
theme_tag <- theme(plot.tag = element_text(size = 11, face = "bold"))

# Y-axis expansion to leave room for p-values
y_expand <- expansion(mult = c(0.05, 0.18))

cat("Generating plots for epigenetic age acceleration", sep = "\n")

# ===========================================================
# 2. Individual plots generation (A to J)
# ===========================================================

# ========== A: Age ==========
cor_test_age <- cor.test(pd$Age, pd$AgeAccelerationResidual)
p_value_age <- cor_test_age$p.value
p_label_age <- if (p_value_age < 0.001) "p < 0.001" else paste0("p = ", format.pval(p_value_age, digits = 3))
y_max_age <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pA <- ggplot(pd, aes(x = Age, y = AgeAccelerationResidual)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = max(pd$Age), y = y_max_age * 1.15,
           label = p_label_age, hjust = 1, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  labs(x = "Age (years)", y = "Age Acc", tag = "A") +
  theme_classic() +
  theme_global +
  theme_tag

# ========== B: Pack-years ==========
cor_test_pack <- cor.test(pd$Pack_year, pd$AgeAccelerationResidual)
p_value_pack <- cor_test_pack$p.value
p_label_pack <- if (p_value_pack < 0.001) "p < 0.001" else paste0("p = ", format.pval(p_value_pack, digits = 3))
y_max_pack <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pB <- ggplot(pd, aes(x = Pack_year, y = AgeAccelerationResidual)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = max(pd$Pack_year), y = y_max_pack * 1.15,
           label = p_label_pack, hjust = 1, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  labs(x = "Pack-years", y = "Age Acc", tag = "B") +
  theme_classic() +
  theme_global +
  theme_tag

# ========== C: Smoking ==========
res_anova_smoking <- aov(AgeAccelerationResidual ~ Smoking, data = pd)
p_val_smoking <- summary(res_anova_smoking)[[1]][["Pr(>F)"]][1]
p_label_smoking <- ifelse(p_val_smoking < 0.001, "p < 0.001", paste0("p = ", round(p_val_smoking, 3)))
y_max_smoking <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pC <- ggboxplot(pd, x = "Smoking", y = "AgeAccelerationResidual",
                ylab = "Age Acc", xlab = "Smoking status",
                add = "jitter", palette = "jco") +
  annotate("text", x = Inf, y = y_max_smoking * 1.15,
           label = p_label_smoking, color = "black", hjust = 1.05, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  coord_cartesian(clip = "off") +
  labs(tag = "C") +
  theme_pubr() +
  theme_global +
  theme_tag

# ========== D: HbCO concentration ==========
cor_test_hbco <- cor.test(pd$HbCO_concentration, pd$AgeAccelerationResidual)
p_value_hbco <- cor_test_hbco$p.value
p_label_hbco <- if (p_value_hbco < 0.001) "p < 0.001" else paste0("p = ", format.pval(p_value_hbco, digits = 3))
y_max_hbco <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pD <- ggplot(pd, aes(x = HbCO_concentration, y = AgeAccelerationResidual)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = max(pd$HbCO_concentration, na.rm = TRUE), y = y_max_hbco * 1.15,
           label = p_label_hbco, hjust = 1, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  labs(x = "HbCO concentration (%)", y = "Age Acc", tag = "D") +
  theme_classic() +
  theme_global +
  theme_tag

# ========== E: BMI ==========
cor_test_bmi <- cor.test(pd$BMI, pd$AgeAccelerationResidual)
p_value_bmi <- cor_test_bmi$p.value
p_label_bmi <- if (p_value_bmi < 0.001) "p < 0.001" else paste0("p = ", format.pval(p_value_bmi, digits = 3))
y_max_bmi <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pE <- ggplot(pd, aes(x = BMI, y = AgeAccelerationResidual)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = max(pd$BMI, na.rm = TRUE), y = y_max_bmi * 1.15,
           label = p_label_bmi, hjust = 1, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  labs(x = "BMI (kg/m²)", y = "Age Acc", tag = "E") +
  theme_classic() +
  theme_global +
  theme_tag

# ========== F: Fragmentation ==========
cor_test_frag <- cor.test(pd$Fragmentation, pd$AgeAccelerationResidual)
p_value_frag <- cor_test_frag$p.value
p_label_frag <- if (p_value_frag < 0.001) "p < 0.001" else paste0("p = ", format.pval(p_value_frag, digits = 3))
y_max_frag <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pF <- ggplot(pd, aes(x = Fragmentation, y = AgeAccelerationResidual)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = max(pd$Fragmentation, na.rm = TRUE), y = y_max_frag * 1.15,
           label = p_label_frag, hjust = 1, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  labs(x = "DNA fragmentation (%)", y = "Age Acc", tag = "F") +
  theme_classic() +
  theme_global +
  theme_tag

# ========== G: OATS ==========
res_ttest_oats <- t.test(AgeAccelerationResidual ~ OATS, data = pd)
p_val_oats <- res_ttest_oats$p.value
p_label_oats <- ifelse(p_val_oats < 0.001, "p < 0.001", paste0("p = ", round(p_val_oats, 3)))
y_max_oats <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pG <- ggboxplot(pd, x = "OATS", y = "AgeAccelerationResidual",
                ylab = "Age Acc", xlab = "OAT",
                add = "jitter", palette = "jco") +
  annotate("text", x = Inf, y = y_max_oats * 1.15,
           label = p_label_oats, color = "black", hjust = 1.05, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  coord_cartesian(clip = "off") +
  labs(tag = "G") +
  theme_pubr() +
  theme_global +
  theme_tag

# ========== H: Pregnancy outcome ==========
res_ttest_preg <- t.test(AgeAccelerationResidual ~ Pregnancy_outcome, data = pd)
p_val_preg <- res_ttest_preg$p.value
p_label_preg <- ifelse(p_val_preg < 0.001, "p < 0.001", paste0("p = ", round(p_val_preg, 3)))
y_max_preg <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pH <- ggboxplot(pd, x = "Pregnancy_outcome", y = "AgeAccelerationResidual",
                ylab = "Age Acc", xlab = "Pregnancy outcome",
                add = "jitter", palette = "jco") +
  annotate("text", x = Inf, y = y_max_preg * 1.15,
           label = p_label_preg, color = "black", hjust = 1.05, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  coord_cartesian(clip = "off") +
  labs(tag = "H") +
  theme_pubr() +
  theme_global +
  theme_tag

# ========== I: Professional exposure ==========
res_ttest_prof <- t.test(AgeAccelerationResidual ~ Professional_exposure, data = pd)
p_val_prof <- res_ttest_prof$p.value
p_label_prof <- ifelse(p_val_prof < 0.001, "p < 0.001", paste0("p = ", round(p_val_prof, 3)))
y_max_prof <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pI <- ggboxplot(pd, x = "Professional_exposure", y = "AgeAccelerationResidual",
                ylab = "Age Acc", xlab = "Professional exposure",
                add = "jitter", palette = "jco") +
  annotate("text", x = Inf, y = y_max_prof * 1.15,
           label = p_label_prof, color = "black", hjust = 1.05, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  coord_cartesian(clip = "off") +
  labs(tag = "I") +
  theme_pubr() +
  theme_global +
  theme_tag

# ========== J: Residential environment ==========
res_ttest_resid <- t.test(AgeAccelerationResidual ~ Residential_environment, data = pd)
p_val_resid <- res_ttest_resid$p.value
p_label_resid <- ifelse(p_val_resid < 0.001, "p < 0.001", paste0("p = ", round(p_val_resid, 3)))
y_max_resid <- max(pd$AgeAccelerationResidual, na.rm = TRUE)

pJ <- ggboxplot(pd, x = "Residential_environment", y = "AgeAccelerationResidual",
                ylab = "Age Acc", xlab = "Residential environment",
                add = "jitter", palette = "jco") +
  annotate("text", x = Inf, y = y_max_resid * 1.15,
           label = p_label_resid, color = "black", hjust = 1.05, vjust = 0, size = 4) +
  scale_y_continuous(expand = y_expand) +
  coord_cartesian(clip = "off") +
  labs(tag = "J") +
  theme_pubr() +
  theme_global +
  theme_tag

# ===========================================================
# 3. Combine and save all plots
# ===========================================================
cat("Combining all plots", sep = "\n")

# Reassignment of panel letters (A to J)
pA <- pA + labs(tag = "A")   # Age
pE <- pE + labs(tag = "B")   # BMI
pC <- pC + labs(tag = "C")   # Smoking
pB <- pB + labs(tag = "D")   # Pack-years
pD <- pD + labs(tag = "E")   # HbCO
pJ <- pJ + labs(tag = "F")   # Residential environment
pI <- pI + labs(tag = "G")   # Professional exposure
pF <- pF + labs(tag = "H")   # Fragmentation
pG <- pG + labs(tag = "I")   # OATS
pH <- pH + labs(tag = "J")   # Pregnancy outcome

# Layout: 2 columns x 5 rows
combined_plot <- grid.arrange(
  pA, pJ,    # Row 1: A=Age / F=Residential environment
  pE, pI,    # Row 2: B=BMI / G=Professional exposure
  pC, pF,    # Row 3: C=Smoking / H=Fragmentation
  pB, pG,    # Row 4: D=Pack-years / I=OATS
  pD, pH,    # Row 5: E=HbCO / J=Pregnancy outcome
  ncol = 2, nrow = 5,
  heights = unit(rep(1, 5), "null")
)

# Save RDS and PDF outputs
saveRDS(combined_plot, file.path(OUT_DIR, "FigEpigeneticAge_combined.rds"))

ggsave(file.path(OUT_DIR, "FigEpigeneticAge_combined.pdf"), combined_plot,
       width = 10, height = 14)

cat("Done!", sep = "\n")