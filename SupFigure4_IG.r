#!/usr/bin/env Rscript
# Methylation of 24 ICR CpGs vs age, in this study and in the public datasets
# (guide script).
#
# Input: data_plot_icr_cpgs.rds, a data frame with one row per sample and CpG:
#   CpG, dataset ("This study" or the public dataset), Age, Beta, facet_label.

library(dplyr)
library(ggplot2)

DATA_DIR <- "path_to_data"
OUT_DIR  <- file.path(DATA_DIR, "ICR")
dir.create(OUT_DIR, showWarnings = FALSE)

data_plot <- readRDS(file.path(DATA_DIR, "data_plot_icr_cpgs.rds"))

# CpGs in the order of the figure panels (6 rows x 4 columns)
target_cpgs <- c(
  "cg19642877", "cg20622142", "cg26516510", "cg27263998",
  "cg04104132", "cg07595203", "cg08359167", "cg08446215",
  "cg19643803", "cg22940816", "cg23870378", "cg24089209",
  "cg02125271", "cg02305723", "cg04072648", "cg08560373",
  "cg11562309", "cg11985632", "cg13797824", "cg14873490",
  "cg21667116", "cg23156962", "cg25407198", "cg26708559"
)

# "This study" first, then the other datasets in alphabetical order
other_datasets <- sort(setdiff(unique(data_plot$dataset), "This study"))
dataset_levels <- c("This study", other_datasets)

data_filtered <- data_plot %>%
  filter(CpG %in% target_cpgs) %>%
  mutate(
    CpG = factor(CpG, levels = target_cpgs),
    dataset = factor(dataset, levels = dataset_levels)
  ) %>%
  arrange(match(CpG, target_cpgs))

# Panel labels in the order of target_cpgs
facet_order <- data_filtered %>%
  distinct(CpG, facet_label) %>%
  slice(match(target_cpgs, CpG)) %>%
  pull(facet_label)

data_filtered <- data_filtered %>%
  mutate(facet_label = factor(facet_label, levels = facet_order))

# Drawn separately so that "This study" is displayed on top
data_others <- data_filtered %>% filter(dataset != "This study")
data_this   <- data_filtered %>% filter(dataset == "This study")

FigICR <- ggplot(data_filtered, aes(x = Age, y = Beta)) +
  geom_smooth(aes(color = dataset), method = "lm", se = TRUE, alpha = 0.15,
              linewidth = 0.8, na.rm = TRUE, show.legend = FALSE) +
  geom_point(
    data = data_others,
    aes(x = Age, y = Beta, color = dataset),
    size = 1.3, alpha = 0.7, na.rm = TRUE
  ) +
  geom_point(
    data = data_this,
    aes(x = Age, y = Beta, color = dataset),
    size = 1.3, alpha = 0.7, na.rm = TRUE
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3, linetype = 0))) +
  facet_wrap(~ facet_label, ncol = 4, scales = "free_y") +
  labs(
    title = NULL,
    x = "Age (years)",
    y = "Beta value",
    color = "Dataset"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom"
  )

print(FigICR)

# Figure size (inches): adjust to the journal requirements
ggsave(file.path(OUT_DIR, "FigICR_age.pdf"), FigICR, width = 10, height = 12)