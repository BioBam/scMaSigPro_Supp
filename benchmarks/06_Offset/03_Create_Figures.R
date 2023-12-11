# Title: Plot Evaluation Metrics
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggpubr))

# Set Paths relative to project
outPath <- "Figures/SuppData/"

# Load Evaluation
evaluation.frame <- read.table("Tables/04_Normalization.Table.tsv", sep = "\t", header = T)

# Plot all values against zero inflation
long_data <- melt(evaluation.frame, id.vars = c("RSQ", "parameter.value"), measure.vars = c("TPR", "FPR", "Accuracy", "F1_Score")) %>% as.data.frame()

# Create Plot per parameter
performance.list <- lapply(unique(long_data$parameter.value), function(normType) {
  # get subset per parameter
  sub.df <- long_data[long_data$parameter.value == normType, ]

  # Create Performance Plot
  performance.plot <- ggplot(
    sub.df,
    aes(
      x = RSQ,
      y = value,
      group = interaction(parameter.value, variable),
      color = variable
    )
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.8) +
    scale_color_manual(values = colorConesa(6)) +
    labs(
      x = "Varying R-Square", y = "Performance Metric",
      title = paste(
        # "Normalization",
        normType
      )
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2), # Major x-axis breaks
      minor_breaks = seq(0, 1, by = 0.1) # Minor x-axis breaks
    ) +
    scale_y_continuous(
      breaks = c(seq(0, 1, by = 0.2), 0.05), # Major y-axis breaks
      minor_breaks = seq(0, 1, by = 0.1) # Minor y-axis breaks
    ) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "black") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25), # Customize major grid lines
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.1) # Customize minor grid lines
    )

  # Return
  return(performance.plot)
})
names(performance.list) <- paste("zi", unique(long_data$parameter.value), sep = "_")

# Create plot
norms <- ggarrange(
  performance.list$zi_FQ,
  performance.list$zi_libSize,
  performance.list$zi_logLibSize,
  performance.list$zi_rawCounts_Offset,
  performance.list$zi_rawCounts,
  performance.list$zi_SCT,
  labels = c("A.", "B.", "C.", "D.", "E.", "F."),
  common.legend = T, ncol = 3, nrow = 2,
  legend = "bottom"
)
norms

ggsave(norms,
  filename = paste0("Figures/SuppData/06_Offset_Performance.png"),
  dpi = 150, height = 8, width = 16
)

# ROC Curve
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = as.factor(parameter.value))) +
  # geom_text(data = subset(evaluation.frame, parameter.value == 60 & RSQ > 0.3 & RSQ <= 0.80), aes(label = sprintf("%.2f", RSQ)), color = "black", hjust = 1, vjust = 1.7) +
  geom_path(linewidth = 1, alpha = 0.7) +
  geom_point(data = subset(evaluation.frame, RSQ == 0.7), color = "black", shape = 2) +
  geom_point(data = subset(evaluation.frame, RSQ == 0.6), color = "blue", shape = 3) +
  scale_x_continuous(
    breaks = seq(0, 0.10, 0.05),
    limits = c(0, 0.10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    limits = c(0.0, 1)
  ) +
  scale_color_manual(values = colorConesa(9)) +
  labs(
    title = "ROC-curve for Normalization",
    subtitle = "Varying R2",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Amount of Zero-Inflation"
  ) +
  theme_classic(base_size = 15) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 0.01, colour = "grey") + # Highlighted the x-intercept of 0.01
  geom_vline(xintercept = 0.05, colour = "grey") +
  geom_vline(xintercept = 0.1, colour = "darkgrey") +
  guides(color = guide_legend(key_width = unit(3, "cm"), key_height = unit(4, "cm")))

print(roc)

# ggsave(roc,
#        filename = paste0("Figures/SuppData/01_Sim_10_to_90_ZI_ROC.png"),
#        dpi = 600, height = 8, width = 10
# )
