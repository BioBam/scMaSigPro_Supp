# Title: Plot Evaluation Metrics
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggpubr))

# Set paths
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
rdsPath <- paste0(base_string, "benchmarks/01_Sparsity/sim/")
imgPath <- paste0(base_string, "benchmarks/01_Sparsity/img/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Load Plots
umap.plots <- readRDS(paste0(imgPath, "01_Zi_60_90.RDS"))

# Load Evaluation
evaluation.frame <- read.table(paste0(tabPath, "01_ZI_Performance.Table.tsv"), sep = "\t", header = T)

# Plot all values against zero inflation
long_data <- melt(evaluation.frame,
  id.vars = c("RSQ", "parameter.value"),
  measure.vars = c("TPR", "FPR", "Accuracy", "F1_Score")
) %>% as.data.frame()

# Create Plot per parameter
performance.list <- lapply(unique(long_data$parameter.value), function(zi) {
  # get subset per parameter
  sub.df <- long_data[long_data$parameter.value == zi, ]

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
      # title = paste("Zero-inflation level", zi)
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


  # Add Annotation
  if (zi == 60) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.45, xmax = 0.65, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (zi == 70) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.4, xmax = 0.6, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (zi == 80) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.35, xmax = 0.55, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (zi == 90) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.3, xmax = 0.5, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  }

  print(paste("finisjed for ", zi))
  # Return
  return(performance.plot)
})
names(performance.list) <- paste("zi", unique(long_data$parameter.value), sep = "_")

# Create
top <- ggarrange(
  umap.plots$sparsity_60,
  umap.plots$sparsity_70,
  umap.plots$sparsity_80,
  umap.plots$sparsity_90,
  labels = c("A.", "B.", "C.", "D."),
  common.legend = T,
  ncol = 4, nrow = 1,
  legend = "bottom"
)
bottom <- ggarrange(
  performance.list$zi_60,
  performance.list$zi_70,
  performance.list$zi_80,
  performance.list$zi_90,
  labels = c("E.", "F.", "G.", "H."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "bottom"
)
zero_inflation <- ggarrange(top, bottom, nrow = 2, ncol = 1)
zero_inflation

# Save the plot
ggsave(
  plot = zero_inflation,
  filename = paste0(figPath_hd, "01_Sim_60_to_90_ZI_Performance.png"),
  dpi = 600, height = 8, width = 16
)
ggsave(
  plot = zero_inflation,
  filename = paste0(figPath_lr, "01_Sim_60_to_90_ZI_Performance.png"),
  dpi = 150, height = 8, width = 16
)

#
# # ROC Curve
# roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = as.factor(parameter.value))) +
#   # geom_text(data = subset(evaluation.frame, parameter.value == 60 & RSQ > 0.3 & RSQ <= 0.80), aes(label = sprintf("%.2f", RSQ)), color = "black", hjust = 1, vjust = 1.7) +
#   geom_path(linewidth = 1, alpha = 0.7) +
#   geom_point(data = subset(evaluation.frame, RSQ == 0.7), color = "black", shape = 2) +
#   geom_point(data = subset(evaluation.frame, RSQ == 0.6), color = "blue", shape = 3) +
#   scale_x_continuous(
#     breaks = seq(0, 0.10, 0.05),
#     limits = c(0, 0.10)
#   ) +
#   scale_y_continuous(
#     breaks = seq(0, 1, 0.1),
#     limits = c(0.0, 1)
#   ) +
#   scale_color_manual(values = colorConesa(9)) +
#   labs(
#     title = "ROC-curve for Zero-Inflation",
#     subtitle = "Varying R2",
#     x = "False Positive Rate (1-Specificity)",
#     y = "True Positive Rate (Sensitivity)",
#     color = "Amount of Zero-Inflation"
#   ) +
#   theme_classic(base_size = 15) +
#   theme(
#     panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
#     panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
#     legend.position = "bottom"
#   ) +
#   geom_vline(xintercept = 0.01, colour = "grey") + # Highlighted the x-intercept of 0.01
#   geom_vline(xintercept = 0.05, colour = "grey") +
#   geom_vline(xintercept = 0.1, colour = "darkgrey") +
#   guides(color = guide_legend(key_width = unit(3, "cm"), key_height = unit(4, "cm")))
#
# print(roc)
#
# # ggsave(roc,
# #        filename = paste0("/supp_data/Figu res/SuppData/01_Sim_10_to_90_ZI_ROC.png"),
# #        dpi = 600, height = 8, width = 10
# # )
