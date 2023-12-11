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

# Load Plots
phate.plots <- readRDS("/supp_data/benchmarks/02_Skewness/simulated/png/01_skew_0_1.RDS")

# Load Evaluation
evaluation.frame <- read.table("Tables/02_Skew_Performance.Table.tsv", sep = "\t", header = T)
evaluation.frame.split <- read.table("Tables/02_SkewSplit_Performance.Table.tsv", sep = "\t", header = T)

# Plot all values against zero inflation
long_data <- melt(evaluation.frame, id.vars = c("RSQ", "parameter.value"), measure.vars = c("TPR", "FPR", "Accuracy", "F1_Score")) %>% as.data.frame()
long_split_data <- melt(evaluation.frame.split, id.vars = c("RSQ", "parameter.value"), measure.vars = c("TPR", "FPR", "Accuracy", "F1_Score")) %>% as.data.frame()

# Create Plot per parameter
performance.list <- lapply(unique(long_data$parameter.value), function(skew) {
  # get subset per parameter
  sub.df <- long_data[long_data$parameter.value == skew, ]

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
      # title = paste("Zero-inflation level", skew)
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
  if (skew == 0.1) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.5, xmax = 0.7, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (skew == 0) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.5, xmax = 0.7, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (skew == 1) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.6, xmax = 0.8, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (skew == 0.9) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.8, xmax = 1, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  }

  # Return
  return(performance.plot)
})
names(performance.list) <- paste("Skew", unique(long_data$parameter.value), sep = "_")

# Create Plot per parameter
performance.list.split <- lapply(unique(long_split_data$parameter.value), function(skew) {
  # get subset per parameter
  sub.df <- long_split_data[long_split_data$parameter.value == skew, ]

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
      # title = paste("Zero-inflation level", skew)
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
  if (skew == 0.1) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.45, xmax = 0.65, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (skew == 0) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.45, xmax = 0.65, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (skew == 1) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.55, xmax = 0.75, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (skew == 0.9) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.5, xmax = 0.7, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  }


  # Return
  return(performance.plot)
})
names(performance.list.split) <- paste("Skew", unique(long_split_data$parameter.value), sep = "_")

# Create
top <- ggarrange(
  phate.plots$skew_0,
  phate.plots$skew_0.1,
  phate.plots$skew_0.9,
  phate.plots$skew_1,
  labels = c("A.", "B.", "C.", "D."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "bottom"
)
bottom <- ggarrange(
  performance.list$Skew_0,
  performance.list$Skew_0.1,
  performance.list$Skew_0.9,
  performance.list$Skew_1,
  labels = c("E.", "F.", "G.", "H."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "none"
)
bottom.2 <- ggarrange(
  performance.list.split$Skew_0,
  performance.list.split$Skew_0.1,
  performance.list.split$Skew_0.9,
  performance.list.split$Skew_1,
  labels = c("I.", "J.", "K.", "L."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "bottom"
)
skewness <- ggarrange(top, bottom, bottom.2, nrow = 3, ncol = 1)
skewness

ggsave(skewness,
  filename = paste0("Figures/SuppData/02_Sim_0_to_1_skew_Performance.png"),
  dpi = 150, height = 10, width = 16
)
