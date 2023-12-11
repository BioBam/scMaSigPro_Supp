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
phate.plots <- readRDS("/supp_data/benchmarks/03_Different_Length/simulated/png/03_len_400_1400.RDS")

# Load Evaluation
evaluation.frame <- read.table("Tables/03_Length_Performance.Table.tsv", sep = "\t", header = T)

# Plot all values against zero inflation
long_data <- melt(evaluation.frame, id.vars = c("RSQ", "parameter.value"), measure.vars = c("TPR", "FPR", "Accuracy", "F1_Score")) %>% as.data.frame()

# Create Plot per parameter
performance.list <- lapply(unique(long_data$parameter.value), function(len) {
  # get subset per parameter
  sub.df <- long_data[long_data$parameter.value == len, ]

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
  if (len == 100.29) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.6, xmax = 0.8, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (len == 1000.20) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.55, xmax = 0.75, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (len == 1500.15) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.4, xmax = 0.6, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  } else if (len == 500.25) {
    performance.plot <- performance.plot +
      annotate("rect", xmin = 0.55, xmax = 0.75, ymin = 0, ymax = 1, alpha = 0.3, fill = "lightgrey")
  }



  # Return
  return(performance.plot)
})
names(performance.list) <- paste("len", unique(long_data$parameter.value), sep = "_")

# Create
top <- ggarrange(
  phate.plots$len_100_2900,
  phate.plots$len_500_2500,
  phate.plots$len_1000_1000,
  phate.plots$len_1500_1500,
  labels = c("A.", "B.", "C.", "D."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "bottom"
)
bottom <- ggarrange(
  performance.list$len_100.29,
  performance.list$len_500.25,
  performance.list$len_1000.2,
  performance.list$len_1500.15,
  labels = c("E.", "F.", "G.", "G."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "bottom"
)
Length <- ggarrange(top, bottom, nrow = 2, ncol = 1)
Length

ggsave(Length,
  filename = paste0("Figures/SuppData/03_Sim_400_to_1400_length_Performance.png"),
  dpi = 150, height = 8, width = 16
)
