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
rdsPath <- paste0(base_string, "benchmarks/03_Different_Length/sim/")
imgPath <- paste0(base_string, "benchmarks/03_Different_Length/img/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Load Plots
umap.plots <- readRDS(paste0(imgPath, "03_len_100_1500.RDS"))

# Load Evaluation
evaluation.frame <- read.table(paste0(tabPath, "03_Length_Performance.Table.tsv"), sep = "\t", header = T)

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
  umap.plots$len_100_2900,
  umap.plots$len_500_2500,
  umap.plots$len_1000_1000,
  umap.plots$len_1500_1500,
  labels = c("A.", "B.", "C.", "D."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "bottom"
)
bottom <- ggarrange(
  performance.list$len_100.29,
  performance.list$len_500.25,
  performance.list$len_1000.2,
  performance.list$len_1500.15,
  labels = c("E.", "F.", "G.", "H."),
  common.legend = T, ncol = 4, nrow = 1,
  legend = "bottom"
)
Length <- ggarrange(top, bottom, nrow = 2, ncol = 1)
Length

# Save the plot
ggsave(
  plot = Length,
  filename = paste0(figPath_hd, "03_Sim_100_to_1500_length_Performance.png"),
  dpi = 600, height = 8, width = 16
)
ggsave(
  plot = Length,
  filename = paste0(figPath_lr, "03_Sim_100_to_1500_length_Performance.png"),
  dpi = 150, height = 8, width = 16
)
