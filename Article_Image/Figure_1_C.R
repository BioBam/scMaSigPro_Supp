# Title: Plot Evaluation Metrics
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(reshape2))

# Set Paths relative to project
dirPath <- "Article_Image/data/output/"

# Load Evaluation
evaluation.frame <- as.data.frame(read.table(paste0(dirPath, "Performance.Table.tsv"), header = T))

# ROC Curve
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR)) +
    geom_point(data = subset(evaluation.frame,RSQ %in% c(0.6, 0.65, 0.7, 0.75)), color = colorConesa(4, reverse = T, palette = "hot")) +
    geom_text(data = subset(evaluation.frame,RSQ %in% c(0.6, 0.65, 0.7, 0.75)), aes(label = sprintf("%.2f", RSQ)), color = "black", hjust = 1, vjust = -0.6)+
    geom_path(linewidth = 1, alpha = 1, color = colorConesa(1)) +
    scale_x_continuous(breaks = seq(0, 0.10, 0.05), 
                       limits = c(0, 0.10)
    ) +
    scale_y_continuous(breaks = seq(0.3, 1, 0.1),
                       limits = c(0.3, 1)
    ) +
    labs(
        title = "ROC-Curve Simulated Data",
        subtitle = "Varying R-Square",
        x = "False Positive Rate (1-Specificity)",
        y = "True Positive Rate (Sensitivity)") +
    theme_classic(base_size = 15) +
    theme(
        panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
        legend.position = "bottom"
    ) +
    geom_vline(xintercept = 0.01, colour = "grey") +  # Highlighted the x-intercept of 0.01
    geom_vline(xintercept = 0.05, colour = "grey") +
    geom_vline(xintercept = 0.1, colour = "darkgrey") +
    guides(color = guide_legend(key_width = unit(3, "cm"), key_height = unit(4, "cm")))

print(roc)

# Plot all values against zero inflation
long_data <- melt(evaluation.frame, id.vars = c("RSQ"), measure.vars = c("Accuracy", "Precision", "FPR", "TPR", "Recall", "F1_Score"))

# Plot performance
performance <- ggplot(long_data, aes(x = RSQ, y = value, group = interaction(variable), color = variable)) +
    geom_line(linewidth = 0.6) + 
    geom_point(size = 0.8) +
    scale_color_manual(values = colorConesa(6)) +
    labs(x = "Varying R-Square", y = "Performance Metric",
         title = "Performance Metric for different levels of zero-inflation",
         color = "Measure") +
    scale_x_continuous(breaks = seq(0.1, 0.95, 0.2), limits = c(0, 0.95)) +
    theme_minimal(base_size = 10) + 
    theme(legend.position = "bottom")

print(performance)

# Save All RDS
saveRDS(performance, file = paste(dirPath, "Performance_plot.RDS"))
saveRDS(acc, file = paste(dirPath, "Accuracy_Plot.RDS"))
saveRDS(roc, file = paste(dirPath, "ROC_plot.RDS"))

# Save Images
ggsave(acc,
       filename = paste0(dirPath, "Accuracy.png"),
       dpi = 600, height = 6, width = 6
)
ggsave(roc,
       filename = paste0(dirPath, "ROC.png"),
       dpi = 600, height = 6, width = 6
)

ggsave(performance,
       filename = paste0(dirPath, "Performance.png"),
       dpi = 600, height = 6, width = 6
)
