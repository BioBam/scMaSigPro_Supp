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
outPath <- "Figures/SuppData/"

# Load Evaluation
evaluation.frame <- read.table("Tables/01_ZI_Performance.Table.tsv", sep = "\t", header = T)

# ROC Curve
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = as.factor(parameter.value))) +
    #geom_text(data = subset(evaluation.frame, parameter.value == 60 & RSQ > 0.3 & RSQ <= 0.80), aes(label = sprintf("%.2f", RSQ)), color = "black", hjust = 1, vjust = 1.7) + 
    geom_path(linewidth = 1, alpha = 0.7) +
    geom_point(data = subset(evaluation.frame, RSQ == 0.7), color = "black", shape = 2) +
    geom_point(data = subset(evaluation.frame, RSQ == 0.6), color = "blue", shape = 3) +
    scale_x_continuous(breaks = seq(0, 0.20, 0.05), 
                       limits = c(0, 0.20)
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                      limits = c(0.0, 1)
    ) +
    scale_color_manual(values = colorConesa(9))+
    labs(
        title = "ROC-curve for Zero-Inflation",
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
    geom_vline(xintercept = 0.01, colour = "grey") +  # Highlighted the x-intercept of 0.01
    geom_vline(xintercept = 0.05, colour = "grey") +
    geom_vline(xintercept = 0.1, colour = "darkgrey") +
    guides(color = guide_legend(key_width = unit(3, "cm"), key_height = unit(4, "cm")))

print(roc)

# Plot all values against zero inflation
long_data <- melt(evaluation.frame, id.vars = c("RSQ", "parameter.value"), measure.vars = c("TPR", "FPR", "Accuracy", "Recall", "Specificity", "F1_Score"))

# Plot performance
performance <- ggplot(long_data, aes(x = RSQ, y = value, group = interaction(parameter.value, variable), color = variable)) +
    geom_line(linewidth = 0.6) + 
    geom_point(size = 0.8) +
    scale_color_manual(values = colorConesa(7)) +
    facet_wrap(~parameter.value, scales = "free_y", nrow = 3, ncol = 3, 
               labeller = labeller(parameter.value = function(x) paste("Zero-inflation level", x))) +
    labs(x = "Varying R-Square", y = "Performance Metric",
         title = "Performance Metric for different levels of zero-inflation",
         color = "Measure") +
    scale_x_continuous(breaks = seq(0.1, 0.95, 0.2), limits = c(0, 0.95)) +
    theme_minimal(base_size = 15) + 
    theme(legend.position = "bottom")

print(performance)

ggsave(roc,
       filename = paste0("Figures/SuppData/01_Sim_10_to_90_ZI_ROC.png"),
       dpi = 600, height = 8, width = 10
)

ggsave(performance,
       filename = paste0("Figures/SuppData/01_Sim_10_to_90_ZI_Performance.png"),
       dpi = 600, height = 8, width = 10
)
