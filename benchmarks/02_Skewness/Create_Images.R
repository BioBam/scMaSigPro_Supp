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
dirPath <- "benchmarks/02_Skewness/data/output/"

# Load Evaluation
evaluation.frame <- read.table(paste0(dirPath, "Performance.Table.tsv"), header = T)

# ROC Curve
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = as.factor(parameter.value))) +
    #geom_text(data = subset(evaluation.frame, parameter.value == 60 & RSQ > 0.3 & RSQ <= 0.80), aes(label = sprintf("%.2f", RSQ)), color = "black", hjust = 1, vjust = 1.7) + 
    geom_path(linewidth = 1, alpha = 0.7) +
    geom_point(data = subset(evaluation.frame, RSQ == 0.7), color = "black") +
    geom_point(data = subset(evaluation.frame, RSQ == 0.6), color = "blue") +
    scale_x_continuous(breaks = seq(0, 0.20, 0.05), 
                       limits = c(0, 0.20)
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       limits = c(0.0, 1)
    ) +
    scale_color_manual(values = colorConesa(9))+
    labs(
        title = "ROC-curve for Skewness",
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
long_data <- melt(evaluation.frame, id.vars = c("RSQ", "parameter.value"), measure.vars = c("TPR", "FPR", "Accuracy", "Precision", "Recall", "Specificity", "F1_Score"))

# Plot performance
performance <- ggplot(long_data, aes(x = RSQ, y = value, group = interaction(parameter.value, variable), color = variable)) +
    geom_line(linewidth = 0.6) + 
    geom_point(size = 0.8) +
    scale_color_manual(values = colorConesa(7)) +
    facet_wrap(~parameter.value, scales = "free_y", nrow = 3, ncol = 3, 
               labeller = labeller(parameter.value = function(x) paste("Skewness level", x))) +
    labs(x = "Varying R-Square", y = "Performance Metric",
         title = "Performance Metric for different levels of zero-inflation",
         color = "Measure") +
    scale_x_continuous(breaks = seq(0.1, 0.95, 0.2), limits = c(0, 0.95)) +
    theme_minimal(base_size = 20) + 
    theme(legend.position = "bottom")

print(performance)


# Save All RDS
saveRDS(performance, file = paste(dirPath, "Performance_plot.RDS"))
saveRDS(roc, file = paste(dirPath, "ROC_plot.RDS"))

ggsave(roc,
       filename = paste0("Images/Supp_Fig_3_A_Vary_Capture_Bias_ROC.png"),
       dpi = 600, height = 6, width = 6
)

ggsave(performance,
       filename = paste0("Images/Supp_Fig_3_B_Vary_Capture_Bias_Performance.png"),
       dpi = 600, height = 6, width = 6
)
