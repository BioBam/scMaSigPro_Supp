# Title: Plot Evaluation Metrics
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(reshape2))

# Load Evaluation
evaluation.frame.zi <- as.data.frame(read.table("benchmarks/01_Sparsity/data/output/Performance.Table.tsv", header = T))
evaluation.frame.skew <- as.data.frame(read.table("benchmarks/02_Skewness/data/output/Performance.Table.tsv", header = T))
evaluation.frame.len <- as.data.frame(read.table("benchmarks/03_Different_Length/data/output/Performance.Table.tsv", header = T))

# Subset
evaluation.frame.zi.sub <- evaluation.frame.zi[(
    evaluation.frame.zi$parameter == "ZI" & evaluation.frame.zi$parameter.value == 60
    ),]
evaluation.frame.skew.sub <- evaluation.frame.skew[(
    evaluation.frame.skew$parameter == "Skew" & evaluation.frame.skew$parameter.value == 0.9
),]
evaluation.frame.len.sub <- evaluation.frame.len[(
    evaluation.frame.len$parameter == "Length" & evaluation.frame.len$parameter.value == "400_and_2600"
),]

# Combine
evaluation.frame <- rbind(evaluation.frame.zi.sub,
                          evaluation.frame.skew.sub,
                          evaluation.frame.len.sub)

# ROC Curve
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = parameter)) +
    geom_path(linewidth = 1, alpha = 0.7) +
    scale_color_manual(values = colorConesa(3)) +
    scale_x_continuous(breaks = seq(0, 0.10, 0.05), 
                       limits = c(0, 0.10)
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       limits = c(0, 1)
    ) +
    geom_point(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "ZI"),
        color = colorConesa(3)[3],
        size = 2) +
    geom_point(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Skew"),
        color = colorConesa(3)[2],
        size = 2) +
    geom_point(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Length"),
        color = colorConesa(3)[1],
        size = 2) +
    
    
    geom_point(
        data = subset(evaluation.frame, TPR > 0.9 & TPR < 0.92 & FPR < 0.01 & parameter == "Length"),
        color = colorConesa(3)[1],
        size = 2) +
    geom_text(data = subset(evaluation.frame, TPR > 0.9 & TPR < 0.92 & FPR < 0.01 & parameter == "Length"),
              aes(label = sprintf("%.2f", RSQ)), color = colorConesa(3)[1], hjust = -0.3, vjust = -0.6)+
    
    geom_point(
        data = subset(evaluation.frame, TPR >= 0.9 & FPR < 0.04 & parameter == "ZI"),
        color = colorConesa(3)[3],
        size = 2) +
    geom_text(data = subset(evaluation.frame, TPR > 0.9 & FPR < 0.04 & parameter == "ZI"),
              aes(label = sprintf("%.2f", RSQ)), color = colorConesa(3)[3], hjust = -0.3, vjust = -1)+
    
    geom_point(
        data = subset(evaluation.frame, TPR >= 0.9 & FPR < 0.04 & parameter == "Skew"),
        color = colorConesa(3)[2],
        size = 2) +
    geom_text(data = subset(evaluation.frame, TPR > 0.9 & FPR < 0.04 & parameter == "Skew"),
              aes(label = sprintf("%.2f", RSQ)), color = colorConesa(3)[2], hjust = -0.3, vjust = 1)+
    
    
    geom_text(data = subset(evaluation.frame, RSQ == 0.7 & parameter == "ZI"),
              aes(label = sprintf("%.2f", RSQ)), color = colorConesa(3)[3], hjust = -0.3, vjust = 0.6)+
    geom_text(data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Skew"),
              aes(label = sprintf("%.2f", RSQ)), color = colorConesa(3)[2], hjust = -0.3, vjust = 0.6)+
    geom_text(data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Length"),
              aes(label = sprintf("%.2f", RSQ)), color = colorConesa(3)[1], hjust = -0.3, vjust = 0.6)+
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
    geom_hline(yintercept = 0.8, colour = "darkgrey", linetype = "dotted") +
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
