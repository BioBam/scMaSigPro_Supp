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
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = as.factor(skew))) +
    geom_point() +
    geom_text(data = subset(evaluation.frame, skew == 60 & VARIABLE > 0.3 & VARIABLE <= 0.90), aes(label = sprintf("%.2f", VARIABLE)), color = "black", hjust = 1, vjust = 1.7) + 
    geom_path(linewidth = 1, alpha = 0.7) +
    scale_x_continuous(breaks = seq(0, 0.1, 0.01), 
                       limits = c(0, 0.1)
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       limits = c(0, 1)
    ) +
    scale_color_manual(values = colorConesa(10))+
    labs(
        title = "ROC-curve for Skewness",
        subtitle = "Varying R2",
        x = "False Positive Rate (1-Specificity)",
        y = "True Positive Rate (Sensitivity)",
        color = "Skewness"
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


# Accuracy
acc <- ggplot(evaluation.frame, aes(
    x = VARIABLE, y = ACCURACY,
    color = as.factor(skew)
)) +
    geom_point() +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), 
                       limits = c(0, 0.9)
    ) +
    scale_y_continuous(breaks = seq(0.5, 1, 0.1),
                       limits = c(0.5, 1)
    ) +
    geom_path(linewidth = 1, alpha = 0.6) +
    scale_color_manual(values = colorConesa(10))+
    labs(
        title = "Accuracy Against R Square",
        x = "Increasing Value for R-Square",
        y = "Accuracy",
        color = "Amount of Skewness"
    ) +
    theme_classic(base_size = 15) +
    theme(
        panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
        legend.position = "bottom",
    )

print(acc)


# Plot all values against zero inflation
long_data <- melt(evaluation.frame, id.vars = c("VARIABLE", "skew"), measure.vars = c("ACCURACY", "PRECISION", "FPR", "TPR", "FNR"))

# Plot performance
performance <- ggplot(long_data, aes(x = VARIABLE, y = value, group = interaction(skew, variable), color = variable)) +
    geom_line(linewidth = 0.6) + 
    geom_point(size = 0.8) +
    scale_color_manual(values = colorConesa(5)) +
    facet_wrap(~skew, scales = "free_y", nrow = 2, ncol = 5, 
               labeller = labeller(Zi = function(x) paste("Skewness level", x))) +
    labs(x = "Varying R-Square", y = "Performance Metric",
         title = "Performance Metric for different levels of Skewness",
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
