# Title: Evaluate the results of ScMaSigPro on simulated datasets with different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/01_Sparsity/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance.R"))

# Load names of files
dataSets <- list.files(paste0(dirPath))
dataSets <- dataSets[!(dataSets %in% c("Accuracy.png", "ROC.png", "Performance.Table.tsv"))]
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = "\\.", i = 4),
  ".RData"
)

# Zero-Inflation.evaluation
eval.list <- list()

# Set-up a for loop
for (i in names(dataSets)) {
  # Validation
  cat(paste("\nRunning for sparsity:", i))

  # Load
  load(file = paste0(dirPath, dataSets[i]))

  # Extract the Row Data
  row_data <- as.data.frame(
    rowData(scmp.obj@sce)
  )[, c("gene_short_name", "status")]

  # Get the gene-info
  gene.change <- rownames(row_data[row_data$status != "No_Change", ])
  gene.no.change <- rownames(row_data[row_data$status == "No_Change", ])

  # Varying R2
  r2_sequence_value <- seq(0.05, 0.95, 0.01)

  if (i == "80") {
    r2_sequence_value <- seq(0.05, 0.9, 0.05)
  }
  # Get Performance
  performance.measure <- get_performance(
    scmp_obj = scmp.obj,
    gene_change = gene.change,
    gene_no_change = gene.no.change,
    r2_sequence = r2_sequence_value
  )

  # Add Inflation
  performance.measure[["Zi"]] <- i

  # Add to list
  eval.list[[i]] <- performance.measure
}

# Combine
evaluation.frame <- do.call(rbind, eval.list)
write.table(evaluation.frame, paste0(dirPath, "Performance.Table.tsv"),
  sep = "\t",
  row.names = F, quote = F
)

# ROC
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = Zi)) +
    geom_point() +
    geom_text(data = subset(evaluation.frame, Zi == 60 & VARIABLE > 0.3 & VARIABLE <= 0.90), aes(label = sprintf("%.2f", VARIABLE)), color = "black", hjust = 1, vjust = 1.7) + 
    geom_path(linewidth = 1.5, alpha = 0.6) +
    scale_x_continuous(breaks = seq(0, 0.10, 0.01), 
                       limits = c(0, 0.10)
                       ) +
    scale_y_continuous(breaks = seq(0.5, 1, 0.1), 
                       limits = c(0.5, 1)
                       ) +
    scale_color_brewer(palette = "Set1", name = "Amount of Simulated Zero-Inflation") +
    labs(
        #title = "ROC-curve",
        #subtitle = "Varying levels of Coefficient of Determination (R-Square)",
        x = "False Positive Rate (1-Specificity)",
        y = "True Positive Rate (Sensitivity)"
    ) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
        legend.position = "bottom",
        legend.key.size = unit(3, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.7, "cm"),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    ) +
    geom_vline(xintercept = 0.01, colour = "grey") +  # Highlighted the x-intercept of 0.01
    geom_vline(xintercept = 0.05, colour = "grey") +
    geom_vline(xintercept = 0.1, colour = "darkgrey") +
    guides(color = guide_legend(key_width = unit(3, "cm"), key_height = unit(4, "cm")))

print(roc)



# Accuracy
acc <- ggplot(evaluation.frame, aes(
  x = VARIABLE, y = ACCURACY,
  color = Zi
)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.05)) +
  geom_path(linewidth = 1.5, alpha = 0.6) +
  scale_y_continuous(breaks = seq(0.7, 1, 0.05)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Accuracy Against Changing R Square",
    subtitle = "Red Dots: False Negatives",
    x = "Increasing Value for R-Square",
    y = "Accuracy"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(linewidth = 0.7, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.2, color = "grey", linetype = "dotted"),
    axis.text = element_text(size = rel(1.5)),
    legend.key.size = unit(4, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.position = "bottom",
    legend.key.height = unit(1, "cm"),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
  )


ggsave(acc,
  filename = paste0(dirPath, "Accuracy.png"),
  dpi = 600, height = 8, width = 8
)
ggsave(roc,
  filename = paste0(dirPath, "ROC.png"),
  dpi = 1200, height = 8, width = 8
)

sensitivity_specificity_plot <- ggplot(evaluation.frame, aes(x=VARIABLE)) +
geom_line(aes(y=TPR, color="Sensitivity", group = Zi)) +
geom_line(aes(y=TNR, color="Specificity", group = Zi)) +
scale_color_manual(values=c("blue", "green"), name="Metric") +
labs(x = "Threshold", y = "Rate") +
theme_minimal()
print(sensitivity_specificity_plot)


# Calculate precision
evaluation.frame$Precision <- with(evaluation.frame, TP / (TP + FP))

# Precision-Recall Plot
precision_recall_plot <- ggplot(evaluation.frame, aes(x=TPR, y=Precision, color=as.factor(Zi))) +
    geom_line() +
    geom_point() +
    scale_color_brewer(palette = "Set1", name = "Amount of Simulated Zero-Inflation") +
    labs(x = "Recall", y = "Precision") +
    theme_minimal()

print(precision_recall_plot)


library(ggplot2)
library(reshape2)

# Assuming 'evaluation.frame' is your dataframe and 'VARIABLE' is your varying parameter

# First, you might want to reshape your data frame to a long format.
long_data <- melt(evaluation.frame, id.vars = "VARIABLE", measure.vars = c("ACCURACY", "Precision", "F1", "Zi"))

# Now create the plot
p <- ggplot(long_data, aes(x = VARIABLE, y = value, color = variable)) +
    geom_path(aes(group = variable)) +
    facet_wrap(~variable, scales = "free_y") +
    labs(x = "VARIABLE", y = "Performance Metric") +
    theme_minimal()

print(p)



library(ggplot2)
library(reshape2)

# Reshape the data to a long format for ggplot
long_data <- melt(evaluation.frame, id.vars = c("VARIABLE", "Zi"), measure.vars = c("ACCURACY", "Precision", "F1", "FPR", "TPR", "FNR"))

# Create a line plot with VARIABLE on the x-axis and the performance metrics on the y-axis
# Facet by 'Zi' to create a separate plot for each zero inflation group
performance_plot <- ggplot(long_data, aes(x = VARIABLE, y = value, group = interaction(Zi, variable), color = variable)) +
    geom_line(linewidth = 1) +
    facet_wrap(~Zi, scales = "free_y") +
    labs(x = "VARIABLE", y = "Performance Metric") +
    theme_minimal() #+
    #scale_color_manual(
       # values = c("red", "green", "blue", "purple", "orange","pink" ), labels = c("Accuracy", "Precision", "F1", "FPR", "TPR","FNR"))

# Print the plot
print(performance_plot)

