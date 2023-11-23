# Title: Evaluate the results of ScMaSigPro on simulated datasets with different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "Article_Image/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Load
load(file = paste0(dirPath, "scmp.obj.article.image.RData"))

# Extract the Row Data
row_data <- as.data.frame(
  rowData(scmp.obj@sce)
)[, c("gene_short_name", "status")]

# Set binary labels
gene.change <- rep(1, length(rownames(row_data[row_data$status != "No_Change", ])))
gene.no.change <- rep(0, length(rownames(row_data[row_data$status == "No_Change", ])))

# Add names
names(gene.change) <- rownames(row_data[row_data$status != "No_Change", ])
names(gene.no.change) <- rownames(row_data[row_data$status == "No_Change", ])

# Ground truth
groundTruth <- c(gene.change, gene.no.change)

# Get Performance
performance.measure <- get_performance_ROCR(
  scmpObj = scmp.obj,
  groundTruth = groundTruth,
  r2_sequence = seq(0.00, 0.95, 0.05),
  include_influ = FALSE
)

# Write
write.table(performance.measure, paste0(dirPath, "Performance.Table.tsv"),
  sep = "\t",
  row.names = F, quote = F
)
