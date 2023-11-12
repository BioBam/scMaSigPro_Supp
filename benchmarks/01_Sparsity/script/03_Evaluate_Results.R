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
#source(paste0(helpScriptsDir, "get_performance.R"))
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Load names of files
dataSets <- list.files(paste0(dirPath))
dataSets <- dataSets[!(dataSets %in% c("Accuracy.png", "ROC.png", "Performance.Table.tsv"))]
names(dataSets) <- str_remove(
  str_remove(dataSets, pattern = "scmp.obj.sparsity."),
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
      rowData(scmp.obj@sce))[, c("gene_short_name", "status")]
  
  # Set binary labels
  gene.change <- rep(1, length(rownames(row_data[row_data$status != "No_Change", ])))
  gene.no.change <-rep(0, length(rownames(row_data[row_data$status == "No_Change", ])))
  
  # Add names
  names(gene.change) <- rownames(row_data[row_data$status != "No_Change", ])
  names(gene.no.change) <- rownames(row_data[row_data$status == "No_Change", ])
  
  # Ground truth
  groundTruth <- c(gene.change, gene.no.change)
  
  # Get Performance
  performance.measure <- as.data.frame(get_performance_ROCR(
      scmpObj = scmp.obj,
      groundTruth = groundTruth,
      r2_sequence = seq(0.00, 0.95, 0.05),
      include_influ = TRUE
  ))
  
  # Add to list
  performance.measure[["parameter"]] <- "ZI"
  performance.measure[["parameter.value"]] <- as.numeric(i)
  eval.list[[i]] <- performance.measure
}

# Combine
evaluation.frame <- do.call(rbind, eval.list)

# Write
write.table(evaluation.frame, paste0(dirPath, "Performance.Table.tsv"),
  sep = "\t",
  row.names = F, quote = F
)
