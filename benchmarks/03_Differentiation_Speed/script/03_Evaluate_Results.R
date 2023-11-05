# Title: Evaluate the results of ScMaSigPro on simulated datasets with of lenEq
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/03_Differentiation_Speed/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance.R"))

# Load names of files
dataSets <- list.files(paste0(dirPath))
dataSets <- dataSets[!(dataSets %in% c("Accuracy.png", "ROC.png", "Performance.Table.tsv"))]
names(dataSets) <- str_remove(
  str_remove(dataSets, pattern = "scmp.obj.Arm."),
  ".RData"
)

# Zero-Inflation.evaluation
eval.list <- list()

# Set-up a for loop
for (i in names(dataSets)) {
  # Validation
  cat(paste("\nRunning for lenEq:", i))

  # Load
  load(file = paste0(dirPath, dataSets[i]))

  # Extract the Row Data
  row_data <- as.data.frame(
    rowData(scmp.obj@sce)
  )[, c("gene_short_name", "status")]

  # Get the gene-info
  gene.change <- rownames(row_data[row_data$status != "No_Change", ])
  gene.no.change <- rownames(row_data[row_data$status == "No_Change", ])


  r2_sequence_value <- seq(0.05, 0.95, 0.05)

  # Get Performance
  performance.measure <- get_performance(
    scmp_obj = scmp.obj,
    gene_change = gene.change,
    gene_no_change = gene.no.change,
    r2_sequence = r2_sequence_value
  )

  # Add Inflation
  performance.measure[["Arm"]] <- i

  # Add to list
  eval.list[[i]] <- performance.measure
}

# Combine
evaluation.frame <- do.call(rbind, eval.list)

# Write
write.table(evaluation.frame, paste0(dirPath, "Performance.Table.tsv"),
  sep = "\t",
  row.names = F, quote = F
)