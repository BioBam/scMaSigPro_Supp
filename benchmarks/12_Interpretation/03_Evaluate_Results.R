# Title: Evaluate the results of ScMaSigPro on simulated datasets with of lenEq
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/10_Split_bin_Parameter/output/"
outPath <- "Tables/"
helpScriptsDir <- "R_Scripts/helper_function/"


# Load helper functions
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
  str_remove(dataSets, pattern = "scmp.obj.skew."),
  ".RData"
)

# Zero-Inflation.evaluation
eval.list <- list()

# Set-up a for loop
for (i in names(dataSets)) {
  # Validation
  cat(paste("\nRunning for Skewness:", i))

  # Load
  load(file = paste0(inPath, dataSets[i]))

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
  performance.measure <- as.data.frame(get_performance_ROCR(
    scmpObj = scmp.obj,
    groundTruth = groundTruth,
    r2_sequence = seq(0.00, 0.95, 0.05),
    include_influ = TRUE
  ))

  # Add to list
  performance.measure[["parameter"]] <- "Skew"
  performance.measure[["parameter.value"]] <- i
  eval.list[[i]] <- performance.measure
}

# Combine
evaluation.frame <- do.call(rbind, eval.list)

# Write
write.table(evaluation.frame, paste0(outPath, "011_Skew_SplitBins_Performance.Table.tsv"),
  sep = "\t",
  row.names = F, quote = F
)
