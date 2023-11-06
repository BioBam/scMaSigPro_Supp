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
source(paste0(helpScriptsDir, "get_performance.R"))

# Load
load(file = paste0(dirPath, "scmp.obj.article.image.RData"))
    
# Extract the Row Data
row_data <- as.data.frame(
    rowData(scmp.obj@sce))[, c("gene_short_name", "status")]
    
# Get the gene-info
gene.change <- rownames(row_data[row_data$status != "No_Change", ])
gene.no.change <- rownames(row_data[row_data$status == "No_Change", ])
    
# Varying R2
r2_sequence_value <- seq(0.05, 0.95, 0.05)
    
# Get Performance
performance.measure <- get_performance(
    scmp_obj = scmp.obj,
    gene_change = gene.change,
    gene_no_change = gene.no.change,
    r2_sequence = r2_sequence_value
    )

# Write
write.table(performance.measure, paste0(dirPath, "Performance.Table.tsv"),
            sep = "\t",
            row.names = F, quote = F
)
