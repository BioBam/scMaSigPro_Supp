# Title: Prepare Input Data for ScMaSigPro Evaluation
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(tidyverse))

# Load Functions
suppressPackageStartupMessages(source("old_Scripts/entropy_discretize.R"))
suppressPackageStartupMessages(source("old_Scripts/make_bulk_design.R"))
suppressPackageStartupMessages(source("old_Scripts/make_bulk_counts.R"))
suppressPackageStartupMessages(source("old_Scripts/create_range.R"))

# Set Paths relative to project
inPath <- "supp/07_Different_Length_of_Pseudotime/data/input/inputFiles/"
outPath <- "supp/07_Different_Length_of_Pseudotime/data/input/PseudobulkInput/"

# Load filenames
fileNames <- list.files(paste0(inPath,"rawMeta"))

# Split String
fileNames <- sapply(str_split(fileNames, "_"), function(x){return(paste(c(x[1], x[2], x[3], x[4]), collapse = "_"))})

# Set Names
names(fileNames) <- fileNames

## Create a list of parameters
for (i in names(fileNames)) {
    
    # Load all the files
    rawMeta <- read.table(paste0(inPath, "rawMeta/",
                                 paste0(i, "_rawMeta.tsv")),
                          row.names = 1, header = 1)
    rawCounts <- as.matrix(read.table(paste0(inPath, "rawCounts/",
                                 paste0(i, "_rawCounts.tsv")),
                          row.names = 1, header = 1))
    rawGT <- read.table(paste0(inPath, "rawGT/",
                                 paste0(i, "_rawGT.tsv")),
                          row.names = 1, header = 1)
    
    # String Split
    str_list <- str_split(i, "_")
    
    # Name of the List
    cat(paste("\nPreparing for Length:", str_list[[1]][1], str_list[[1]][2], str_list[[1]][3],
              str_list[[1]][4]))
    
    # Add Compression Columns for Pseudobulking
    pbMeta <- as.data.frame(entropy_discretize(rawMeta, time_col = "Pseudotime"))
    
    # bulk-the cell Metadata
    bulk_cell_metadata <- make_pseudobulk_design(design.file = pbMeta,
                                                 paths.vector = c("Path1", "Path2"),
                                                 binnedCol = "binnedTime")
    
    # Sum, Mean, Average, of Counts 
    bulk_counts_list <- make_bulk_counts(counts = rawCounts, cluster_member_col = "cluster.members", 
                                         bin_col = "bin", pseudo_bulk_profile = bulk_cell_metadata,
                                         round = T)
    
    # Extract Sum Counts
    bulk_counts <- bulk_counts_list[["pseudo_bulk_counts_sum"]]
    
    # Extract MetaData
    bulk_cell_metadata <- bulk_cell_metadata[,c("binnedTime", "Path1", "Path2"), drop = F]
    bulk_cell_metadata <- as.data.frame(bulk_cell_metadata %>% 
                                            mutate(Replicate = data.table::rleid(
                                                binnedTime, Path1, Path2)))
    bulk_cell_metadata <- bulk_cell_metadata[, c("binnedTime", "Replicate", "Path1", "Path2"),drop = F]
    colnames(bulk_cell_metadata) <- c("Pooled_Pseudotime", "Replicate", "Path1", "Path2")
    
    # Create an Object
    pseudobulk.data.input <- list(pbCounts = as.matrix(bulk_counts),
                                  pbMeta = bulk_cell_metadata,
                                  gt = rawGT)
    
    # Write
    dir.create(outPath,showWarnings = F, recursive = T)
    saveRDS(pseudobulk.data.input,
            paste0(outPath, str_list[[1]][1], "_", str_list[[1]][2],
                   "_", str_list[[1]][3],
                   "_", str_list[[1]][4], ".RDS"))
    }

