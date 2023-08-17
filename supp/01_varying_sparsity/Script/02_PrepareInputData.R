# Title: Prepare Input Data for ScMaSigPro Evaluation
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(tidyverse))

# Set Paths relative to project
inPath <- "supp/01_varying_sparsity/data/input/sceObjects/"
outPath <- "supp/01_varying_sparsity/data/input/inputFiles/"

# Load filenames
fileNames <- list.files(inPath)

# Set Names
names(fileNames) <- fileNames

## Create a list of parameters
for (i in names(fileNames)) {
    
    #stop("Expected Stop")
    
    # String Split
    str_list <- strsplit(i, "_")
    
    # Name of the List
    cat(paste("\nPreparing for Zero-Inflation:", str_list[[1]][3]))
    
    # Load the object
    sce.obj <- readRDS(paste0(inPath, i))
    
    # Extract the Raw counts
    raw_counts <- as.matrix(sce.obj@assays@data@listData$counts)
    
    # Write the Raw Counts
    raw_counts_file_name <- paste0(outPath, "rawCounts/",str_list[[1]][1], "_",
                                   str_list[[1]][3], "_", str_list[[1]][4], "_", 
                                   str_list[[1]][5], "_",str_list[[1]][6], "_",
                                   str_remove(str_list[[1]][7], ".RDS"), "_rawCounts.tsv")
    dir.create(showWarnings = F, recursive = T, paste0(outPath, "rawCounts/"))

    # Extract Metadata
    cell_metadata <- as.data.frame(colData(sce.obj))
    
    # Select Columns of Interest
    cell_metadata <- cell_metadata[, c("Cell","Group", "Step"),drop = F]
    
    # Add Factor Design
    cell_metadata$Path1 <- ifelse(cell_metadata$Group == "Path1", 1,0)
    cell_metadata$Path2 <- ifelse(cell_metadata$Group == "Path2", 1,0)
    
    # Add Replicate Information
    cell_metadata <- as.data.frame(cell_metadata %>% 
                                       mutate(Replicate = data.table::rleid(
                                           Step, Path1, Path2)))
    
    # Select Columns of Interest
    cell_metadata <- cell_metadata[, c("Step","Replicate","Path1", "Path2"),drop = F]
    colnames(cell_metadata) <- c("Pseudotime","Replicate","Path1", "Path2")
    
    # Save data
    cell_meta_file_name <- paste0(outPath, "rawMeta/",str_list[[1]][1], "_",
                                  str_list[[1]][3], "_", str_list[[1]][4], "_", 
                                  str_list[[1]][5], "_",str_list[[1]][6], "_",
                                  str_remove(str_list[[1]][7], ".RDS"), "_rawMeta.tsv")
    dir.create(showWarnings = F, recursive = T, paste0(outPath, "rawMeta/"))
    
    # Extract the Ground Truth
    gene.meta <- as.data.frame(rowData(sce.obj))
    
    # Extract the columns of Interest
    gene.meta <- gene.meta[, !(colnames(gene.meta) %in% c("DEFacPath1","DEFacPath2",
                                                          "Gene", "gene_short_name",
                                                          "PathAExp","PathBExp")),]
    # Additional Annotation
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange == 0 & gene.meta$BaseToPathBChange == 0),"No_Change", "Change")
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange > 0 & gene.meta$BaseToPathBChange > 0),"UpBoth",gene.meta$status3)
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange < 0 & gene.meta$BaseToPathBChange < 0),"DownBoth",gene.meta$status3)
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange < 0 & gene.meta$BaseToPathBChange > 0),"Opposite",gene.meta$status3)
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange > 0 & gene.meta$BaseToPathBChange < 0),"Opposite",gene.meta$status3)
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange == 0 & gene.meta$BaseToPathBChange < 0),"DownPath2",gene.meta$status3)
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange == 0 & gene.meta$BaseToPathBChange > 0),"UpPath2",gene.meta$status3)
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange < 0 & gene.meta$BaseToPathBChange == 0),"DownPath1",gene.meta$status3)
    gene.meta$status3 <- ifelse((gene.meta$BaseToPathAChange > 0 & gene.meta$BaseToPathBChange == 0),"UpPath1",gene.meta$status3)
    
    # Subset and Rename
    gene.meta <- gene.meta[, c("BaseGeneMean", "status2", "status3")]
    colnames(gene.meta) <- c("base", "fc", "status")
    gene.meta$cobraTruth <- ifelse(gene.meta$status != "No_Change", 1,0)
    
    # Create Files
    gene_meta_file_name <- paste0(outPath, "rawGT/",str_list[[1]][1], "_",
                                  str_list[[1]][3], "_", str_list[[1]][4], "_", 
                                  str_list[[1]][5], "_",str_list[[1]][6], "_",
                                  str_remove(str_list[[1]][7], ".RDS"), "_rawGT.tsv")
    dir.create(showWarnings = F, recursive = T, paste0(outPath, "rawGT/"))
    
    # Sort Before Write
    cell_metadata <- cell_metadata[mixedorder(cell_metadata$Pseudotime),,drop = F]
    raw_counts <- raw_counts[, rownames(cell_metadata),drop = F]
    
    # Writting Files
    write.table(raw_counts, raw_counts_file_name, row.names = T, sep = "\t", col.names = NA)
    write.table(gene.meta, gene_meta_file_name, row.names = T, sep = "\t", col.names = NA)
    write.table(cell_metadata, cell_meta_file_name, row.names = T, sep = "\t", col.names = NA)
}

