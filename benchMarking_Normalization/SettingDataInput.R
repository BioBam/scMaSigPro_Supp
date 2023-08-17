# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtools))

# Load Prefix
parentFix <- "benchMarking_Normalization/" 
scriptFix <- paste0("old_Scripts/")
inPrefix <- paste0(parentFix, "data/output/SimObjectsSCE/")
outPrefix <- paste0(parentFix, "data/output/")

# Load File names
#all_file_names <- list.files("benchMarking_Compression/Automated/rdsObjects/")
all_file_names <- list.files(inPrefix)
all_file_names <- str_remove(all_file_names, ".RDS")
split_list <- str_split(all_file_names, "_")
names(split_list) <- list.files(inPrefix)
#split_list <- split_list[names(split_list) %in% "zi_shape_-0.2.RDS"]

# Run for Loop one-by-one
for (i in names(split_list)){
    
    # Set Variables
    file_name <- i
    
    # Name of Test
    test_name <- split_list[[i]][1]
    
    # Test_mid
    test_mid <- split_list[[i]][2]
    
    # Test Value
    test_value <- split_list[[i]][3]
    
    # Validity Run
    cat(paste("\n\nRunning for", test_name, "at", test_value, sep = " "))
    
    # Load data
    sce_obj_file_name <- paste0(inPrefix, i) 
    sce.obj <- readRDS(sce_obj_file_name)
    #print(sce.obj)
    
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
    dir_path <- paste0(outPrefix, "RawMeta/")
    dir.create(dir_path, showWarnings = F)
    cell_meta_file_name <- paste0(dir_path, test_name, "_",test_mid, "_", test_value, "_", "RawMeta.tsv")
    cell_metadata <- cell_metadata[mixedorder(cell_metadata$Pseudotime),, drop = F]
    write.table(cell_metadata, cell_meta_file_name, row.names = T, sep = "\t", col.names = NA)
    
    # Extract the Raw counts
    raw_counts <- as.matrix(sce.obj@assays@data@listData$TrueCounts)
    
    # Write the Raw Counts
    dir_path <- paste0(outPrefix, "RawCounts/")
    dir.create(dir_path, showWarnings = F)
    raw_counts_file_name <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "TrueCounts.tsv")
    raw_counts <- raw_counts[, rownames(cell_metadata),drop = F]
    write.table(raw_counts, raw_counts_file_name, row.names = T, sep = "\t", col.names = NA)
    
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
    dir_path <- paste0(outPrefix, "GroundTruth/")
    dir.create(dir_path, showWarnings = F)
    gene_meta_file_name <- paste0(dir_path,test_name, "_", test_mid, "_",test_value, "_", "GroundTruth.tsv")
    write.table(gene.meta, gene_meta_file_name, row.names = T, sep = "\t", col.names = NA)
}  

cat(paste0("\n","Script Finished"))