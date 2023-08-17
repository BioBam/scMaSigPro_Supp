# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtools))

# Load Prefix
parentFix <- "benchMarking_Normalization/" 
scriptFix <- paste0("old_Scripts/")
inPrefix <- paste0(parentFix, "data/output/")
outPrefix <- paste0(parentFix, "data/output/")

# Load Custim Functions
load <- function(){
    source(paste0(scriptFix, "plot_simulations().R"))
    source(paste0(scriptFix, "add_gene_anno().R"))
    source(paste0(scriptFix, "calc_bin_size.R"))
    source(paste0(scriptFix, "calcNormCounts.R"))
    source(paste0(scriptFix, "entropy_discretize.R"))
    source(paste0(scriptFix, "make_bulk_design.R"))
    source(paste0(scriptFix, "make_bulk_counts.R"))
    source(paste0(scriptFix, "create_range.R"))
}

# Load File names
all_file_names <- list.files(inPrefix)
# all_file_names <- all_file_names[all_file_names %in% c(
#     "CLRNorm", "CPMNorm", "FQNorm", "LogNorm",
#     "RawCounts", "RcNorm", "SctNorm", "RawMeta")]
all_file_names <- all_file_names[all_file_names %in% c("RawCounts")]
names(all_file_names) <- all_file_names


# Load files in the list
count_name_list <- list()
for (i in all_file_names){
    count_table_names <- list.files(paste0(inPrefix, i, "/"))
    count_table_names <- count_table_names[!grepl(pattern = ".png", count_table_names)]
    count_name_list[[i]] <- count_table_names
}

count_name_list[["RawCounts"]] <- count_name_list[["RawCounts"]][count_name_list[["RawCounts"]] == "zi_shape_-0.2_TrueCounts.tsv"]

load()

# Run per count
for (i in names(count_name_list)){
    
    print(paste0("Running: ", i))
    file_names <- count_name_list[[i]]
    
    for (j in file_names){
        print(paste0("Reading: ", j))
        count_table <- as.matrix(read.table(
            file = paste0(inPrefix, i,"/", j),
            header = T, sep = "\t", row.names = 1))
        
        # Load the corresponding Meta
        cell_meta_name <-  str_split(j, "_")
        cell_meta_name <- paste0(outPrefix, "RawMeta/",
                                 paste0(cell_meta_name[[1]][1], "_",
                                        cell_meta_name[[1]][2], "_",
                                        cell_meta_name[[1]][3], "_RawMeta.tsv"))
        
        print(paste0("Reading: ", str_sub(cell_meta_name, start = 48)))
                     
        rawMeta <- as.data.frame(read.table(
            file = cell_meta_name,
            header = T, sep = "\t", row.names = 1))
        
        
        # Pseudobulk the Metadata
        #Add Compression Columns for Pseudobulking
        bulk_discretize <- entropy_discretize(
            rawMeta, time_col = "Pseudotime", drop.fac = 0.6)
        
        # bulk-the cell Metadata
        bulk_cell_metadata <- make_pseudobulk_design(design.file = bulk_discretize,
                                                         paths.vector = c("Path1", "Path2"),
                                                         binnedCol = "binnedTime")
        
        bulk_counts_list <- make_bulk_counts(counts = count_table, round = T,
                                        cluster_member_col = "cluster.members",
                                        bin_col = "bin", 
                                        pseudo_bulk_profile = bulk_cell_metadata)
        
        bulk_cell_metadata <- bulk_cell_metadata[,c("binnedTime", "Path1", "Path2"), drop = F]
        bulk_cell_metadata <- as.data.frame(bulk_cell_metadata %>% mutate(Replicate = data.table::rleid(
                                                       binnedTime, Path1, Path2)))
        bulk_cell_metadata <- bulk_cell_metadata[, c("binnedTime", "Replicate", "Path1", "Path2"),drop = F]
        colnames(bulk_cell_metadata) <- c("Pooled_Pseudotime", "Replicate", "Path1", "Path2")
        
        if(all(c(all(colnames(bulk_counts_list[["pseudo_bulk_counts_sum"]]) == rownames(bulk_cell_metadata)),
        all(colnames(bulk_counts_list[["pseudo_bulk_counts_sum"]]) == rownames(bulk_cell_metadata)),
        all(colnames(bulk_counts_list[["pseudo_bulk_counts_sum"]]) == rownames(bulk_cell_metadata)))) ==F){
            stop(paste("Mismatch for", j))
        }
        
        data.list <- list(count_table = list(sum = bulk_counts_list[["pseudo_bulk_counts_sum"]],
                                             mean = bulk_counts_list[["pseudo_bulk_counts_mean"]], 
                                             median = bulk_counts_list[["pseudo_bulk_counts_median"]]),
                          cell_meta = bulk_cell_metadata)
        
        fName <- paste0("benchMarking_Normalization/data/input/PB_Input/",
                        str_remove(j, ".tsv"), ".RDS")
        saveRDS(data.list,fName)
        
        x <- gc()
        x <- NULL
    }
    x <- gc()
    x <- NULL
}

cat(paste0("\n","Script Finished"))
