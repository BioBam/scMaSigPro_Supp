# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtools))

pltColSums <- function(countTable){
    countTable <- as.matrix(countTable)
    total_counts_per_cell <- colSums(countTable)
    df <- data.frame(total_counts = total_counts_per_cell)
    p <- ggplot(df, aes(x = total_counts)) +
        geom_histogram(bins = 30, fill = "lightblue", color = "black") +
        labs(title = "Histogram of total counts per cell",
             x = "Total counts", 
             y = "Number of cells") +
        theme_minimal() + theme(axis.text = element_text(size = rel(1.5)), 
                                plot.title = element_text(size = rel(3)),
                                axis.title = element_text(size = rel(2)))
    return(p)
}


# Load Prefix
parentFix <- "benchMarking_Normalization/" 
scriptFix <- paste0("old_Scripts/")
inPrefix <- paste0(parentFix, "data/output/RawCounts/")
outPrefix <- paste0(parentFix, "data/output/")

# Load Custim Functions
load <- function(){
    source(paste0(scriptFix, "plot_simulations().R"))
    source(paste0(scriptFix, "add_gene_anno().R"))
    source(paste0(scriptFix, "calc_bin_size.R"))
    source(paste0(scriptFix, "calcNormCounts.R"))
}

# Load File names
all_file_names <- list.files(inPrefix)
all_file_names <- str_remove(all_file_names, "_RawCounts.tsv")
split_list <- str_split(all_file_names, "_")
names(split_list) <- list.files(inPrefix)

load()

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
    rawCount <- as.matrix(read.table(
        file = paste0(inPrefix, file_name), 
        header = T, sep = "\t", row.names = 1))
    
    p <- pltColSums(rawCount)
    fname <- paste0(paste0(inPrefix,str_remove(i, ".tsv"), ".png"))
    ggsave(p, dpi = 1200, filename = fname, width = 9, height = 7)
    
    # Log-Normalize
    dir_path <- paste0(outPrefix, "LogNorm/")
    dir.create(dir_path, showWarnings = F)
    normCountsName <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "LogNorm.tsv")
    normCounts <- as.matrix(calcNormCounts(rawCounts = rawCount, cat = "logLibSize", size_fac = 10000))
    write.table(normCounts, normCountsName, row.names = T, sep = "\t", col.names = NA)
    p <- pltColSums(normCounts)
    fname <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "LogNorm.png")
    ggsave(p, dpi = 1200, filename = fname, width = 9, height = 7)
    
    # CLR
    dir_path <- paste0(outPrefix, "FQNorm/")
    dir.create(dir_path, showWarnings = F)
    normCountsName <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "FQNorm.tsv")
    normCounts <- as.matrix(calcNormCounts(rawCounts = rawCount, cat = "FQNorm", size_fac = 10000))
    write.table(normCounts, normCountsName, row.names = T, sep = "\t", col.names = NA)
    p <- pltColSums(normCounts)
    fname <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "FQNorm.png")
    ggsave(p, dpi = 1200, filename = fname, width = 9, height = 7)
    
    # Relative Counts
    dir_path <- paste0(outPrefix, "RcNorm/")
    dir.create(dir_path, showWarnings = F)
    normCountsName <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "RcNorm.tsv")
    normCounts <- as.matrix(calcNormCounts(rawCounts = rawCount, cat = "libSize", size_fac = 10000))
    write.table(normCounts, normCountsName, row.names = T, sep = "\t", col.names = NA)
    p <- pltColSums(normCounts)
    fname <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "RcNorm.png")
    ggsave(p, dpi = 1200, filename = fname, width = 9, height = 7)
    
    # SCT Transform
    dir_path <- paste0(outPrefix, "SctNorm/")
    dir.create(dir_path, showWarnings = F)
    normCountsName <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "SctNorm.tsv")
    normCounts <- as.matrix(calcNormCounts(rawCounts = rawCount, cat = "sctransform", size_fac = 10000))
    write.table(normCounts, normCountsName, row.names = T, sep = "\t", col.names = NA)
    p <- pltColSums(normCounts)
    fname <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "SctNorm.png")
    ggsave(p, dpi = 1200, filename = fname, width = 9, height = 7)
    
    # CPM
    dir_path <- paste0(outPrefix, "CLRNorm/")
    dir.create(dir_path, showWarnings = F)
    normCountsName <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "CLRNorm.tsv")
    normCounts <- as.matrix(calcNormCounts(rawCounts = rawCount, cat = "CLR", size_fac = 10000))
    write.table(normCounts, normCountsName, row.names = T, sep = "\t", col.names = NA)
    p <- pltColSums(normCounts)
    fname <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "CLRNorm.png")
    ggsave(p, dpi = 1200, filename = fname, width = 9, height = 7)
    
    # Shifted Logs
    dir_path <- paste0(outPrefix, "CPMNorm/")
    dir.create(dir_path, showWarnings = F)
    normCountsName <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "CPMNorm.tsv")
    normCounts <- as.matrix(calcNormCounts(rawCounts = rawCount, cat = "libSize", size_fac = 1000000))
    write.table(normCounts, normCountsName, row.names = T, sep = "\t", col.names = NA)
    p <- pltColSums(normCounts)
    fname <- paste0(dir_path, test_name, "_", test_mid, "_",test_value, "_", "CPMNorm.png")
    ggsave(p, dpi = 1200, filename = fname, width = 9, height = 7)
    
}  

cat(paste0("\n","Script Finished"))