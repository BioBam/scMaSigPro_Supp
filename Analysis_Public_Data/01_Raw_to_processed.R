#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Pre-Process and CC Remove ##
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(parallel))

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Read BioMart info
biomart.anno <- readRDS(paste(dirPath, "cell_cycle_data.mart", sep = "/"))
reg.out <- c(unique(biomart.anno$SYMBOL))

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
names(rep_vec) <- rep_vec

# Run lapply
umaps.list <- lapply(rep_vec, function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
    # Step-1: Add Annotation for donors
    if (rep_i == "rep1") {
        individual <- "Donor-1"
        age <- "35"
        sex <- "Male"
    } else if (rep_i == "rep2") {
        individual <- "Donor-2"
        age <- "28"
        sex <- "Female"
    } else if (rep_i == "rep3") {
        individual <- "Donor-3"
        age <- "19"
        sex <- "Female"
    }
    
    # Step-2: Load the filtered Matrix
    filt_mat <- Read10X_h5(
        paste(inPath, rep_i, "filtered_feature_bc_matrix.h5", sep = "/")
        )
    # Renaming the features
    filt_mat@Dimnames[[1]] <- str_replace(filt_mat@Dimnames[[1]], "_", "-")
    
    # Create Raw Seurat Object
    sob.raw <- CreateSeuratObject(
        counts = filt_mat, min.cells = 100,
        min.features = 100, project = rep_i
    )
    
    # Calculate Percentage of Mitochondrial Reads
    sob.raw[["percent.mt"]] <- PercentageFeatureSet(sob.raw, pattern = "^MT-")
    
    # Remove Cells 
    sob.sub <- subset(sob.raw, subset = nFeature_RNA > 100 & nCount_RNA < 40000 & percent.mt < 10)
    
    # Normalize
    sob.prs <- NormalizeData(sob.sub, verbose = F)
    
    # Find Variable features
    sob.prs <- FindVariableFeatures(sob.prs, selection.method = "dispersion", verbose = F)
    
    # Check how many genes are present in the dataset
    indata <- rownames(sob.prs)[rownames(sob.prs) %in% reg.out]
    
    # Keep only those which are present in data
    biomart.anno <- biomart.anno[biomart.anno$SYMBOL %in% indata, ]
    
    # For seurat we need to divide genes into vectors
    s.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0006260"), "SYMBOL"])
    g2m.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0000087", "GO:0000279", "GO:0007059", "GO:0048285"), "SYMBOL"])
    
    # Cell Cycle Scores
    sob.prs <- CellCycleScoring(sob.prs,
                               s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = TRUE
    )
    
    # Regress Cell Cycle Score
    sob.prs <- ScaleData(sob.prs,
                        vars.to.regress = c("S.Score", "G2M.Score"),
                        features = rownames(sob.prs), verbose = T,
    )
    
    # Save 
    file_name <- paste(dirPath, rep_i, paste(rep_i, "prs", sep = "_"), sep = "/")
    SaveH5Seurat(
        object = sob.prs, filename = file_name,
        overwrite = T, verbose = FALSE
    )
    
    # Return UMAP
    return(NULL)
    
})