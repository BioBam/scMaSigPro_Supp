##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: Monocle3 Analysis #############################
##########################################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(mclust))

# Load helper script
source("R_Scripts/helper_function/detect_root_nodes().R")
source("R_Scripts/helper_function/FQnorm.R")

# Prefix
prefixIn <- "Analysis_Public_Data/data/"
prefixOut <- "Analysis_Public_Data/data/"

# Get folder names
rep_vec <- list.dirs(prefixIn, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "Setty_et_al_2019_Integrated_sob.h5seurat", "Human_Cell_Atlas"))]
names(rep_vec) <- rep_vec

# Create the CDS per donor
sce.list <- mclapply(rep_vec, FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
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
    
    # Step-2: Load the Seurat Object
    sob.prs <- LoadH5Seurat(
        paste0(paste(inPath, rep_i,sep = "/"), "/",rep_i, "_azimuth.h5seurat"),
        verbose = FALSE
    )
    
    # Run PCA
    sob.prs <- RunPCA(sob.prs,
                      features = VariableFeatures(sob.prs),
                      verbose = F, ndims.print = 0, nfeatures.print = 0
    )
    
    # Compute UMAP
    sob.prs <- FindNeighbors(sob.prs, verbose = F)
    sob.prs <- FindClusters(sob.prs, verbose = F)
    sob.prs <- RunUMAP(sob.prs, dims = 1:50, verbose = F)
    sob <- sob.prs
    
    # Create SCE
    sce <- SingleCellExperiment(assays = List(counts = sob@assays$RNA@counts, # raw Counts
                                              norm = FQnorm(sob@assays$RNA@counts)),
                                colData = DataFrame(sob@meta.data),
                                rowData = DataFrame(
                                    data.frame(
                                        gene_short_name = rownames(sob@assays$RNA@counts),
                                        row.names = rownames(sob@assays$RNA@counts)
                                        )))
    
    # Dimension reduction
    reducedDims(sce)$UMAP <- sob@reductions$umap@cell.embeddings
    
    # Clustring
    sce@colData$seurat_clusters <- as.factor(as.character(sob@meta.data$seurat_clusters))
    
    # Return CDS
    return(sce)
    
}, mc.cores = detectCores(), mc.set.seed = 123)


# Basic Pre-Processing and log Normalization
sce.plot.list <- mclapply(sce.list, FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
    # Plot UMAP and PCA
    umap_cell = plotUMAP(rep_i, colour_by = "predicted.celltype.l2")
    umap_cluster = plotUMAP(rep_i, colour_by = "seurat_clusters")
    umap <- umap_cell + umap_cluster
    
    # Return
    return(list(
        UMAP = umap
        ))
}, mc.cores = detectCores(), mc.set.seed = 123)

# Run slingshot
sling.sce.list <-  mclapply(sce.list, FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    rep_i <- slingshot(rep_i, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP',
                       start.clus =  find_root_pp(rep_i,factor_col = "seurat_clusters",
                                                  cell_col = "predicted.celltype.l2",
                                                  objType = "sling")$root,
                       end.clus = find_root_pp(rep_i,factor_col = "seurat_clusters",
                                               cell_col = "predicted.celltype.l2",
                                               objType = "sling")$end
                       )
}, mc.cores = detectCores(), mc.set.seed = 123)

# Plot trajectory
sling.plot.list <-  mclapply(names(sling.sce.list), FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
    slingsce <- sling.sce.list[[rep_i]]
    
    df <- as.data.frame(reducedDim(slingsce))
    df <- cbind(df, data.frame(seurat_clusters = slingsce@colData$seurat_clusters,
                               azimuth_cells = slingsce@colData$predicted.celltype.l2))
    curve.df <- slingCurves(slingsce, as.df = TRUE)
    
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(fill = azimuth_cells), col = "grey70", shape = 21) + 
        geom_path(data = curve.df, aes(x = UMAP_1, y = UMAP_2))+
        theme_classic()
    
    ggsave(p,
           filename = paste0(paste(inPath, rep_i ,sep = "/"), "/",rep_i, "_Slingshot.png"),
           dpi = 800, limitsize = FALSE, width = 11, height = 6
    )
    
    saveRDS(slingsce,
            file = paste0(paste(inPath, rep_i ,sep = "/"), "/",rep_i, "_Slingshot.RDS")
            )
}, mc.cores = detectCores(), mc.set.seed = 123)
