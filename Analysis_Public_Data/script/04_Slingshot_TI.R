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
    sob <- LoadH5Seurat(
        paste0(paste(inPath, rep_i,sep = "/"), "/",rep_i, "_azimuth.h5seurat"),
        verbose = FALSE
    )
    
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
    sce <- runTSNE(sce, exprs_values = "norm")
    
    # Clustring
    sce@colData$mclusters <- as.factor(as.character(Mclust(reducedDim(sce))$classification))
    
    # Return CDS
    return(sce)
    
}, mc.cores = detectCores(), mc.set.seed = 123)


# Basic Pre-Processing and log Normalization
sce.plot.list <- lapply(sce.list, FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
    # Plot UMAP and PCA
    #pca <- plotPCA(rep_i, colour_by = "predicted.celltype.l2")
    #umap <- plotUMAP(rep_i, colour_by = "predicted.celltype.l2")
    tsne_cell = plotTSNE(rep_i, colour_by = "predicted.celltype.l2")
    tsne_cluster = plotTSNE(rep_i, colour_by = "mclusters")
    tsne <- tsne_cell + tsne_cluster
    
    # Return
    return(list(
        TSNE = tsne
        ))
    
    #}, mc.cores = detectCores(), mc.set.seed = 123)
})

# Run slingshot
sling.sce.list <-  mclapply(sce.list, FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    rep_i <- slingshot(rep_i, clusterLabels = 'mclusters', reducedDim = 'TSNE',
                       start.clus = find_root_pp(rep_i,
                                                 objType = "sling")$root,
                       end.clus = find_root_pp(rep_i,
                                               objType = "sling")$end
                       )
}, mc.cores = detectCores(), mc.set.seed = 123)


# Plot trajectory
sling.plot.list <-  mclapply(sling.sce.list, FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
    df <- as.data.frame(reducedDim(rep_i))
    df <- cbind(df, data.frame(mclusters = rep_i@colData$mclusters,
                               azimuth_cells = rep_i@colData$predicted.celltype.l2))
    curve.df <- slingCurves(rep_i, as.df = TRUE)
    
    ggplot(df, aes(x = TSNE1, y = TSNE2)) +
        geom_point(aes(fill = azimuth_cells), col = "grey70", shape = 21) + 
        geom_path(data = curve.df, aes(x = TSNE1, y = TSNE2))+
        theme_classic()
    
}, mc.cores = detectCores(), mc.set.seed = 123)



df <- as.data.frame(reducedDim(sling.sce.list$rep1))
df <- cbind(df, data.frame(mclusters = sling.sce.list$rep1@colData$mclusters,
                           azimuth_cells = sling.sce.list$rep1@colData$predicted.celltype.l2))
curve.df <- slingCurves(sling.sce.list$rep1, as.df = TRUE)


ggplot(df, aes(x = TSNE1, y = TSNE2)) +
    geom_point(aes(fill = mclusters), col = "grey70", shape = 21) + 
    geom_path(data = curve.df, aes(x = TSNE1, y = TSNE2))+
    theme_classic()

ggplot(df, aes(x = TSNE1, y = TSNE2)) +
    geom_point(aes(fill = azimuth_cells), col = "grey70", shape = 21) + 
    geom_path(data = curve.df, aes(x = TSNE1, y = TSNE2))+
    theme_classic()


Slomgh# Basic Pre-Processing and log Normalization
cds.umap.list <- lapply(cds.list, FUN = function(cds, inPath = prefixIn, outPath = prefixOut){
    
    # Plot the Graph
    pseudotime <- plot_cells(cds, color_cells_by = "pseudotime", alpha = 0.8,
                             cell_stroke = 0.5,  trajectory_graph_segment_size = 1.5,
                             label_cell_groups = F, label_branch_points = F, label_leaves = F,
                             label_principal_points = F, label_roots = T) +
        ggtitle("Inferred Pseudotime") + theme(legend.position = "bottom")
    
    # Cell Type
    cell_type <- plot_cells(cds, color_cells_by = "predicted.celltype.l2", alpha = 0.8,
                            cell_stroke = 0.5,  trajectory_graph_segment_size = 1.5,
                            label_cell_groups = F, label_branch_points = T, label_leaves = F,
                            label_principal_points = F, label_roots = T) +
        ggtitle("Cell Type") + theme(legend.position = "bottom")
    
    # Rep
    rep_i <- levels(unique(cds@colData@listData$orig.ident))
    
    # Save
    save(cds, file = paste0(paste(inPath, rep_i ,sep = "/"), "/",rep_i, "_processed.RData"))
    
    # Join plots
    umap <- ggarrange(pseudotime, cell_type)
    
    # Save Images
    ggsave(umap,
           filename = paste0(paste(inPath, rep_i ,sep = "/"), "/",rep_i, "_Monocle3.png"),
           dpi = 800, limitsize = FALSE, width = 11, height = 6
    )
    
    return(umap)
})

#---------------------------------------

# Plot Markers

marker_list <- c("CD34", "MPO", "GATA2", "GATA1", "IRF8", "EPOR", "APOBEC3C", "ACY3", "RUNX1")

# Marker list
marker.list <- lapply(cds.list, FUN = function(cds, inPath = prefixIn, outPath = prefixOut, markerList = marker_list){
    
    rep_plots <- lapply(markerList, function(gene_i){
        gene_i_plot <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == gene_i, ],
                                                color_cells_by = "predicted.celltype.l2")
        return(gene_i_plot)
    })
    names(rep_plots) <- markerList
    
    # Rep
    rep_i <- levels(unique(cds@colData@listData$orig.ident))
    
    # Combined Plot
    combined_plot <- ggarrange(rep_plots[[1]], rep_plots[[2]], rep_plots[[3]], 
                               rep_plots[[4]], rep_plots[[5]], rep_plots[[6]],
                               rep_plots[[7]], rep_plots[[8]], rep_plots[[9]],
                               ncol = 3, nrow = 3)
    
    # Save
    ggsave(file = paste0(paste(inPath, rep_i ,sep = "/"), "/",rep_i, "_Markers.png"),
           plot = combined_plot, width = 12, height = 8, dpi = 800, limitsize = FALSE)
    
    return(rep_plots)
})


# HPO in Monocyte Like cells
ggsave(p,
       filename = paste0(prefixOut, i, "_MPO_Trend.png"),
       dpi = 600, limitsize = FALSE
)
