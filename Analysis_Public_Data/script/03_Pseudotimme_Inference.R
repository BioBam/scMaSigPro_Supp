##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: Monocle3 Analysis #############################
##########################################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(parallel))

# Load helper script
source("R_Scripts/helper_function/detect_root_nodes().R")

# Prefix
prefixIn <- "Analysis_Public_Data/data/"
prefixOut <- "Analysis_Public_Data/data/"

# Get folder names
rep_vec <- list.dirs(prefixIn, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "Setty_et_al_2019_Integrated_sob.h5seurat", "Human_Cell_Atlas"))]
names(rep_vec) <- rep_vec

# Create the CDS per donor
cds.list <- mclapply(rep_vec, FUN = function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
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
    
    # Create CDS
    cds <- new_cell_data_set(
        expression_data = sob@assays$RNA@counts,
        cell_metadata = DataFrame(sob@meta.data),
        gene_metadata = data.frame(
            gene_short_name = rownames(sob@assays$RNA@counts),
            row.names = rownames(sob@assays$RNA@counts)
        )
    )
    
    # Return CDS
    return(cds)
    
}, mc.cores = detectCores(), mc.set.seed = 123)

# Basic Pre-Processing and log Normalization
cds.list <- lapply(cds.list, FUN = function(cds, inPath = prefixIn, outPath = prefixOut){
    
    # Pre-process
    cds <- preprocess_cds(cds)
    
    # Reduce dimensions
    cds <- reduce_dimension(cds)
    
    # Cluster the Data
    cds <- cluster_cells(cds)
    
    # Learn Graph
    cds <- learn_graph(cds, use_partition = T,
                       close_loop = F,
                       list(ncenter=50, prune_graph =F))
    
    # Calculate the Pseudotime
    cds <- order_cells(cds,
                       reduction_method = "UMAP",
                       root_pr_nodes = find_root_pp(cds))
    
    # Return
    return(cds)
    
#}, mc.cores = detectCores(), mc.set.seed = 123)
})

# Basic Pre-Processing and log Normalization
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
    combined_plot <- ggarrange( = rep_plots, ncol = 3, nrow = 3)
    
    print(combined_plot)
    
    # Save
    ggsave(file = paste0(paste(inPath, rep_i ,sep = "/"), "/",rep_i, "_Markers.png"),
           plot = combined_plot, width = 16, height = 16, dpi = 1600, limitsize = FALSE)
    
    return(rep_plots)
})




# HPO in Monocyte Like cells
ggsave(p,
  filename = paste0(prefixOut, i, "_MPO_Trend.png"),
  dpi = 600, limitsize = FALSE
)
