# Load libraries
library(tidyverse)
library(monocle3)

# Load helper script
source("R_Scripts/helper_function/detect_root_nodes().R")

# Load And Create Monocle 3 objects

lapply(c("don1", "don2", "don3"), function(don){
    
    # Load
    don_Counts <- readRDS(paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/", don, ".counts.RDS"))
    don_CellMeta <- readRDS(paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/", don, ".meta.RDS"))
    
    # Create CDS
    cds <- new_cell_data_set(
        expression_data = don_Counts,
        cell_metadata = don_CellMeta,
        gene_metadata = data.frame(
            row.names = rownames(don_Counts),
            gene_short_name = rownames(don_Counts)
        )
    )
    
    # Pre-process
    cds <- preprocess_cds(cds,
                          method = "PCA",
                          num_dim = 25,
                          norm_method = "log",
                          scaling = T)
    
    # Compute UMAP
    cds <- reduce_dimension(
        cds = cds,
        max_components = 2,
        umap.min_dist = 0.1, # Increase this to relax global structure
        reduction_method = "UMAP",
        preprocess_method = "PCA"
    )
    
    # Compute Clusters
    cds <- cluster_cells(
        cds = cds,
        resolution = 0.01)
    
    # Visulizse Clusters
    clusters <- plot_cells(cds, color_cells_by = "cluster",cell_size = 2) +
        theme(legend.position = "bottom")
    
    if(don == "don1"){
        cds <- cds[,!(clusters(cds) %in% c(20, 8))]
    }
    
    # Drop genes
    keep <- (rownames(
        cds@assays@data@listData$counts)[rowSums(
            cds@assays@data@listData$counts
        ) >= 100]
    )
    cds <- cds[keep,]
    
    # Pre-process
    cds <- preprocess_cds(cds,
                          method = "PCA",
                          num_dim = 10,
                          norm_method = "log",
                          scaling = T)
    
    # Compute UMAP
    cds <- reduce_dimension(
        cds = cds,
        max_components = 2,
        umap.min_dist = 0.1, # Increase this to relax global structure
        reduction_method = "UMAP",
        preprocess_method = "PCA"
    )
    
    # Compute Clusters
    cds <- cluster_cells(
        cds = cds)
    
    # Learn graph
    cds <- learn_graph(cds, use_partition = T,
                       close_loop = F,
                       list(prune_graph =F))
    
    # Calculate the Pseudotime
    cds <- order_cells(cds,
                       reduction_method = "UMAP",
                       root_pr_nodes = 
                           find_root_pp(cds, cell_col = "cell_type", cell = "HMP"))
    
    # Save For ScMaSigPro anaLysis
    saveRDS(cds, 
            paste0(
                "Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/Monocle3_Processed_",don,".RDS")
            )
})
