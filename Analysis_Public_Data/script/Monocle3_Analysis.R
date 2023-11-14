# Load libraries
library(tidyverse)
library(monocle3)

# Load helper script
source("R_Scripts/helper_function/detect_root_nodes().R")

# Load data for donor-3
don1Counts <- readRDS("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/don1.counts.RDS")
don1CellMeta <- readRDS("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/don1.meta.RDS")

# Create CDS
cds <- new_cell_data_set(
    expression_data = don1Counts,
    cell_metadata = don1CellMeta
    #gene_metadata = geneMeta
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

# Explore Data
cell_type <- plot_cells(cds, color_cells_by = "cell_type",cell_size = 2) +
    theme(legend.position = "bottom")

# Compute Clusters
cds <- cluster_cells(
    cds = cds,
    resolution = 0.01)

# Visulizse Clusters
clusters <- plot_cells(cds, color_cells_by = "cluster",cell_size = 2) +
    theme(legend.position = "bottom")

cell_type + clusters

# Drop cluster 31, 14, 33, (Donor1)
cds <- cds[,!(clusters(cds) %in% c(20, 8))]

# Drop genes
keep <- (rownames(
    cds@assays@data@listData$counts)[rowSums(
        cds@assays@data@listData$counts
    ) >= 10]
)
cds <- cds[keep,]

# Compute UMAP and PCA again
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
    umap.min_dist = 0.05, # Increase this to relax global structure
    reduction_method = "UMAP",
    preprocess_method = "PCA"
)

# Explore Data
cell_type <- plot_cells(cds, color_cells_by = "cell_type",cell_size = 2) +
    theme(legend.position = "bottom")

# Compute Clusters
cds <- cluster_cells(
    cds = cds)

# Visulizse Clusters
clusters <- plot_cells(cds, color_cells_by = "cluster",cell_size = 2) +
    theme(legend.position = "bottom")

cell_type + clusters

# Learn graph
cds <- learn_graph(cds, use_partition = T,
                   close_loop = F,
                   list(prune_graph =F))

# Calculate the Pseudotime
cds <- order_cells(cds,
                   reduction_method = "UMAP",
                   root_pr_nodes = 
                       find_root_pp(cds, cell_col = "cell_type", cell = "HMP"))

pTime <-  plot_cells(cds, color_cells_by = "pseudotime",cell_size = 2) +
    theme(legend.position = "bottom")

cell_type + clusters + pTime


# Save For ScMaSigPro anaLysis
saveRDS(cds, "Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/Monocle3_Processed_Donor1.RDS")
