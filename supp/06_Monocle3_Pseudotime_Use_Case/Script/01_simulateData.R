# Title: Analyze 60% Zero-Inflated data with Monocle3
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(monocle3))

# Set Paths relative to project
inPath <- "supp/06_Monocle3_Pseudotime_Use_Case/data/input/sceObjects/"
outPath <- "supp/06_Monocle3_Pseudotime_Use_Case/data/output/cell_data_set/"
dir.create(outPath, showWarnings = F, recursive = T)

# Load cds
sce <- readRDS(file = paste0(inPath, "zi_mid_60_mid_0_shape_0.25.RDS"))

# Convert to cell dataset ojet
cds <- new_cell_data_set(sce@assays@data@listData$counts,
                         cell_metadata = colData(sce),
                         gene_metadata = rowData(sce))

# Pre-process
cds <- preprocess_cds(cds, num_dim = 10, method = "PCA", norm_method = "log")

# Reduce Dimensions
cds <- reduce_dimension(cds, max_components = 2, reduction_method = "UMAP",
                        preprocess_method = "PCA")

# Plot the structure with simulated pseudotime
plot_cells(cds, show_trajectory_graph = F, color_cells_by = "Step",
           cell_size = 2)

# Run Clustering
cds <- cluster_cells(cds, resolution = 0.001)

# Plot with inferred clusters and partitions
plot_cells(cds, show_trajectory_graph = F, color_cells_by = "cluster",
           cell_size = 2)
plot_cells(cds, show_trajectory_graph = F, color_cells_by = "partition",
           cell_size = 2)

# Learn Graph
cds <- learn_graph(cds)
plot_cells(cds, show_trajectory_graph = T, color_cells_by = "cluster",
           cell_size = 2)

# Order cells
cds <- order_cells(cds, root_cells = rownames(cds@colData[cds@colData$Step %in% c(1:10),]))

# View
plot_cells(cds, show_trajectory_graph = T, color_cells_by = "pseudotime",
           cell_size = 2)

# # Detect Differential genes
# cds_pgraph <- graph_test(cds,
#                          neighbor_graph="principal_graph", 
#                          expression_family = "negbinomial"
#                          )

cds_reg <- fit_models(cds,
                         model_formula_str = "~Step",cores = 1, 
                         expression_family = "negbinomial"
)



