# Title: Analyze 60% Zero-Inflated data with Monocle3
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(igraph))

# Set Paths relative to project
inPath <- "supp/06_Monocle3_Pseudotime_Use_Case/data/input/sceObjects/"
outPath <- "supp/06_Monocle3_Pseudotime_Use_Case/data/output/cell_data_set/"
dir.create(outPath, showWarnings = F, recursive = T)

# Load cds
sce <- readRDS(file = paste0(inPath, "zi_mid_40_mid_0_shape_-0.15.RDS"))

# Convert to cell dataset ojet
cds <- new_cell_data_set(sce@assays@data@listData$counts,
  cell_metadata = colData(sce),
  gene_metadata = rowData(sce)
)

# Pre-process
cds <- preprocess_cds(cds, num_dim = 5, method = "PCA", norm_method = "log")

# Reduce Dimensions
cds <- reduce_dimension(cds,
  max_components = 2, reduction_method = "UMAP",
  preprocess_method = "PCA"
)

# Plot the structure with simulated pseudotime
plot_cells(cds,
  show_trajectory_graph = F, color_cells_by = "Batch",
  cell_size = 1.5, label_cell_groups = F
)

# Run Clustering
cds <- cluster_cells(cds, resolution = 0.001)

# Plot with inferred clusters and partitions
plot_cells(cds,
  show_trajectory_graph = F, color_cells_by = "cluster",
  cell_size = 2
)
plot_cells(cds,
  show_trajectory_graph = F, color_cells_by = "partition",
  cell_size = 2
)

# Learn Graph
cds <- learn_graph(cds)
plot_cells(cds,
  show_trajectory_graph = T, color_cells_by = "cluster",
  cell_size = 2, trajectory_graph_segment_size = 3
)

# Order cells
cds <- order_cells(cds, root_cells = rownames(cds@colData[cds@colData$Step %in% c(1:10), ]))

# View
plot_cells(cds,
  show_trajectory_graph = T, color_cells_by = "pseudotime",
  cell_size = 2, label_roots = T, label_principal_points = T,
  trajectory_graph_segment_size = 3, trajectory_graph_color = "black"
)
# View
plot_cells(cds,
  show_trajectory_graph = T, color_cells_by = "partition",
  cell_size = 2, label_roots = T, label_principal_points = T
)

# Extract Prinicipal graph
y_to_cells <- principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- paste0("Y_", y_to_cells$V1)
y_to_cells$V1 <- NULL

# Get the root vertices
# It is the same node as above
root <- cds@principal_graph_aux$UMAP$root_pr_nodes
mst <- principal_graph(cds)$UMAP

# Get the end points
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

cellWeights <- lapply(endpoints, function(endpoint) {
  # We find the path between the endpoint and the root
  path <- paste("Y", as.character(
    shortest_paths(mst, root, endpoint)$vpath[[1]]
  ),
  sep = "_"
  )

  # Get the cells along the path
  df <- y_to_cells[y_to_cells$Y %in% path, ]

  # Binarize per path
  df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>%
  do.call(what = "cbind", args = .) %>%
  as.matrix()
rownames(cellWeights) <- colnames(cds)

# Get the Pseudotime per path
pseudotime <- matrix(pseudotime(cds),
  ncol = ncol(cellWeights),
  nrow = ncol(cds), byrow = FALSE
)

# Check for universal pseudotime
# Assuming 'df' is your dataframe
result <- all(apply(pseudotime, 1, function(x) length(unique(x)) == 1))

# Universal Pseudotime exist
path_lengths <- apply(cellWeights, 2, FUN = function(x) {
  length(x[x == 1])
})

cellWeights <- as.data.frame(cellWeights)
# Plot Each length
cell.meta <- as.data.frame(colData(cds))
cell.meta[cell.meta$Cell %in% rownames(cellWeights[cellWeights$Y_5 == 1, ]), "Lineage"] <- "Y_5"
cell.meta[cell.meta$Cell %in% rownames(cellWeights[cellWeights$Y_53 == 1, ]), "Lineage"] <- "Y_53"
cell.meta[cell.meta$Cell %in% rownames(cellWeights[cellWeights$Y_68 == 1, ]), "Lineage"] <- "Y_68"
cell.meta[cell.meta$Cell %in% rownames(cellWeights[cellWeights$Y_78 == 1, ]), "Lineage"] <- "Y_78"
cell.meta[cell.meta$Cell %in% rownames(cellWeights[cellWeights$Y_83 == 1, ]), "Lineage"] <- "Y_83"
cell.meta[cell.meta$Cell %in% rownames(cellWeights[cellWeights$Y_92 == 1, ]), "Lineage"] <- "Y_92"

cell.meta$Pseudotime <- pseudotime[, 1]
colData(cds) <- DataFrame(cell.meta)

# Major lineage Detected are 92 and 78
# Subset
cds_lin <- cds[, colnames(cds) %in% cell.meta[cell.meta$Lineage %in% c("Y_92", "Y_78"), "Cell"]]
plot_cells(cds_lin,
  show_trajectory_graph = T, color_cells_by = "Lineage",
  cell_size = 2, label_roots = T, label_principal_points = T
)
plot_cells(cds_lin,
  show_trajectory_graph = T, color_cells_by = "pseudotime",
  cell_size = 2, label_roots = T, label_principal_points = T
)

# get the cell metadata
cell.meta <- as.data.frame(colData(cds_lin))

# Subset the cell Metadata
cell.meta <- cell.meta[, c("Lineage", "Pseudotime")]

# Separate bu lineage
cell.meta$Lin1 <- ifelse(cell.meta$Lineage == "Y_92", 1, 0)
cell.meta$Lin2 <- ifelse(cell.meta$Lineage == "Y_78", 1, 0)
# Remove cells that does not belong to any lineage
cell.meta <- cell.meta[!(cell.meta$Lin1 == 0 & cell.meta$Lin2 == 0), ]

# Make Pseudobulk expression
source("old_Scripts/entropy_discretize.R")
pbProfile <- entropy_discretize(cell.meta, time_col = "Pseudotime", drop.fac = 0.7)

# Create Bin counts
source("old_Scripts/calc_bin_size.R")
bulk_cell_metadata <- make_pseudobulk_design(
  design.file = pbProfile, addBinSize = T,
  paths.vector = c("Lin1", "Lin2"),
  binnedCol = "binnedTime"
)

# Derive counts
bulk_counts_list <- make_bulk_counts(
  counts = as.matrix(cds_lin@assays@data@listData$counts),
  cluster_member_col = "cluster.members",
  bin_col = "bin", pseudo_bulk_profile = bulk_cell_metadata,
  round = T
)

# Create New SCE Object
meanCounts <- as.matrix(bulk_counts_list$pseudo_bulk_counts_mean)
sumCounts <- as.matrix(bulk_counts_list$pseudo_bulk_counts_sum)
medianCounts <- as.matrix(bulk_counts_list$pseudo_bulk_counts_median)

# Import scMaSigPro
library(maSigPro)

# MasigPro Daat
bulk_cell_metadata
newDesign <- bulk_cell_metadata[, c(1, 5, 6, 7)]
design <- make.design.matrix(newDesign,
  degree = 2, time.col = 1,
  repl.col = 2
)

library(MASS)
p.vec <- p.vector(sumCounts,
  design = design,
  counts = T
)
t.fit <- T.fit(p.vec,
  design = design$dis,
  family = p.vec$family
)

sigs <- get.siggenes(t.fit, vars = "groups", rsq = 0.6)
suma2Venn(sigs$summary[, c(1:2)])

PlotGroups(
  data = sumCounts[91, ], edesign = design$edesign, show.fit = T,
  dis = design$dis, groups.vector = design$groups.vector
)
