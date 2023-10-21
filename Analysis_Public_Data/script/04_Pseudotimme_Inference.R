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

# Prefix
prefixIn <- "Analysis_Public_Data/data/"
prefixOut <- "Analysis_Public_Data/data/"

#  the Processed Seurat Object
sob <- LoadH5Seurat(file = paste(prefixIn, "Setty_et_al_2019_Integrated_sob.h5seurat", sep = "/"), verbose = F)

# Construct Cell Dataset object
cds <- new_cell_data_set(
  expression_data = sob@assays$integrated@scale.data,
  cell_metadata = DataFrame(sob@meta.data),
  gene_metadata = data.frame(
    gene_short_name = rownames(sob@assays$integrated@scale.data),
    row.names = rownames(sob@assays$integrated@scale.data)
  )
)


# Pre-process
cds <- preprocess_cds(cds, num_dim = 50, scaling = F, norm_method = "none")

# Reduce dimensions
cds <- reduce_dimension(cds,
  max_components = 3,
  preprocess_method = "PCA"
)

# Cluster the Data
cds <- cluster_cells(cds)

# Learn Graph
cds <- learn_graph(cds)

# Cell metadata
cell_metadata <- as.data.frame(colData(cds))

# Select barcodes
root_barcodes <- rownames(cell_metadata[cell_metadata$azimuth_celltype == "HSC", ])

# Calculate the Pseudotime
cds <- order_cells(cds,
  reduction_method = "UMAP",
  root_cells = root_barcodes
)

# Plot the Graph
plot_cells(cds, color_cells_by = "pseudotime", alpha = 0.8,
                              cell_stroke = 0.5,  trajectory_graph_segment_size = 1.5,
                label_cell_groups = F, label_branch_points = F, label_leaves = F,
                label_principal_points = F, label_roots = T) +
    ggtitle("Inferred Pseudotime") + theme(legend.position = "bottom")

f <- plot_cells(cds, color_cells_by = "pseudotime", alpha = 0.8,
                cell_stroke = 0.5,  trajectory_graph_segment_size = 1.5,
                label_cell_groups = F, label_branch_points = F, label_leaves = F,
                label_principal_points = T, label_roots = F) +
    ggtitle("Inferred Pseudotime") + theme(legend.position = "bottom")

combined.map <- ggarrange(
    a,b,c,d,e,f, ncol = 3, nrow = 2
)

combined.map.2 <- ggarrange(
   d, e,f, ncol = 3, nrow = 1
)

ggsave(combined.map.2,
  filename = paste0(prefixOut, i, "_combined_2map.png"),
  dpi = 800, limitsize = FALSE, width = 16, height = 6
)

# Save
save(cds, file = paste0(prefixOut, "monocle3_inferred_pseudotime.RData"))

# Plotting Marker Genes

# HPO in Monocyte Like cells
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "MPO", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
  filename = paste0(prefixOut, i, "_MPO_Trend.png"),
  dpi = 600, limitsize = FALSE
)

# GATA1 in Dendritic like cells
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "GATA1", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
  filename = paste0(prefixOut, i, "_GATA1_Trend.png"),
  dpi = 600, limitsize = FALSE
)

# GATA2
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "GATA2", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
  filename = paste0(prefixOut, i, "_GATA2_Trend.png"),
  dpi = 600, limitsize = FALSE
)

# CD34 in stem cell like cells
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "CD34", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
  filename = paste0(prefixOut, i, "_CD34_Trend.png"),
  dpi = 600, limitsize = FALSE
)

#  IRF8
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "IRF8", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
  filename = paste0(prefixOut, i, "_IRF8_Trend.png"),
  dpi = 600, limitsize = FALSE
)

# EPOR in erythocyte
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "EPOR", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
  filename = paste0(prefixOut, i, "_EPOR_Trend.png"),
  dpi = 600, limitsize = FALSE
)

# Down in GMP
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "APOBEC3C", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
  filename = paste0(prefixOut, i, "_APOBEC3C_Trend.png"),
  dpi = 600, limitsize = FALSE
)


# UP in Lymhpoid
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "ACY3", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
       filename = paste0(prefixOut, i, "_ACY3_Trend.png"),
       dpi = 600, limitsize = FALSE
)

# UP in Lymhpoid
p <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name == "ACY3", ], color_cells_by = "predicted.celltype.l2")
ggsave(p,
       filename = paste0(prefixOut, i, "_ACY3_Trend.png"),
       dpi = 600, limitsize = FALSE
)
