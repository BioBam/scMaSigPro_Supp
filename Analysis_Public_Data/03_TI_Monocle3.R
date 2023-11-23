##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: TI with Monocle3 ######
##################################

set.seed(007)

suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(monocle3))

# Detect Root cells
source("R_Scripts/helper_function/detect_root_nodes().R")

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
names(rep_vec) <- rep_vec

# Run lapply
umaps.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # Load seurat object
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_anno.RDS"))

  # Create cds
  cds <- new_cell_data_set(
    expression_data = sob@assays$RNA@data,
    cell_metadata = sob@meta.data,
    gene_metadata = data.frame(
      row.names = rownames(sob),
      gene_short_name = rownames(sob)
    )
  )

  # No Normalization
  cds <- preprocess_cds(cds, norm_method = "none")

  # Compute UMAP
  cds <- reduce_dimension(cds)

  # We will use UMAP and clusters from seurat
  new_umap <- as.matrix(sob@reductions$umap@cell.embeddings)
  colnames(new_umap) <- NULL
  reducedDims(cds)[["UMAP"]] <- new_umap

  # Compute clusters and use single partition
  cds <- cluster_cells(cds, resolution = 0.5)
  
  # sing.partition <- rep(1, length(cds@clusters$UMAP$partitions))
  # names(sing.partition) <- names(cds@clusters$UMAP$partitions)
  # cds@clusters$UMAP$partitions <- as.factor(sing.partition)

  # Learn graph
  cds <- learn_graph(cds)

  # Order cells
  cds@colData$cell_type <- cds@colData$fine_labels
  cds <- order_cells(cds,
    root_pr_nodes = find_root_pp(cds,
      cell = "Hematopoietic stem cells_CD133+ CD34dim",
      cell_col = "cell_type"
    )[1]
  )
  
  # Save
  file_name <- paste(outPath, rep_i, paste(rep_i, "cds.RDS", sep = "_"), sep = "/")
  saveRDS(
      object = cds, file = file_name)

  # Plot
  pseudotime <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1.5)+
      theme(legend.position = "bottom") + ggtitle(paste(rep_i))
  cell_type <- plot_cells(cds, color_cells_by = "fine_labels", cell_size = 1.5,
                          label_cell_groups = F)+
      theme(legend.position = "bottom") + ggtitle(paste(rep_i))

  return(list(pseudotime = pseudotime,
              cell_type = cell_type))
})

bottom <- ggarrange(umaps.list$rep1$pseudotime,
          umaps.list$rep2$pseudotime,
          umaps.list$rep3$pseudotime,
          labels = c("D.", "E.", "F."),nrow = 1,
          common.legend = F, legend = "bottom")
top <- ggarrange(umaps.list$rep1$cell_type,
                    umaps.list$rep2$cell_type,
                    umaps.list$rep3$cell_type,
                    labels = c("A.", "B.", "C."),nrow = 1,
                    common.legend = T, legend = "bottom")

combined_plot <- ggarrange(top, bottom, nrow = 2)
ggsave(combined_plot,
       filename = paste0("Figures/SuppData/01_Real_Data.png"),
       dpi = 600, height = 8, width = 16
)
