##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: TI with Monocle3 ######
##################################

set.seed(007)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(parallel))

# Prefix
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
dirPath <- paste0(base_string, "analysis_public_data/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")


# Detect Root cells
source(paste0(helpScriptsDir, "detect_root_nodes().R"))

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec

# Run lapply
umaps.list <- mclapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # rep_i = "rep3"

  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    loop <- TRUE
    graph_prune <- FALSE
    root_cell <- "HSC"
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    graph_prune <- TRUE
    loop <- TRUE
    root_cell <- "LMPP"
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    loop <- FALSE
    graph_prune <- TRUE
    root_cell <- "HSC"
  }

  sob <- readRDS(file = paste0(dirPath, rep_i, "/", rep_i, "subSampled.RDS"))

  # Get var features
  var_features <- VariableFeatures(sob)

  # Subset
  counts <- sob@assays$RNA$scale.data

  # Subset counts
  counts <- counts[var_features, , drop = FALSE]

  cell_metadata <- sob@meta.data
  cell_metadata <- cell_metadata[colnames(counts), ]
  gene_metadata <- data.frame(
    row.names = rownames(sob),
    gene_short_name = rownames(sob)
  )
  gene_metadata <- gene_metadata[rownames(counts), , drop = F]
  # Create cds
  cds <- new_cell_data_set(
    expression_data = counts,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
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
  cds <- cluster_cells(cds)

  # Learn graph
  cds <- learn_graph(cds,
    close_loop = F,
    use_partition = F,
    learn_graph_control = list(
      prune_graph = F
    )
  )

  # Order cells
  cds <- order_cells(cds,
    root_pr_nodes = find_root_pp(cds,
      cell = root_cell,
      cell_col = "cell_type"
    )
  )

  # Plot
  pseudotime <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1.5) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8)) + ggtitle(paste(
      individual, "| Age:", age,
      "| sex:", sex
    )) + xlab("UMAP-1") + ylab("UMAP-2")
  cell_type <- plot_cells(cds,
    color_cells_by = "cell_type", cell_size = 1.5,
    label_cell_groups = F
  ) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8)) + ggtitle(paste(
      individual, "| Age:", age,
      "| sex:", sex
    )) + xlab("UMAP-1") + ylab("UMAP-2")
  pseudotime
  cell_type

  # Save
  file_name <- paste0(outPath, rep_i, "/", paste(rep_i, "cds.RDS", sep = "_"))
  # saveRDS(object = cds, file = file_name)

  return(list(
    pseudotime = pseudotime,
    cell_type = cell_type
  ))
}, mc.cores = 3)

bottom <- ggarrange(umaps.list$rep1$pseudotime,
  umaps.list$rep2$pseudotime,
  umaps.list$rep3$pseudotime,
  labels = c("D.", "E.", "F."), nrow = 1,
  common.legend = F, legend = "bottom"
)
bottom
top <- ggarrange(umaps.list$rep1$cell_type,
  umaps.list$rep2$cell_type,
  umaps.list$rep3$cell_type,
  labels = c("A.", "B.", "C."), nrow = 1,
  common.legend = F, legend = "bottom"
)

combined_plot <- ggarrange(top, bottom, nrow = 2)
combined_plot

ggsave(
  plot = combined_plot, filename = "05_Real_Data_TI.png",
  path = figPath_hd,
  width = 14, dpi = 600, height = 8,
  background = "white"
)

ggsave(
  plot = combined_plot, filename = "05_Real_Data_TI.png",
  path = figPath_lr,
  width = 14, dpi = 150, height = 8,
  background = "white"
)
