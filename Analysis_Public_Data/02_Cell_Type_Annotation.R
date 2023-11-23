###########################################
## Author: Priyansh Srivastava ############
## Email: spriyansh29@gmail.com ###########
## Script: Automatic Cell type Inference ##
###########################################

set.seed(007)

suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))

# Get reference GSE24759 (Novershtern et al. 2011).
HSPCs <- NovershternHematopoieticData()

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
names(rep_vec) <- rep_vec

# Run lapply
umaps.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # Load seurat object
  sob <- LoadH5Seurat(file = paste0(inPath, rep_i, "/", rep_i, "_prs.h5seurat"))

  # Create a single cell experiment object
  sce <- SingleCellExperiment(list(logcounts = sob@assays$RNA@scale.data))

  # Assign labels
  pred.sce <- SingleR(
    test = sce, ref = HSPCs,
    assay.type.test = "logcounts",
    labels = HSPCs$label.main
  )

  # Add the annotations to the seurat object
  sob@meta.data$main_labels <- pred.sce$labels
  sob@meta.data$main_label_scores <- pred.sce$scores

  # Keep the annotation with high score
  sob <- subset(sob, main_labels %in% c("CMPs", "Dendritic cells", "Erythroid cells", "GMPs", "HSCs", "Megakaryocytes", "MEPs", "Monocytes"))

  # Rerun single R with fine labels
  sce <- SingleCellExperiment(list(logcounts = sob@assays$RNA@scale.data))
  # Assign fine labels
  pred.sce <- SingleR(
    test = sce, ref = HSPCs,
    assay.type.test = "logcounts",
    labels = HSPCs$label.fine
  )

  # Add the annotations to the seurat object
  sob@meta.data$fine_labels <- pred.sce$labels
  sob@meta.data$fine_label_scores <- pred.sce$scores

  # Keep the annotation with with Cd34+ cells
  sob <- subset(sob, fine_labels %in% c(
    "Common myeloid progenitors", "Early B cells",
    "Erythroid_CD34+ CD71+ GlyA-",
    "Hematopoietic stem cells_CD133+ CD34dim",
    "Hematopoietic stem cells_CD38- CD34+",
    "Megakaryocyte/erythroid progenitors",
    "Megakaryocytes", "Monocytes",
    "Myeloid Dendritic Cells",
    "Plasmacytoid Dendritic Cells", "Pro B cells "
  ))

  # Keep cells with high scores
  sob <- subset(sob, fine_label_scores >= -0.05)

  # Downstream processibg
  sob <- RunPCA(sob, features = VariableFeatures(object = sob), verbose = F)
  sob <- FindNeighbors(sob, dims = 1:10, verbose = F)
  sob <- FindClusters(sob, resolution = 1, verbose = F)
  sob <- RunUMAP(sob, dims = 1:10, verbose = F)

  # Save
  file_name <- paste(outPath, rep_i, paste(rep_i, "anno", sep = "_"), sep = "/")
  SaveH5Seurat(
    object = sob, filename = file_name,
    overwrite = T, verbose = FALSE
  )

  # Return UMAP
  return(NULL)
})
