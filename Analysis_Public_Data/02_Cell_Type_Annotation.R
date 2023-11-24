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
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(ggpubr))

# Get reference GSE24759 (Novershtern et al. 2011).
HSPCs <- NovershternHematopoieticData()

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Run lapply
azimuth.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # Load seurat object
  sob <- LoadH5Seurat(file = paste0(inPath, rep_i, "/", rep_i, "_prs.h5seurat"))
  
  sob <- RunAzimuth(
      query = sob,
      reference = paste0(inPath, "Azimuth_Human_BoneMarrow"
      )
  )
  
  # Add column
  sob@meta.data$cell_type <- sob@meta.data$predicted.celltype.l2
  return(sob)
})

# Run lapply
null.list <- lapply(azimuth.list, function(sob, inPath = dirPath, outPath = dirPath) {
    
  # Select cells with relatively high score
  sob.sub <- subset(sob, predicted.celltype.l2.score > .2)
  
  # Keep cell type of interest
  sob.sub <- subset(sob.sub, cell_type %in% c("CLP", "EMP", "GMP", "HSC",
                               "LMPP",
                               "pre-mDC", "pre-pDC", "Prog Mk"))
  
  sob.sub <- RunPCA(sob.sub, features = VariableFeatures(object = sob.sub), verbose = F)
  sob.sub <- FindNeighbors(sob.sub, dims = 1:5, verbose = F)
  sob.sub <- FindClusters(sob.sub, resolution = 1, verbose = F)
  sob.sub <- RunUMAP(sob.sub, dims = 1:5, verbose = F, n.neighbors = 25,
                     min.dist = 0.0001)
  
  # Save
  # file_name <- paste(outPath, rep_i, paste(rep_i, "anno.RDS", sep = "_"), sep = "/")
  # saveRDS(
  #   object = sob.sub, file = file_name)

  # Return UMAP
  return(DimPlot(sob.sub, group.by = "cell_type"))
})

ggarrange(null.list[[1]], null.list[[1]], null.list[[1]], ncol =1)
# SingleR Clusters
#   # Create a single cell experiment object
#   sce <- SingleCellExperiment(list(logcounts = sob@assays$RNA@scale.data))
# 
#   # Assign labels
#   pred.sce <- SingleR(
#     test = sce, ref = HSPCs,
#     assay.type.test = "logcounts",
#     labels = HSPCs$label.main
#   )
# 
#   # Add the annotations to the seurat object
#   sob@meta.data$main_labels <- pred.sce$labels
#   sob@meta.data$main_label_scores <- pred.sce$scores
# 
#   # Keep the annotation with high score
#   sob <- subset(sob, main_labels %in% c("CMPs", "Dendritic cells", "Erythroid cells", "GMPs", "HSCs", "Megakaryocytes", "MEPs", "Monocytes"))

# # Rerun single R with fine labels
# sce <- SingleCellExperiment(list(logcounts = sob@assays$RNA@scale.data))
# # Assign fine labels
# pred.sce <- SingleR(
#   test = sce, ref = HSPCs,
#   assay.type.test = "logcounts",
#   labels = HSPCs$label.fine
# )
# 
# # Add the annotations to the seurat object
# sob@meta.data$fine_labels <- pred.sce$labels
# sob@meta.data$fine_label_scores <- pred.sce$scores
# 
# # Keep the annotation with with Cd34+ cells
# sob <- subset(sob, fine_labels %in% c(
#   "Common myeloid progenitors", 
#   #"Early B cells",
#   "Erythroid_CD34+ CD71+ GlyA-",
#   "Hematopoietic stem cells_CD133+ CD34dim",
#   "Hematopoietic stem cells_CD38- CD34+",
#   "Megakaryocyte/erythroid progenitors",
#   #"Megakaryocytes", 
#   "Monocytes",
#   "Myeloid Dendritic Cells"#,
#   #"Plasmacytoid Dendritic Cells",
#   #"Pro B cells"
# ))
# 
# # Keep cells with high scores
# #sob <- subset(sob, fine_label_scores >= -0.05)