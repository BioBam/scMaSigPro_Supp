##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: Azimuth Annotation + Monocle3 Input Prepare ###
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

# Prefix
prefixIn <- "benchmarks/11_RealDataSmall/data/results/"
prefixOut <- "benchmarks/11_RealDataSmall/data/results/"

# Set Variables
i <- "rep1"
individual <- "1"
age <- "35"
sex <- "Male"

#  the Processed Seurat Object
sob <- LoadH5Seurat(file = paste0(prefixIn, "rep1_Processed_sob.h5seurat"), verbose = F)

# Azimuth Annotations
bm <- RunAzimuth(
  query = sob,
  reference = paste0(
    str_sub(prefixIn, end = -9),
    "input/Azimuth_Human_BoneMarrow"
  )
)

# Plot the Annotations
p <- DimPlot(bm, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_hue(l = 50) + theme(legend.position = "bottom")

ggsave(p,
  filename = paste0(prefixOut, i, "_All_Anno_Azimuth.png"),
  dpi = 1400, limitsize = FALSE, width = 8, height = 8
)

# Subset the Scores
bm.score <- subset(bm, subset = predicted.celltype.l2.score >= 0.6)

# Plot the Annotations
p <- DimPlot(bm.score, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_hue(l = 50) + theme(legend.position = "bottom")

# Save
ggsave(p,
  filename = paste0(prefixOut, i, "_Score_Sub_Anno_Azimuth.png"),
  dpi = 1400, limitsize = FALSE, width = 8, height = 8
)

# Subset the Cell types
bm.cell <- subset(bm.score,
  subset = predicted.celltype.l2 %in% c(
    "Early Eryth",
    "GMP", "HSC",
    "LMPP", "Late Eryth",
    "CLP", "EMP"
  )
)

# Plot the Annotations
p <- DimPlot(bm.cell, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_hue(l = 50) + theme(legend.position = "bottom")

ggsave(p,
  filename = paste0(prefixOut, i, "_Cell_Sub_Anno_Azimuth.png"),
  dpi = 1400, limitsize = FALSE, width = 8, height = 8
)

# Plot Clusters
p <- DimPlot(bm.cell, reduction = "umap", pt.size = 1, group.by = "seurat_clusters") +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_hue(l = 50) + theme(legend.position = "bottom")

ggsave(p,
  filename = paste0(prefixOut, i, "_Seurat_Clusters.png"),
  dpi = 1400, limitsize = FALSE, width = 8, height = 8
)

# Write Seurat H5
file_name <- paste0(prefixOut, i, "_Annotated_sob")
SaveH5Seurat(
  object = bm.cell, filename = file_name,
  overwrite = T, verbose = FALSE
)
