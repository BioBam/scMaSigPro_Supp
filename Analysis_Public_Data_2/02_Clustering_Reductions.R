#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Cluster and Dimensions #####
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(diffusionMap))
suppressPackageStartupMessages(library(phateR))
suppressPackageStartupMessages(library(reticulate))

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data_2"

# Load Python
use_python("/tmp/RtmpmKsIgq/rstudio/terminal/python3", required = TRUE)

# Load leidenalg
leidenalg <- import("leidenalg")
phate <- import("phate")

# Load Seurat Object
sob_prs <- LoadH5Seurat(paste(dirPath,"input_seurat_object.h5seurat",sep="/"))

## Plot elbow
ElbowPlot(sob_prs, ndims = 100)
# Selct first 50 componenets

## Get Variance Explained by top 50 components
mat <- GetAssayData(sob_prs, layer = "RNA", slot = "scale.data")
pca <- sob_prs[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
sum(varExplained[c(1:50)] * 100)

# free up space
mat <- NULL
gc()

# Create Graph
sob_prs <- FindNeighbors(sob_prs, dims = 1:50)

# Compute Clusters with 
sob_prs <- FindClusters(sob_prs, resolution = 0.1)

# Compute tSNE
#sob_prs <- RunTSNE(sob_prs, dims = 1:50, perplexity = 500)

#DimPlot(sob_prs, reduction = "tsne") + DimPlot(sob_prs, reduction = "tsne", group.by = "sample")

# # Now we will compute PHATE dimensions over PCA
pca_matrix <- sob_prs@reductions$pca@cell.embeddings
pca_matrix <- pca_matrix[,c(1:5)]

# Compute distance matrix
pca_matrix_dist <- dist(pca_matrix)
sob_prs <-NULL
gc()

# Compute DiffusionMap
diffusion_ob <- diffuse(pca_matrix_dist, maxdim = 2, neigen = 2,
                        delta = 10^-1)
# # 
# # Compute Phate
# phate_ob <- phate(pca_matrix, n.jobs=20)
# 
# # Extract Phate embeddings
# phate_matrix <- phate_ob$embedding

