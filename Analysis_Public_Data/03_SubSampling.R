###########################################
## Author: Priyansh Srivastava ############
## Email: spriyansh29@gmail.com ###########
## Script: Subsampling ####################
###########################################

set.seed(007)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
# suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(parallel))

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec

# Load data
azimuth.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # rep_i <- "rep3"
  # inPath = dirPath
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    dims <- c(1:50)
    seed.use <- 123
    min.dist <- 0.1
    n.neighbors <- 200
    spread <- 10
    nfeatures <- 6000
    npcs <- 100
    select_cells <- c(
      "GMP", "Early Eryth", "HSC", "LMPP", "EMP", "Prog Mk", "CLP",
      "pre-pDC", "pre-mDC", "pro B"
    )
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    dims <- c(1:50)
    seed.use <- 123
    min.dist <- 0.01
    n.neighbors <- 200
    spread <- 5
    npcs <- 100
    nfeatures <- 6000
    select_cells <- c("LMPP", "pre-pDC", "pre-mDC", "CLP")
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    dims <- c(1:25)
    min.dist <- 0.1
    n.neighbors <- 200
    spread <- 10
    nfeatures <- 6000
    npcs <- 100
    select_cells <- c(
      "GMP", "Early Eryth", "HSC", "LMPP", "EMP", "Prog Mk", "CLP",
      "pre-pDC", "pre-mDC", "pro B"
    )
  }
  print(rep_i)
  # # Load seurat object
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_azimuth.RDS"))
  saveRDS(rownames(sob), paste0(outPath, rep_i, "/", rep_i, "background.RDS"))

  # # Subset
  sob.sub <- subset(sob, cell_type %in% select_cells)

  # # Find to 15000 Varible features
  sob.sub <- FindVariableFeatures(sob.sub,
    nfeatures = 6000,
    verbose = F
  )
  # # Compute PCA
  sob.sub <- RunPCA(sob.sub, npcs = 100, verbose = F)
  #
  # # Compute tsne
  sob.sub <- RunUMAP(sob.sub,
    verbose = T,
    reduction = "pca",
    dims = dims,
    seed.use = 123,
    min.dist = min.dist,
    n.neighbors = n.neighbors,
    spread = spread
  )
  #
  # # Plot
  plt <- DimPlot(sob.sub, group.by = "cell_type", reduction = "umap") + xlab("UMAP-1") + ylab("UMAP-2") + ggtitle(paste(
    individual, "| Age:", age,
    "| sex:", sex
  )) + theme(
    legend.position = "bottom", legend.justification = "center",
    legend.text = element_text(size = 6)
  )
  plt
  #
  # # Compute
  sob <- RunPCA(sob, verbose = F)
  sob <- RunUMAP(sob, dim = c(1:10))
  # Plot
  plt.all <- DimPlot(sob, group.by = "cell_type", reduction = "umap") + xlab("UMAP-1") + ylab("UMAP-2") + ggtitle(paste(
    individual, "| Age:", age,
    "| sex:", sex
  )) + theme(
    legend.position = "bottom", legend.justification = "center",
    legend.text = element_text(size = 6)
  )
  plt.all

  file_name <- paste0(outPath, rep_i, "/", rep_i, "subSampled.RDS")
  saveRDS(sob.sub, file_name)

  # Return
  return(list(
    plt = plt,
    plt.all = plt.all
  ))
})
names(azimuth.list) <- rep_vec

# Create plot
bottom <- ggarrange(azimuth.list$rep1$plt,
  azimuth.list$rep2$plt,
  azimuth.list$rep3$plt,
  labels = c("D.", "E.", "F."), nrow = 1,
  common.legend = F, legend = "bottom"
)
bottom

# Create plot
top <- ggarrange(azimuth.list$rep1$plt.all,
  azimuth.list$rep2$plt.all,
  azimuth.list$rep3$plt.all,
  labels = c("A.", "B.", "C."), nrow = 1,
  common.legend = F, legend = "bottom"
)
top

combined <- ggarrange(top, bottom,
  nrow = 2
)
combined

ggsave(combined,
  filename = paste0("/supp_data/Figures/SuppData/05_Real_Data_SubSampling.png"),
  dpi = 300, height = 12, width = 16
)
