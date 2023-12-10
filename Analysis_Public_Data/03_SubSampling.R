###########################################
## Author: Priyansh Srivastava ############
## Email: spriyansh29@gmail.com ###########
## Script: Subsampling ####################
###########################################

set.seed(007)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(parallel))

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec[3]

# Load data
azimuth.list <- mclapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
    #rep_i <- "rep3"
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    dims = c(1:50)
    seed.use = 123
    min.dist = 0.1
    n.neighbors = 200
    spread = 10
    nfeatures = 6000
    npcs =100
    select_cells = c("GMP", "Early Eryth", "HSC", "LMPP", "EMP", "Prog Mk", "CLP",
                     "pre-pDC", "pre-mDC", "pro B")
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    dims = c(1:50)
    seed.use = 123
    min.dist = 0.01
    n.neighbors = 200
    spread = 5
    npcs = 100
    nfeatures = 6000
    select_cells = c("LMPP", "pre-pDC", "pre-mDC","CLP")
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    dims = c(1:25)
    min.dist = 0.1
    n.neighbors = 200
    spread = 10
    nfeatures = 6000
    npcs =100
    select_cells = c("GMP", "Early Eryth", "HSC", "LMPP", "EMP", "Prog Mk", "CLP",
                     "pre-pDC", "pre-mDC", "pro B")
  }
  # Load seurat object
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_azimuth.RDS"))
  saveRDS(rownames(sob), paste0(outPath, rep_i, "/", rep_i, "background.RDS"))

  # Subset
  sob.sub <- subset(sob, cell_type %in% select_cells)
  
  # Find to 15000 Varible features
  sob.sub <- FindVariableFeatures(sob.sub, 
                                  nfeatures = 6000, 
                                  verbose = F)
  # Compute PCA
  sob.sub <- RunPCA(sob.sub, npcs = 100, verbose = F)
  
  # Compute tsne
  sob.sub <- RunUMAP(sob.sub, verbose = T,
                     reduction = "pca",
                     dims = dims,
                     seed.use = 123,
                     min.dist = min.dist,
                     n.neighbors = n.neighbors,
                     spread = spread
                     )
  
  # Plot
  plt <- DimPlot(sob.sub, group.by = "cell_type", reduction = "umap") + xlab("UMAP-1") + ylab("UMAP-2")+ ggtitle(paste(
      individual, "| Age:", age,
      "| sex:", sex
  )) + theme(legend.position = "bottom", legend.justification = "center")
  plt
  
  file_name <- paste0(outPath, rep_i, "/", rep_i, "subSampled.RDS")
  saveRDS(sob.sub, file_name)
  
  # Return
  return(plt)
}, mc.cores = 16)
names(azimuth.list) <- rep_vec

# Create plot
top <- ggarrange(azimuth.list$rep1,
                 azimuth.list$rep2,
                 azimuth.list$rep3,
                    labels = c("A.", "B.", "C."), nrow = 1,
                    common.legend = F, legend = "bottom"
)
top

ggsave(top,
       filename = paste0("Figures/SuppData/05_Real_Data_SubSampling.png"),
       dpi = 150, height = 8, width = 12
)
