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
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Load data
azimuth.list <- mclapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    #subset.vector <- c("HSC", "LMPP", "CLP", "GMP")
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    #subset.vector <- c("HSC", "EMP", "GMP", "LMPP", "CLP")
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    #subset.vector <- c("HSC", "CLP", "GMP", "LMPP", "pre-pDC", "pre-mDC")
  }
  # Load seurat object
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_azimuth.RDS"))

  # Subset
  sob.sub <-subset(sob, cell_type %in% c(
      "HSC", "GMP", "EMP", "CLP", "LMPP", "Prog Mk", "Early Eryth", "pDc", "pre B",
      "pre-mDC", "pre-pDC", "pro B"))

  # Recompute
  sob.sub <- RunPCA(sob.sub, features = VariableFeatures(object = sob.sub), verbose = F)
  sob.sub <- FindNeighbors(sob.sub, verbose = F)
  sob.sub <- FindClusters(sob.sub, resolution = 1, verbose = F)

  # Compute UMAP
  sob.sub <- RunUMAP(sob.sub,
                     #umap.method = "umap-learn",
                     min.dist = 0.01,
                     n.neighbors = 200, 
    verbose = F, dims = c(1:5) #features = VariableFeatures(sob.sub)
  )

  # Plot
  plt <- DimPlot(sob, group.by = "cell_type") + ggtitle(paste(
      individual, "| Age:", age,
      "| sex:", sex
  )) + theme(legend.position = "bottom", legend.justification = "center")
  
  plt.sub <- DimPlot(sob.sub, group.by = "cell_type") + ggtitle(paste(
    individual, "| Age:", age,
    "| sex:", sex
  )) + theme(legend.position = "bottom", legend.justification = "center")
  plt.sub

  file_name <- paste0(outPath, rep_i, "/", rep_i, "subSampled.RDS")
  saveRDS(sob.sub, file_name)
  
  # Return
  return(list(all = plt,
              sub = plt.sub))
}, mc.cores = 24)
names(azimuth.list) <- rep_vec

top <- ggarrange(azimuth.list$rep1$all,
                 azimuth.list$rep2$all,
                 azimuth.list$rep3$all,
                    labels = c("A.", "B.", "C."), nrow = 1,
                    common.legend = F, legend = "bottom"
)
top
bottom <- ggarrange(azimuth.list$rep1$sub,
                    azimuth.list$rep2$sub,
                    azimuth.list$rep3$sub,
                 labels = c("D.", "E.", "F."), nrow = 1,
                 common.legend = F, legend = "bottom"
)

combined_plot <- ggarrange(top, bottom, nrow = 2)
combined_plot

ggsave(combined_plot,
       filename = paste0("Figures/SuppData/05_Real_Data-SubSampling.png"),
       dpi = 150, height = 8, width = 12
)
