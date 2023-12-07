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
suppressPackageStartupMessages(library(diffusionMap))
suppressPackageStartupMessages(library(destiny))

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
    subset.vector <- c("EMP", "Early Eryth", "Prog Mk")
    min_dist = 0.001
    neigh =  200
    dim = c(1:10)
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    subset.vector <- c("GMP", "CLP", "HSC", "LMPP")
    min_dist = 0.1
    neigh =  10
    dim = c(1:5)
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    subset.vector <- c("HSC", "CLP", "GMP", "LMPP")
    min_dist = 0.001
    neigh =  200
    dim = c(1:5)
  }
  # Load seurat object
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_azimuth.RDS"))

  # Subset
  sob <- subset(sob, cell_type %in% c("pre-pDC", "GMP", "Early Eryth", "HSC",
                                      "LMPP", "EMP", "pDC", "Late Eryth", 
                                      "Prog Mk", "CLP", "cDC2", "pre-mDC", "ASDC"))
  sob.sub <- subset(sob, predicted.celltype.l2.score >= 0.2)
  
  # Compute PCA
  sob.sub <- RunPCA(sob.sub, npcs = 100, verbose = F)
  
  # Compute tsne
  sob.sub <- RunTSNE(sob.sub,
                     reduction = "pca",
                     dim.embed = 2)
  
  # Plot
  plt <- DimPlot(sob.sub, group.by = "cell_type", reduction = "tsne") + xlab("tSNE-1") + ylab("tSNE-2")+ ggtitle(paste(
      individual, "| Age:", age,
      "| sex:", sex
  )) + theme(legend.position = "bottom", legend.justification = "center")

  plt
  
  file_name <- paste0(outPath, rep_i, "/", rep_i, "subSampled.RDS")
  saveRDS(sob.sub, file_name)
  
  # Return
  return(list(all = plt))
}, mc.cores = 1)
names(azimuth.list) <- rep_vec

# Create plot
top <- ggarrange(azimuth.list$rep1$all,
                 azimuth.list$rep2$all,
                 azimuth.list$rep3$all,
                    labels = c("A.", "B.", "C."), nrow = 1,
                    common.legend = F, legend = "bottom"
)
top
# bottom <- ggarrange(azimuth.list$rep1$sub,
#                     azimuth.list$rep2$sub,
#                     azimuth.list$rep3$sub,
#                  labels = c("D.", "E.", "F."), nrow = 1,
#                  common.legend = F, legend = "bottom"
# )
# 
# bottom2 <- ggarrange(azimuth.list$rep1$sub.sub,
#                     azimuth.list$rep2$sub.sub,
#                     azimuth.list$rep3$sub.sub,
#                     labels = c("G.", "H.", "I."), nrow = 1,
#                     common.legend = F, legend = "bottom"
# )
# 
# combined_plot <- ggarrange(top, bottom, bottom2, nrow = 3)
# combined_plot

ggsave(top,
       filename = paste0("Figures/SuppData/05_Real_Data-SubSampling.png"),
       dpi = 150, height = 8, width = 12
)
