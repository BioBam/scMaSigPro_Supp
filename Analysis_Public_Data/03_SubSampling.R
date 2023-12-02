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

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Load data
azimuth.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
  }
  # Load seurat object
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_azimuth.RDS"))


  # Subset
  sob.sub <- subset(sob, main_labels %in% c(
      "HSC", "GMP", "EMP", "CLP", "LMPP", "Prog Mk", "Early Eryth", "pDc", "pre B",
      "pre-mDC", "pre-pDC", "pro B"))

  # Recompute
  sob.sub <- RunPCA(sob.sub, features = VariableFeatures(object = sob.sub), verbose = F)
  sob.sub <- FindNeighbors(sob.sub, verbose = F)
  sob.sub <- FindClusters(sob.sub, resolution = 1, verbose = F)

  # Compute UMAP
  sob.sub <- RunUMAP(sob.sub,
    verbose = F, features = VariableFeatures(sob.sub)
  )

  # Plot
  plt <- DimPlot(sob.sub, group.by = "cell_type") + ggtitle(paste(
    individual, "| Age:", age,
    "| sex:", sex
  )) + theme(legend.position = "bottom", legend.justification = "center")

  # Return
  return(list(
    sob.sub = sob.sub,
    plt = plt
  ))
})

# Save
nullList <- lapply(names(azimuth.list), function(rep_i) {
  ob.list <- azimuth.list[[rep_i]]$sob.sub
  file_name <- paste0(dirPath, rep_i, "/", rep_i, "subSampled.RDS")
  saveRDS(ob.list, file_name)
})
