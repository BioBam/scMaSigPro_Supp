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
  sob.sub <- sob

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
  )) +theme(legend.position = "bottom", legend.justification = "center")

  # Return
  return(list(sob.sub = sob.sub,
              plt = plt))
})

# Subsample for donor-1
all.donor.subSample.list <- lapply(
  names(azimuth.list),
  function(rep_i, inPath = dirPath, outPath = dirPath) {
      
      sob <- azimuth.list[[rep_i]][["sob.sub"]]
      
      # Step-1: Add Annotation for donors
      if (rep_i == "rep1") {
          individual <- "Donor-1"
          age <- "35"
          sex <- "Male"
          min_dist = 0.1
          n_neighbors = 20
          sp = 1
          
          # Keep Cells
          subSample <- c("ASDC", "HSC","EMP", "GMP", "Macrophage", "cDC2", "pre-mDC", "pre-pDC",
                         "EMP", "Prog Mk", "Platelet", "Early Eryth", "pro B")
          
      } else if (rep_i == "rep2") {
          individual <- "Donor-2"
          age <- "28"
          sex <- "Female"
          n_neighbors = 20
          min_dist = 0.1
          sp = 1
          
          # Keep Cells
          subSample <- c("ASDC", "HSC","EMP", "GMP", "Macrophage", "cDC2", "pre-mDC", "pre-pDC",
                         "EMP", "Early Eryth")
          
      } else if (rep_i == "rep3") {
          individual <- "Donor-3"
          age <- "19"
          sex <- "Female"
          min_dist = 0.2
          n_neighbors = 50
          sp = 0.5
          
          # Keep Cells
          subSample <- c("ASDC", "HSC","EMP", "GMP", "Macrophage", "cDC2", "pre-mDC", "pre-pDC",
                         "EMP", "Prog Mk", "Platelet", "Early Eryth", "pro B")
          
      }

    # Get all cells
    all_cells <- unique(sob@meta.data$cell_type)

      # Keep cells
      keep <- all_cells[all_cells %in% subSample]

      # Subset
      sob.sub <- subset(sob, cell_type %in% keep)
      sob.sub <- subset(sob.sub, subset = nFeature_RNA > 50)

      # Recompute
      sob.sub <- RunPCA(sob.sub, features = VariableFeatures(object = sob.sub), verbose = F)
      sob.sub <- FindNeighbors(sob.sub, verbose = F)
      sob.sub <- FindClusters(sob.sub, resolution = 1, verbose = F)

      # Compute UMAP
      sob.sub <- RunUMAP(sob.sub,
        verbose = F, features = VariableFeatures(sob.sub),
        min.dist = min_dist,
        n.neighbors = n_neighbors,
        spread = sp
      )

      # Plot
      plt <- DimPlot(sob.sub, group.by = "cell_type") + ggtitle(paste(
        individual, "| Age:", age,
        "| sex:", sex
      ))+theme(legend.position = "bottom", legend.justification = "center")

    # Return UMAP
    return(list(obj = sob.sub,
                plot = plt))
  }
)

names(all.donor.subSample.list) <-  names(azimuth.list)

# Plot
sub_samples <- ggarrange(
    azimuth.list$rep1$plt,
    azimuth.list$rep2$plt,
    azimuth.list$rep3$plt,
    all.donor.subSample.list$rep1$plot,
    all.donor.subSample.list$rep2$plot,
    all.donor.subSample.list$rep3$plot,
    nrow=2, ncol = 3,common.legend = T,
    
legend = "bottom",
    labels = c("A.", "B.", "C.", "D.", "E.", "F.")
    )
sub_samples


ggsave(sub_samples,
       filename = "05_RealData-SubSampling.png",
       path = "Figures/SuppData",
       width = 12, dpi = 150, height = 8)

# Save
nullList <- lapply(names(all.donor.subSample.list), function(rep_i){
    ob.list <- all.donor.subSample.list[[rep_i]]$obj
        file_name <- paste0(dirPath, rep_i,"/", rep_i,"subSampled.RDS")
        saveRDS(ob.list, file_name)
})
