###########################################
## Author: Priyansh Srivastava ############
## Email: spriyansh29@gmail.com ###########
## Script: Automatic Cell type Inference ##
###########################################

set.seed(007)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(Azimuth))

# Prefix
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
dirPath <- paste0(base_string, "analysis_public_data/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec

# Run lapply
azimuth.list <- mclapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # Load seurat object
  sob <- readRDS(paste0(inPath, rep_i, "/", paste0(rep_i, "_prs.RDS")))

  sob <- RunAzimuth(
    query = sob,
    reference = paste0(inPath, "Azimuth_Human_BoneMarrow")
  )

  # Add column
  sob@meta.data$cell_type <- sob@meta.data$predicted.celltype.l2

  # Save
  file_name <- paste0(outPath, rep_i, "/", paste(rep_i, "azimuth.RDS", sep = "_"))
  saveRDS(
    object = sob, file = file_name
  )

  return(sob)
}, mc.cores = 1)
