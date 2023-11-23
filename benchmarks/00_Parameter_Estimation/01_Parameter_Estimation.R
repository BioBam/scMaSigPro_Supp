# Title: Estimate Parameters with SplatEstimate
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Load packages
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(Seurat))

# Set paths
inPath <- "/supp_data/benchmarks/00_Parameter_Estimation/input/"
outPath <- "/supp_data/benchmarks/00_Parameter_Estimation/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load scripts
source(paste0(helpScriptsDir, "plot_simulations().R"))

# Load Counts donor-1
cr_counts <- Read10X_h5(paste0(inPath, "filtered_feature_bc_matrix.h5"))

# Create SeuratObject
sob <- CreateSeuratObject(counts = cr_counts,
                          project = "setty_et_all_don1",
                          min.cells = 2000, min.features = 2500)
# Basic QC
sob[["percent.mt"]] <- PercentageFeatureSet(sob, pattern = "^MT-")
sob <- subset(sob, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

# Convert to matrix for splatter
cr_counts_matrix <- as.matrix(sob@assays$RNA@counts)

# Estimate parameters
set.seed(2023)
splatParamObj <- splatEstimate(cr_counts_matrix)

# Save Object
saveRDS(object = splatParamObj, file = paste0(outPath, "setty_et_al_d1_splatEstimates.RDS"))