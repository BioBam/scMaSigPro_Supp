# Title: Simulate 4 Datasets with Different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(phateR))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(Seurat))

# Set paths
paramEstimates <- readRDS("/supp_data/benchmarks/00_Parameter_Estimation/output/setty_et_al_d1_splatEstimates.RDS")
outDir <- "/supp_data/benchmarks/06_Offset/simulated/sce/"
helpScriptsDir <- "R_Scripts/helper_function/"
imgPath <- paste0(outDir, "png/")
sce_path <- paste0(outDir, "sce/")

# Create Directories
dir.create(outDir, showWarnings = F, recursive = T)
dir.create(imgPath, recursive = T, showWarnings = F)
dir.create(sce_path, recursive = T, showWarnings = F)

# Load Custom Functions
source(paste0(helpScriptsDir, "plot_simulations().R"))
source(paste0(helpScriptsDir, "add_gene_anno().R"))
source(paste0(helpScriptsDir, "calc_bin_size.R"))
source(paste0(helpScriptsDir, "calcNormCounts.R"))

# Zero-Inflation
"sparsity_60" <- -0.5

# Create Base parameters/ Same for All groups
params.groups <- newSplatParams(
  batch.rmEffect = TRUE, # No Batch affect
  batchCells = 3000, # Number of Cells
  nGenes = 2000, # Number of Genes
  seed = 2022, # Set seed
  mean.rate = paramEstimates@mean.rate,
  mean.shape = paramEstimates@mean.shape,
  lib.scale = paramEstimates@lib.scale,
  lib.loc = paramEstimates@lib.loc,
  bcv.common = paramEstimates@bcv.common,
  bcv.df = paramEstimates@bcv.df,
  dropout.type = "experiment",
  group.prob = c(0.5, 0.5),
  path.from = c(0, 0),
  de.prob = 0.3,
  de.facLoc = 1,
  out.facLoc = paramEstimates@out.facLoc,
  dropout.mid = paramEstimates@dropout.mid,
  out.facScale = paramEstimates@out.facScale,
  out.prob = paramEstimates@out.prob,
  path.skew = c(0.5, 0.5),
  path.nSteps = c(1500, 1500)
)

# Simulate Object
sim.sce <- splatSimulate(
  params = params.groups,
  method = "paths",
  verbose = F,
  dropout.shape = -0.5 # 60%
)

# Sparsity values
trueSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
simulatedSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
totSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)

# Add gene Info
gene.info <- add_gene_anno(sim.sce = sim.sce)
gene.info <- gene.info[mixedsort(gene.info$gene_short_name), ]

# Update the SCE Simulated Object
rowData(sim.sce) <- DataFrame(gene.info)

# SaveRDS
obj.path <- paste0(sce_path, paste0("zi_", totSparsity, ".RData"))
save(sim.sce, file = obj.path)

# Names
label_vector <- c(
  "Total_Sparsity" = totSparsity,
  "True_Sparsity" = trueSparsity,
  "Simulated_Sparsity" = simulatedSparsity,
  "Filename" = paste0("zi_", totSparsity, ".RData")
)

# Create Normlization vector
norm <- c(
  # "trueCounts",
  "rawCounts",
  "SCT",
  "libSize",
  "logLibSize",
  "rawCounts_Offset",
  # "FQ_Offset",
  "FQ"
)

# Generate Datasets
parameter.list <- mclapply(norm, function(norm_type, sce.path = sce_path,
                                          sim_sce = sim.sce) {
  # Extract Counts
  rawCounts <- sim_sce@assays@data@listData$counts
  trueCounts <- sim_sce@assays@data@listData$TrueCounts

  # Check norm type
  if (norm_type == "SCT") {
    cat(paste("Performing Seurat's SCT"))

    # Compute
    prs_counts <- calcNormCounts(
      rawCounts = rawCounts,
      cat = "sctransform"
    )
  } else if (norm_type == "libSize") {
    cat(paste("Performing Seurat's RC"))

    # Compute
    prs_counts <- calcNormCounts(
      rawCounts = rawCounts,
      cat = "libSize"
    )
  } else if (norm_type == "logLibSize") {
    cat(paste("Performing Seurat's RC"))

    # Compute
    prs_counts <- calcNormCounts(
      rawCounts = rawCounts,
      cat = "logLibSize"
    )
  } else if (norm_type == "rawCounts" | norm_type == "rawCounts_Offset") {
    cat(paste("Using Raw Counts"))
    # Compute
    prs_counts <- rawCounts
  } else if (norm_type == "FQ") {
    cat(paste("Performing FQ"))
    # Compute
    prs_counts <- calcNormCounts(
      rawCounts = rawCounts,
      cat = "FQNorm"
    )
  }

  # Update the Sce.object
  sim_sce@assays@data@listData$counts <- as(as.matrix(prs_counts), "dgCMatrix")

  # Save RDS
  file_name <- paste0(outDir, "norm_", norm_type, ".RData")
  save(sim_sce, file = file_name)


  # Return
  return(NULL)
}, mc.cores = 16)
