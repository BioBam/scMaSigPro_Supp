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
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(Seurat))

# Set paths
paramEstimates <- readRDS("/supp_data/benchmarks/00_Parameter_Estimation/output/setty_et_al_d1_splatEstimates.RDS")
outDir <- "/supp_data/benchmarks/01_Sparsity/simulated/"
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

# Zero-Inflation
zi <- list(
  #  "sparsity_50" = -2, # 50 % #
  #  "sparsity_55" = -1, # 55 % #
  "sparsity_60" = -0.5, # 60 % #
  #  "sparsity_65" = -0.2, # 65 % #
  "sparsity_70" = 0.02, # 70 % #
  #  "sparsity_75" = 0.3, # 75 % #
  "sparsity_80" = 0.7, # 80 %#
  #  "sparsity_85" = 1.2, # 85 % #
  "sparsity_90" = 2.5 # 90 % #
)

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

# Generate Datasets
parameter.list <- mclapply(names(zi), function(dropout_shape, params_groups = params.groups,
                                               sce.path = sce_path) {
  # Get Variables
  total_sparsity <- str_remove(pattern = "sparsity_", dropout_shape)
  dropout_shape_value <- zi[[dropout_shape]]

  # Simulate Object
  sim.sce <- splatSimulate(
    params = params_groups,
    method = "paths",
    verbose = F,
    dropout.shape = dropout_shape_value
  )

  # Sparsity values
  trueSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
  simulatedSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
  totSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)

  # Add gene Info
  gene.info <- add_gene_anno(sim.sce = sim.sce)
  gene.info <- gene.info[mixedsort(gene.info$gene_short_name), ]

  # print(nrow(gene.info[((gene.info[["DEFacPath1"]] != 1 | gene.info[["DEFacPath2"]] != 1) & gene.info[["BaseGeneMean"]] >= 1), , drop=FALSE]))

  # Update the SCE Simulated Object
  rowData(sim.sce) <- DataFrame(gene.info)

  # SaveRDS
  obj.path <- paste0(sce.path, paste0("zi_", totSparsity, ".RData"))
  save(sim.sce, file = obj.path)

  # Names
  label_vector <- c(
    "Total_Sparsity" = totSparsity,
    "True_Sparsity" = trueSparsity,
    "Simulated_Sparsity" = simulatedSparsity,
    "Filename" = paste0("zi_", totSparsity, ".RData")
  )

  # Compute UMAP Dimensions
  sob <- CreateSeuratObject(
    counts = sim.sce@assays@data@listData$counts,
    meta.data = as.data.frame(sim.sce@colData)
  )
  sob <- NormalizeData(sob,
    normalization.method = "LogNormalize",
    scale.factor = 10000, verbose = F
  )
  sob <- FindVariableFeatures(sob,
    selection.method = "vst", nfeatures = 2000,
    verbose = F
  )
  sob <- ScaleData(sob, verbose = F)
  sob <- RunPCA(sob, features = VariableFeatures(object = sob), verbose = F)
  sob <- RunUMAP(sob, dims = 1:10, verbose = F)

  # Create Plotting frame for PHATE
  plt.data <- data.frame(
    UMAP_1 = sob@reductions[["umap"]]@cell.embeddings[, 1],
    UMAP_2 = sob@reductions[["umap"]]@cell.embeddings[, 2],
    Simulated_Steps = sim.sce@colData$Step,
    Path = sim.sce@colData$Group
  )

  # Plot PHATE dimensions
  plt <- ggplot(plt.data) +
    geom_point(
      aes(
        x = UMAP_1,
        y = UMAP_2,
        color = Simulated_Steps,
        shape = Path
      ),
      size = 1.5
    ) +
    theme_minimal(base_size = 12) +
    scale_color_viridis(option = "C") +
    ggtitle(
      paste("Total Sparsity:", totSparsity)
    )

  # Return
  return(list(
    parameters = label_vector,
    plots = plt
  ))
}, mc.cores = 1)
# Set names
names(parameter.list) <- names(zi)

# Extract Parameters
parameters <- lapply(parameter.list, function(i) {
  return(i[["parameters"]])
})

# Convert to dataframe
parameter.frame <- do.call("rbind", parameters)

# Save in text files
write.table(parameter.frame,
  file = "/supp_data/Tables/01_ZI_Parameter.Table.tsv",
  sep = "\t", quote = F, row.names = F
)

# Extract Plots
plots <- lapply(parameter.list, function(i) {
  return(i[["plots"]])
})

# Save
saveRDS(plots, paste0(imgPath, "01_Zi_60_90.RDS"))
