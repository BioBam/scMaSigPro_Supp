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

# Set paths
paramEstimates <- readRDS("/supp_data/benchmarks/00_Parameter_Estimation/output/setty_et_al_d1_splatEstimates.RDS")
outDir <- "/supp_data/benchmarks/03_Different_Length/simulated/"
helpScriptsDir <- "R_Scripts/helper_function/"
imgPath <- paste0(outDir, "png/")
sce_path <- paste0(outDir, "sce/")
tab_path <- paste0("/supp_data/Tables/")

# Create Directories
dir.create(outDir, showWarnings = F, recursive = T)
dir.create(imgPath, recursive = T, showWarnings = F)
dir.create(sce_path, recursive = T, showWarnings = F)
dir.create(tab_path, recursive = T, showWarnings = F)

# Load Custom Functions
source(paste0(helpScriptsDir, "plot_simulations().R"))
source(paste0(helpScriptsDir, "add_gene_anno().R"))
source(paste0(helpScriptsDir, "calc_bin_size.R"))

# Different length
len <- list(
  "len_100_2900" = c(100, 2900),
  "len_500_2500" = c(500, 2500),
  "len_1000_1000" = c(1000, 2000),
  "len_1500_1500" = c(1500, 1500)
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
  dropout.shape = -0.5
)

# Generate Datasets
parameter.list <- mclapply(names(len), function(length, params_groups = params.groups,
                                                sce.path = sce_path) {
  # Get Variables
  length_value <- len[[length]]

  # Simulate Object
  sim.sce <- splatSimulate(
    params = params_groups,
    method = "paths",
    verbose = F,
    group.prob = c(1500 / length_value[1], 1500 / length_value[2]), # "Override"
    path.nSteps = c(length_value[1], length_value[2])
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
  obj.path <- paste0(sce.path, paste0("len_", paste(length_value, collapse = "."), ".RData"))
  save(sim.sce, file = obj.path)

  # Names
  label_vector <- c(
    "Total_Sparsity" = totSparsity,
    "True_Sparsity" = trueSparsity,
    "Simulated_Sparsity" = simulatedSparsity,
    "Filename" = paste0("len_", paste(length_value, collapse = "."), ".RData")
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
      paste("Skew:", paste(length_value, collapse = "."))
    )

  # Return
  return(list(
    parameters = label_vector,
    plots = plt
  ))
}, mc.cores = 20)
# Set names
names(parameter.list) <- names(len)

# Extract Parameters
parameters <- lapply(parameter.list, function(i) {
  return(i[["parameters"]])
})

# Convert to dataframe
parameter.frame <- do.call("rbind", parameters)

# Save in text files
write.table(parameter.frame,
  file = paste0(tab_path, "03_len_Parameter.Table.tsv"),
  sep = "\t", quote = F, row.names = F
)

# Extract Plots
plots <- lapply(parameter.list, function(i) {
  return(i[["plots"]])
})

# Save
saveRDS(plots, paste0(imgPath, "03_len_100_1500.RDS"))
