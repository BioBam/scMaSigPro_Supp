# Title: Simulate BigData for evaluation
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
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(viridis))

# Set paths
paramEstimates <- readRDS("/supp_data/benchmarks/00_Parameter_Estimation/output/setty_et_al_d1_splatEstimates.RDS")
outDir <- "/supp_data/ComparisonWithTradeSeq/simulated/"
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

# Create Base parameters/ Same for All groups
params.groups <- newSplatParams(
  batch.rmEffect = TRUE, # No Batch affect
  batchCells = 20000, # 50k Cells
  nGenes = 10000, # 5k Genes
  seed = 2022, # Set seed
  mean.rate = paramEstimates@mean.rate,
  mean.shape = paramEstimates@mean.shape,
  lib.scale = paramEstimates@lib.scale,
  lib.loc = paramEstimates@lib.loc,
  bcv.common = paramEstimates@bcv.common,
  bcv.df = paramEstimates@bcv.df,
  dropout.type = "experiment",
  group.prob = c(0.6, 0.4),
  path.from = c(0, 0),
  de.prob = 0.3,
  de.facLoc = 1,
  path.nonlinearProb = 0.3,
  path.sigmaFac = 0.5,
  out.facLoc = paramEstimates@out.facLoc,
  dropout.mid = paramEstimates@dropout.mid,
  out.facScale = paramEstimates@out.facScale,
  out.prob = paramEstimates@out.prob,
  path.skew = c(0.5, 0.5),
  dropout.shape = -0.5,
  path.nSteps = c(5000, 5000)
)

# Simulate Object
sim.sce <- splatSimulate(
  params = params.groups,
  method = "paths",
  verbose = F
)

# Proportion of true Sparsity
trueSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
simulatedSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
totSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)

cat(paste("\nTotal:", totSparsity))
cat(paste("\nsimulatedSparsity:", simulatedSparsity))
cat(paste("\ntrueSparsity:", trueSparsity))

# Add gene Info
gene.info <- add_gene_anno(sim.sce = sim.sce)
gene.info <- gene.info[mixedsort(gene.info$gene_short_name), ]

# Create Bar
bar.df <- gene.info[, c("status", "status2")]
colnames(bar.df) <- c("DE", "Fold_Change")
bar.df <- as.data.frame(table(bar.df[, c("DE", "Fold_Change")]))
bar.df <- bar.df[bar.df$Freq != 0, ]
bar <- ggplot(bar.df, aes(x = DE, y = Freq, fill = Fold_Change)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 3) +
  theme_minimal() +
  ggtitle("Number of genes in the simulated data",
    subtitle = "Category-wise distribution"
  ) +
  labs(x = "DE", y = "Frequency", fill = "Fold Change") +
  theme(legend.position = "bottom")

# Update the SCE Simulated Object
rowData(sim.sce) <- DataFrame(gene.info)

# SaveRDS
obj.path <- paste0(sce_path, paste0("time_20k_cells.RData"))
save(sim.sce, file = obj.path)

# Compute UMAP Dimensions
sob <- CreateSeuratObject(
  counts = sim.sce@assays@data@listData$counts,
  meta.data = as.data.frame(sim.sce@colData)
)