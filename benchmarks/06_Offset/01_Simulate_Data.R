# Title: Simulate one Dataset with different Normalization
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))

outDir <- "/supp_data/benchmarks/06_Normalization/simulated/"
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

# Create Base parameters/ Same for All groups
params.groups <- newSplatParams(
  batch.rmEffect = TRUE, # No Batch affect
  batchCells = 3000, # Number of Cells
  nGenes = 5000, # Number of Genes
  seed = 2022, # Set seed
  mean.rate = 0.3, mean.shape = 5, lib.scale = 0.2,
  lib.loc = 12, dropout.type = "experiment",
  group.prob = c(0.5, 0.5), path.from = c(0, 0),
  de.prob = 0.3, de.facLoc = 1, path.nonlinearProb = 0,
  path.sigmaFac = 0,
  path.nSteps = c(1500, 1500),
  dropout.mid = 0,
  dropout.shape = 0.03,
  path.skew = c(0.5, 0.5)
)

# Simulate Object
sim.sce <- splatSimulate(params.groups,
  method = "paths",
  verbose = F
)

# Add gene Info
gene.info <- add_gene_anno(sim.sce = sim.sce)
gene.info <- gene.info[mixedsort(gene.info$gene_short_name), ]

# Norm Method
norm_methods <- c(
  "true",
  "raw",
  "offset",
  "cpm",
  "clr",
  "logN",
  "rc",
  "fquant",
  "sct"
)
names(norm_methods) <- norm_methods
img.list <- list()

# Run
for (i in names(norm_methods)) {
  cat(paste("\nRunning for:", i))

  tryCatch(
    expr = {
      if (i == "true") {
        sce.obj <- SingleCellExperiment(list(counts = sim.sce@assays@data@listData$TrueCounts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "raw") {
        sce.obj <- SingleCellExperiment(list(counts = sim.sce@assays@data@listData$counts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "offset") {
        sce.obj <- SingleCellExperiment(list(counts = sim.sce@assays@data@listData$counts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "cpm") {
        tmpCounts <- calcNormCounts(
          rawCounts = sim.sce@assays@data@listData$counts,
          cat = "libSize", size_fac = 1000000
        )
        sce.obj <- SingleCellExperiment(list(counts = tmpCounts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "clr") {
        tmpCounts <- calcNormCounts(
          rawCounts =
            as.matrix(sim.sce@assays@data@listData$counts),
          cat = "CLR", size_fac = 10000
        )
        sce.obj <- SingleCellExperiment(list(counts = tmpCounts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "logN") {
        tmpCounts <- calcNormCounts(
          rawCounts =
            as.matrix(sim.sce@assays@data@listData$counts),
          cat = "logLibSize", size_fac = 10000
        )
        sce.obj <- SingleCellExperiment(list(counts = tmpCounts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "rc") {
        tmpCounts <- calcNormCounts(
          rawCounts =
            as.matrix(sim.sce@assays@data@listData$counts),
          cat = "libSize", size_fac = 10000
        )
        sce.obj <- SingleCellExperiment(list(counts = tmpCounts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "fquant") {
        tmpCounts <- calcNormCounts(
          rawCounts =
            as.matrix(sim.sce@assays@data@listData$counts),
          cat = "FQNorm", size_fac = 10000
        )
        sce.obj <- SingleCellExperiment(list(counts = tmpCounts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      } else if (i == "sct") {
        tmpCounts <- calcNormCounts(
          rawCounts =
            as.matrix(sim.sce@assays@data@listData$counts),
          cat = "sctransform", size_fac = 10000
        )
        sce.obj <- SingleCellExperiment(list(counts = tmpCounts))
        sce.obj <- runPCA(sce.obj, exprs_values = "counts")
        sce.obj@colData <- sim.sce@colData
        pca <- plotPCA(sce.obj, colour_by = "Step") + ggtitle(paste("Transformation on count", i))
        img.list[[i]] <- pca
      }

      # Add Row and Col Data
      rowData(sce.obj) <- DataFrame(gene.info)
      colData(sce.obj) <- colData(sim.sce)

      # Plot Data
      countLDImgPath <- paste0(imgPath, "PCA_topology/", "type_", i, ".png")
      countPCA.plot <- plot_simulations(sce.obj,
        assay_type = "counts",
        plot3d = F, plot2d = T, frame = 1,
        title.2d = paste("Counts:", i, "; Sparsity: 54; Simulated: 50")
      )
      ggsave(filename = countLDImgPath, plot = countPCA.plot, dpi = 600)

      # SaveRDS
      obj.path <- paste0(sce_path, paste0("type_", i, ".RData"))
      save(sce.obj, file = obj.path)
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}

#
# # Plot for supplemnetary Material
combined_pplot <- ggarrange(
  img.list[["true"]],
  img.list[["raw"]],
  img.list[["offset"]],
  img.list[["cpm"]],
  img.list[["clr"]],
  img.list[["logN"]],
  img.list[["rc"]],
  img.list[["fquant"]],
  img.list[["sct"]],
  ncol = 3, nrow = 3,
  labels = c("A.", "B.", "C.", "D.", "E.", "F.", "G.", "H.", "I.")
)
#
ggsave(
  filename = "Figures/SuppData/06_Sim_Norm.png",
  plot = combined_pplot, dpi = 600, height = 8, width = 14
)
