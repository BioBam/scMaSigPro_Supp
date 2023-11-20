# Title: Simulate Datasets with skewness
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

# Set Paths relative to project
outDir <- "/supp_data/benchmarks/12_Interpretation/simulated/"
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
        dropout.shape = 0.03
    )
    
    # Simulate Object
    sim.sce <- splatSimulate(params.groups,
                             method = "paths",
                             verbose = F,
                             path.skew = c(0.5, 0.5)
    )
    
    # Proportion of true Sparsity
    trueSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
    simulatedSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
    totSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)
    
    cat(paste("\nTotal:",totSparsity))
    cat(paste("\nsimulatedSparsity:", simulatedSparsity))
    cat(paste("\ntrueSparsity:", trueSparsity))
    
    # Add gene Info
    gene.info <- add_gene_anno(sim.sce = sim.sce)
    gene.info <- gene.info[mixedsort(gene.info$gene_short_name), ]
    
    # Update the SCE Simulated Object
    rowData(sim.sce) <- DataFrame(gene.info)

    # SaveRDS
    obj.path <- paste0(sce_path, paste0("Sparsity.50.RData"))
    save(sim.sce, file = obj.path)
