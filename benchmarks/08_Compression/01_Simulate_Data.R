# Title: Simulate 9 Datasets with different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr)) # Set Paths relative to project

# Set paths
outDir <- "/supp_data/benchmarks/08_Compression/simulated/"
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
# "Inflation Level" = c(dropout.mid, dropout.shape)
zi <- list(
  "sparsity_50" = c(0, 0.03)
)

# List of Images
img.list <- list()

## Create a list of parameters
for (i in names(zi)) {
  # stop("Expected Stop")

  # Get Sparsity Level
  sparsityLevel <- str_remove(pattern = "sparsity_", i)

  # Name of the List
  cat(paste("\nSimulating for Zero-Inflation level:", sparsityLevel))

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
    path.skew = c(0.5, 0.5),
    path.nSteps = c(1500, 1500)
  )

  # Simulate Object
  sim.sce <- splatSimulate(params.groups,
    method = "paths",
    verbose = F,
    dropout.mid = as.numeric(zi[[i]][1]),
    dropout.shape = as.numeric(zi[[i]][2])
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

  # Update the SCE Simulated Object
  rowData(sim.sce) <- DataFrame(gene.info)

  # Plot Base Gene Mean Histogram
  histImgName <- paste0(imgPath, "base_expression_gene_hist/", "zi_", sparsityLevel, ".png")

  # Create Histogram
  base.exp.hist <- ggplot(gene.info, aes(x = BaseGeneMean)) +
    geom_histogram(fill = "blue", color = "black", alpha = 0.6) +
    labs(
      title = paste0("Sparsity: ", sparsityLevel),
      subtitle = paste("Sparsity:", totSparsity, "Simulated:", simulatedSparsity),
      x = "Base Gene Mean",
      y = "Count"
    ) +
    theme_minimal()

  ggsave(filename = histImgName, plot = base.exp.hist, dpi = 600)

  # Plotting True Trajectory Topology
  truTopImgName <- paste0(imgPath, "true_topology_pca_step/", "zi_", sparsityLevel, ".png")
  truTopImg.plot <- plot_simulations(sim.sce,
    assay_type = "TrueCounts",
    plot3d = F, plot2d = T, frame = 2,
    title.2d = paste("Sparsity:", totSparsity, "Simulated:", simulatedSparsity)
  )
  ggsave(filename = truTopImgName, plot = truTopImg.plot, dpi = 600)

  # Plot Simulated Topology
  simTopImgName <- paste0(imgPath, "sim_topology_pca_step/", "zi_", sparsityLevel, ".png")
  simTopImg.plot <- plot_simulations(sim.sce,
    assay_type = "counts",
    plot3d = F, plot2d = T, frame = 2,
    title.2d = paste("Sparsity:", totSparsity, "Simulated:", simulatedSparsity)
  )
  ggsave(filename = simTopImgName, plot = simTopImg.plot, dpi = 600)

  img.list[[sparsityLevel]] <- simTopImg.plot

  # Plotting True Trajectory Topology Group
  truTopImgNameGroup <- paste0(imgPath, "true_topology_pca_group/", "zi_", sparsityLevel, ".png")
  truTopImgGroup.plot <- plot_simulations(sim.sce,
    assay_type = "TrueCounts",
    plot3d = F, plot2d = T, frame = 1,
    title.2d = paste("Sparsity:", totSparsity, "Simulated:", simulatedSparsity)
  )
  ggsave(filename = truTopImgNameGroup, plot = truTopImgGroup.plot, dpi = 600)

  # Plot Simulated Topology
  simTopImgNameGroup <- paste0(imgPath, "sim_topology_pca_group/", "zi_", sparsityLevel, ".png")
  simTopImgGroup.plot <- plot_simulations(sim.sce,
    assay_type = "counts",
    plot3d = F, plot2d = T, frame = 1,
    title.2d = paste("Sparsity:", totSparsity, "Simulated:", simulatedSparsity)
  )
  ggsave(filename = simTopImgNameGroup, plot = simTopImgGroup.plot, dpi = 600)

  # Extract Cell Metadata Information
  cell.meta <- as.data.frame(colData(sim.sce))

  # Select Columns
  plt.table <- cell.meta[, c("Cell", "Step", "Group")]

  # Group by
  plt.table <- plt.table %>%
    group_by(Step, Group) %>%
    summarise(cluster.members = paste0(Cell, collapse = "|"))
  plt.table$Num <- apply(plt.table, 1, calc_bin_size)

  # Select Columns
  plt.table <- plt.table[, !(colnames(plt.table) %in% "cluster.members")]

  # Plot Cell Association

  cellAssociation <- paste0(imgPath, "cellAssociation/", "zi_", sparsityLevel, ".png")
  p <- ggplot(plt.table, aes(x = Num)) +
    geom_histogram(
      binwidth = 0.5,
      color = "#f68a53", fill = "#f68a53", alpha = 0.5
    ) +
    geom_vline(aes(xintercept = mean(Num)), linetype = "dashed", color = "#139289") +
    theme_classic() +
    theme(
      legend.position = "none", strip.text = element_text(size = rel(2)),
      axis.text = element_text(size = rel(1)),
      panel.grid.major = element_line(linewidth = 0.7, linetype = "dotted"),
      panel.grid.minor = element_line(linewidth = 0.2)
    ) +
    ggtitle("Distribution of cells per Time-point",
      subtitle = paste("Sparsity:", totSparsity, "Simulated:", simulatedSparsity)
    ) +
    facet_wrap(~Group) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    ylab("Number of cells") +
    xlab("Number of cells associations")

  ggsave(filename = cellAssociation, plot = p, dpi = 600, height = 5, width = 6)

  # SaveRDS
  obj.path <- paste0(sce_path, paste0("sparsity_", sparsityLevel, ".RData"))
  save(sim.sce, file = obj.path)
}
