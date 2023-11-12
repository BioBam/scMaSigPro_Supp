# Title: Simulate Datasets with different lenghts
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))

# Set Paths relative to project
dirPath <- "benchmarks/03_Different_Length/data/simulated/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Create Path
imgPath <- paste0(dirPath, "png/")
sce_path <- paste0(dirPath, "sce/")
dir.create(imgPath, recursive = T, showWarnings = F)
dir.create(sce_path, recursive = T, showWarnings = F)

# Load Custom Functions
source(paste0(helpScriptsDir, "plot_simulations().R"))
source(paste0(helpScriptsDir, "add_gene_anno().R"))
source(paste0(helpScriptsDir, "calc_bin_size.R"))

# Unequal Arms
path_a_length <- seq(100, 2900, 100)
path_b_length <- seq(2900, 100, -100)
time_length <- list()
for (i in 1:length(path_a_length)) {
  sublist <- c(path_a_length[i], path_b_length[i])
  time_length[[i]] <- sublist
}
names(time_length) <- lapply(time_length, FUN = function(i){
    return(paste(i, collapse = "_"))
})
time_length <- time_length[seq(2,15,2)]

# List of Images
img.list <- list()

## Create a list of parameters
for (i in time_length) {
  # Get Sparsity Level
  ArmLength <- paste(i, collapse = "_and_")

  # Name of the List
  cat(paste("\nSimulating for Time-Series Length:", ArmLength))
  
  # Adjust group probablity according to path length
  path1_length <- str_split_i(ArmLength, pattern = "_", 1)
  path2_length <- str_split_i(ArmLength, pattern = "_", 3)
  
  if(path2_length > path1_length){
      path1.prob <- 0.3
      path2.prob <- 0.7
  }else{
      path2.prob <- 0.3
      path1.prob <- 0.7
  }

  # Create Base parameters/ Same for All groups
  params.groups <- newSplatParams(
    batch.rmEffect = TRUE, # No Batch affect
    batchCells = 3000, # Number of Cells
    nGenes = 5000, # Number of Genes
    seed = 2022, # Set seed
    mean.rate = 0.3, mean.shape = 5, lib.scale = 0.2,
    lib.loc = 12, dropout.type = "experiment",
    group.prob = c(path1.prob, path2.prob), path.from = c(0, 0),
    de.prob = 0.3, de.facLoc = 1, path.nonlinearProb = 0,
    path.sigmaFac = 0,
    dropout.mid = 1.7, dropout.shape = -1,
    path.skew = c(0.5, 0.5)
  )

  # Simulate Object
  sim.sce <- splatSimulate(params.groups,
    method = "paths",
    verbose = F,
    path.nSteps = c(i[1], i[2])
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

  # Plot Base Gene Mean Histogram
  histImgName <- paste0(imgPath, "base_expression_gene_hist/", "Arm_", ArmLength, ".png")

  # Create Histogram
  base.exp.hist <- ggplot(gene.info, aes(x = BaseGeneMean)) +
    geom_histogram(fill = "blue", color = "black", alpha = 0.6) +
    labs(
      title = paste0("Path Length: ", ArmLength),
      subtitle = paste("Total Sparsity of Dataset:", totSparsity, "Paths: Equal"),
      x = "Base Gene Mean",
      y = "Count"
    ) +
    theme_minimal()

  ggsave(filename = histImgName, plot = base.exp.hist, dpi = 600)

  # Plotting True Trajectory Topology
  truTopImgName <- paste0(imgPath, "true_topology_pca_step/", "Arm_", ArmLength, ".png")
  truTopImg.plot <- plot_simulations(sim.sce,
    assay_type = "TrueCounts",
    plot3d = F, plot2d = T, frame = 2,
    title.2d = paste("PathLength:", ArmLength, totSparsity, "Simulated:", simulatedSparsity)
  )
  ggsave(filename = truTopImgName, plot = truTopImg.plot, dpi = 600)
  
  img.list[[ArmLength]] <- truTopImg.plot

  # Plot Simulated Topology
  simTopImgName <- paste0(imgPath, "sim_topology_pca_step/", "Arm_", ArmLength, ".png")
  simTopImg.plot <- plot_simulations(sim.sce,
    assay_type = "counts",
    plot3d = F, plot2d = T, frame = 2,
    title.2d = paste("PathLength:", totSparsity, "Simulated:", simulatedSparsity)
  )
  ggsave(filename = simTopImgName, plot = simTopImg.plot, dpi = 600)

  # Plotting True Trajectory Topology Group
  truTopImgNameGroup <- paste0(imgPath, "true_topology_pca_group/", "Arm_", ArmLength, ".png")
  truTopImgGroup.plot <- plot_simulations(sim.sce,
    assay_type = "TrueCounts",
    plot3d = F, plot2d = T, frame = 1,
    title.2d = paste("PathLength:", totSparsity, "Simulated:", simulatedSparsity)
  )
  ggsave(filename = truTopImgNameGroup, plot = truTopImgGroup.plot, dpi = 600)

  # Plot Simulated Topology
  simTopImgNameGroup <- paste0(imgPath, "sim_topology_pca_group/", "Arm_", ArmLength, ".png")
  simTopImgGroup.plot <- plot_simulations(sim.sce,
    assay_type = "counts",
    plot3d = F, plot2d = T, frame = 1,
    title.2d = paste("PathLength:", totSparsity, "Simulated:", simulatedSparsity)
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
  cellAssociation <- paste0(imgPath, "cellAssociation/", "Arm_", ArmLength, ".png")
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
      subtitle = paste("PathLength:", totSparsity, "Simulated:", simulatedSparsity)
    ) +
    facet_wrap(~Group) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    ylab("Number of cells") +
    xlab("Number of cells associations")

  ggsave(filename = cellAssociation, plot = p, dpi = 600, height = 5, width = 6)

  # SaveRDS
  obj.path <- paste0(sce_path, paste0("Arm_", ArmLength, ".RData"))
  save(sim.sce, file = obj.path)
}

combined_pplot <- ggarrange(img.list[["300_and_2700"]],
                            img.list[["600_and_2400"]],
                            img.list[["900_and_2100"]],
                            img.list[["1200_and_1800"]],
                            img.list[["1500_and_1500"]],
                            ncol = 3, nrow = 3,
                            labels = c("A.","B.","C.","D.","E."))

ggsave(filename = "Images/Supp_Fig_4_Vary_Branch_Length.png", plot = combined_pplot, dpi = 600, height = 8, width = 14)
