# Title: Simulate Datasets with different lengths of pseudotime in branches
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(SingleCellExperiment))

# Set Paths relative to project
inPath <- "supp/07_Different_Length_of_Pseudotime/Script/"
outPath <- "supp/07_Different_Length_of_Pseudotime/Script/"

# Load Custom Functions
source(paste0(inPath, "helpScripts/help_plot_simulations().R"))
source(paste0(inPath, "helpScripts/help_calc_bin_size().R"))
source(paste0(inPath, "helpScripts/help_add_gene_anno().R"))

# Different combinations of pseudotim length
path_a_length <- seq(100, 2900, 100)
path_b_length <- seq(2900, 100, -100)
time_length <- list()
for (i in 1:length(path_a_length)) {
  sublist <- c(path_a_length[i], path_b_length[i])
  time_length[[i]] <- sublist
}

## Create a list of parameters
for (i in time_length) {
  # Name of the List
  cat(paste("\nSimulating for Time-Series Length:", paste(i, collapse = " and ")))

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
  trueSparsity <- round(coop::sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
  simulatedSparsity <- round(coop::sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
  totSparsity <- round(coop::sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)

  # Add gene Info
  gene.info <- add_gene_anno(sim.sce = sim.sce)
  status <- data.frame(table(gene.info$status))
  status2 <- data.frame(table(gene.info$status2))

  # Plot Base Gene Mean Histogram
  imgPath <- paste0(inPath, "../data/input/imgs/baseMean")
  dir.create(imgPath, recursive = T, showWarnings = F)
  img <- paste0(imgPath, "/pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity, ".png")
  png(img, width = 800, height = 500)
  hist(gene.info$BaseGeneMean,
    xlab = "Base Gene Mean",
    main = paste0("pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity)
  )
  dev.off()

  # Plotting True Trajectory Topology
  imgPath <- paste0(inPath, "../data/input/imgs/TrueTopology")
  dir.create(imgPath, recursive = T, showWarnings = F)
  img <- paste0(imgPath, "/pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity, ".png")
  png(img, width = 500, height = 500)
  plot_simulations(sim.sce,
    assay_type = "TrueCounts",
    plot3d = F, plot2d = T,
    frame = 2,
    title.2d = paste("True:", trueSparsity, ",",
      "Sim:", simulatedSparsity, ",",
      "Tot:", totSparsity,
      sep = " "
    )
  )
  dev.off()

  # Plotting Simulated Trajectory Topology
  imgPath <- paste0(inPath, "../data/input/imgs/SimTopology")
  dir.create(imgPath, recursive = T, showWarnings = F)
  img <- paste0(imgPath, "/pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity, ".png")
  png(img, width = 500, height = 500)
  plot_simulations(sim.sce,
    assay_type = "counts", frame = 2,
    plot3d = F, plot2d = T,
    title.2d = paste("True:", trueSparsity, ",",
      "Sim:", simulatedSparsity, ",",
      "Tot:", totSparsity,
      sep = " "
    )
  )
  dev.off()

  # Save Gene-Level Information
  table_path <- paste0(inPath, "../data/input/GeneInfo")
  dir.create(table_path, recursive = T, showWarnings = F)
  table <- paste0(table_path, "/pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity, ".tsv")
  write.table(gene.info,
    file = table, sep = "\t", col.names = NA,
    row.names = T
  )

  # Update the SCE object
  gene.info <- gene.info[gtools::mixedsort(gene.info$gene_short_name), ]
  gene.info$Gene <- gene.info$gene_short_name
  rowData(sim.sce) <- DataFrame(gene.info)

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
  imgPath <- paste0(inPath, "../data/input/imgs/cellAssociation")
  dir.create(imgPath, recursive = T, showWarnings = F)
  img <- paste0(imgPath, "/pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity, ".png")
  png(img, width = 750, height = 500)
  p <- ggplot(plt.table, aes(x = Num)) +
    geom_histogram(
      binwidth = 0.5, ,
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
      subtitle = paste0("pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity)
    ) +
    facet_wrap(~Group) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    ylab("Number of cells") +
    xlab("Number of cells associations")
  print(p)
  dev.off()

  # SaveRDS
  obj.path <- paste0(inPath, "../data/input/sceObjects")
  dir.create(obj.path, recursive = T, showWarnings = F)
  obj <- paste0(obj.path, "/pathA_", i[1], "_pathB_", i[2], "_zi_", simulatedSparsity, ".RDS")
  saveRDS(sim.sce, file = obj)
}
