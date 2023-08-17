# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))

# Load Prefix
parentFix <- "benchMarking_Normalization/" 
scriptFix <- paste0("old_Scripts/")
inPrefix <- paste0(parentFix, "data/input/")
outPrefix <- paste0(parentFix, "data/output/")

# Load Custim Functions
load <- function(){
    source(paste0(scriptFix, "plot_simulations().R"))
    source(paste0(scriptFix, "add_gene_anno().R"))
    source(paste0(scriptFix, "calc_bin_size.R"))
}

# These will be used for skewness and will model start and ending number of cells
skew <- seq(0, 1, 0.1) # 0.5

# Length of the trajectory
length <- seq(100, 2500, 200) # 1500

# Outlier Probabbility
zi <- round(seq(-1.2, 1, 0.2), 1) # -0.4

# Parameters
run_param <- list(skew = skew, length = length, zi = zi)

load()

## Create a list of parameters
for (i in names(run_param)) {

  # Name of the List
  print(paste("Making Data for", i))

  if (i == "skew") {
    # FileName
    h1 <- "skew"
    h2 <- "density"

    # Create Base parameters/ Same for All groups
    params.groups <- newSplatParams(
      batch.rmEffect = TRUE, # No Batch affect
      batchCells = 3000, # Number of Cells
      nGenes = 5000, # Number of Genes
      seed = 2022, # Set seed
      mean.rate = 0.3, mean.shape = 5, lib.scale = 0.2,
      lib.loc = 12, dropout.mid = 0,
      dropout.type = "experiment",
      group.prob = c(0.5, 0.5), path.from = c(0, 0),
      de.prob = 0.3, de.facLoc = 1, path.nonlinearProb = 0,
      path.sigmaFac = 0,

      # Changing Parameters
      dropout.shape = -0.4,
      path.nSteps = c(1500, 1500)
    )
    
  } else if (i == "length") {
    # FileName
    h1 <- "length"
    h2 <- "strEnd"

    # Create Base parameters/ Same for All groups
    params.groups <- newSplatParams(
      batch.rmEffect = TRUE, # No Batch affect
      batchCells = 3000, # Number of Cells
      nGenes = 5000, # Number of Genes
      seed = 2022, # Set seed
      mean.rate = 0.3, mean.shape = 5, lib.scale = 0.2,
      lib.loc = 12, dropout.mid = 0,
      dropout.type = "experiment",
      group.prob = c(0.5, 0.5), path.from = c(0, 0),
      de.prob = 0.3, de.facLoc = 1, path.nonlinearProb = 0,
      path.sigmaFac = 0,

      # Changing Parameters
      dropout.shape = -0.4,
      path.skew = c(0.5, 0.5)
    )
  } else if (i == "zi") {
    # FileName
    h1 <- "zi"
    h2 <- "shape"

    # Create Base parameters/ Same for All groups
    params.groups <- newSplatParams(
      batch.rmEffect = TRUE, # No Batch affect
      batchCells = 3000, # Number of Cells
      nGenes = 5000, # Number of Genes
      seed = 2022, # Set seed
      mean.rate = 0.3, mean.shape = 5, lib.scale = 0.2,
      lib.loc = 12, dropout.mid = 0,
      dropout.type = "experiment",
      group.prob = c(0.5, 0.5), path.from = c(0, 0),
      de.prob = 0.3, de.facLoc = 1, path.nonlinearProb = 0,
      path.sigmaFac = 0,

      # Changing Parameters
      path.skew = c(0.5, 0.5),
      path.nSteps = c(1500, 1500)
    )
  }

  # Run for individual values
  for (j in run_param[[i]]) {
      
      print(paste("Making Data for", i, "with ", j))
      
    if (i == "zi") {
      sim.sce <- splatSimulate(params.groups,
        method = "paths",
        verbose = F,
        dropout.shape = j
      )
    } else if (i == "skew") {
      sim.sce <- splatSimulate(params.groups,
        method = "paths",
        path.skew = c(j, j),
        verbose = F
      )
    } else if (i == "length") {
      sim.sce <- splatSimulate(params.groups,
        method = "paths",
        verbose = F,
        path.nSteps = c(j, j)
      )
    }

    # Proportion of true Sparsity
    trueSparsity <- round(coop::sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
    simulatedSparsity <- round(coop::sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
    totSparsity <- round(coop::sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)

    # Add gene Info
    gene.info <- add_gene_anno(sim.sce = sim.sce)

    # Plotting Simulations
    dir_path <- paste0(outPrefix, "baseMean/")
    dir.create(dir_path, showWarnings = F)
    image_path <- paste0(dir_path, h1, "_", h2, "_", j, ".png")
    png(image_path, width = 800, height = 500)
    hist(gene.info$BaseGeneMean)
    dev.off()

    # Plotting Simulations
    dir_path <- paste0(outPrefix, "TrueTrends/")
    dir.create(dir_path, showWarnings = F)
    image_path <- paste0(dir_path, h1, "_", h2, "_", j, ".png")
    png(image_path, width = 800, height = 500)
    plot_simulations(sim.sce,
      assay_type = "TrueCounts",
      plot3d = F, plot2d = T,
      title.2d = paste("True:", trueSparsity, ",",
        "Sim:", simulatedSparsity, ",",
        "Tot:", totSparsity,
        sep = " "
      )
    )
    dev.off()
    
    # Simulated Counts
    dir_path <- paste0(outPrefix, "SimTrends/")
    dir.create(dir_path, showWarnings = F)
    image_path <- paste0(dir_path, h1, "_", h2, "_", j, "_.png")
    png(image_path, width = 800, height = 500)
    plot_simulations(sim.sce,
      assay_type = "counts",
      plot3d = F, plot2d = T,
      title.2d = paste("True:", trueSparsity, ",",
        "Sim:", simulatedSparsity, ",",
        "Tot:", totSparsity,
        sep = " "
      )
    )
    dev.off()

    # Save annoated gene info
    dir_path <- paste0(outPrefix, "GeneTables/")
    dir.create(dir_path, showWarnings = F)
    table_path <- paste0(dir_path, h1, "_", h2, "_", j, ".tsv")
    write.table(gene.info, file = table_path, sep = "\t",
                col.names = NA, row.names = T
                )
    
    # Create MaSigPro Table
    gene.info <- gene.info[gtools::mixedsort(gene.info$gene_short_name), ]
    gene.info$Gene <- gene.info$gene_short_name
    rowData(sim.sce) <- DataFrame(gene.info)
    cell.meta <- as.data.frame(colData(sim.sce))

    # Select Columns
    plt.table <- cell.meta[, c("Cell", "Step", "Group")]

    # Group by
    plt.table <- plt.table %>%
      group_by(Step, Group) %>%
      summarise(cluster.members = paste0(Cell, collapse = "|"))
    plt.table$Num <- apply(plt.table, 1, calc_bin_size)

    # Select Columns
    plt.table <- plt.table[, !colnames(plt.table) %in% "cluster.members"]

    # Plot of cell distribution
    dir_path <- paste0(outPrefix, "CellDist/")
    dir.create(dir_path, showWarnings = F)
    image_path2 <- paste0(dir_path, h1, "_", h2, "_", j, ".png")
    png(image_path2, width = 500, height = 500)
    # Plot Data
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
      ggtitle("Distribution of cells per Time-point") +
      facet_wrap(~Group) +
      scale_x_continuous(breaks = seq(0, 30, by = 2)) +
      ylab("Number of cells") +
      xlab("Number of cells associations")

    print(p)
    dev.off()

    # SaveRDS
    dir_path <- paste0(outPrefix, "SimObjectsSCE/")
    dir.create(dir_path, showWarnings = F)
    obj.path <- paste0(dir_path, h1, "_", h2, "_", j, ".RDS")
    saveRDS(sim.sce, file = obj.path)
  }
}
