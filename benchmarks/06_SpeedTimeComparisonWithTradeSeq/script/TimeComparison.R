# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(tidyverse))

# Set paths
dirPath <- "benchmarks/06_SpeedTimeComparisonWithTradeSeq/data/input/"
resPath <- "benchmarks/06_SpeedTimeComparisonWithTradeSeq/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))

# ReadData
load(paste0(dirPath, "sparsity_60.RData"))

# Readuce Dataset
keepGenes <- sample(rownames(rowData(sim.sce)), size = 2500, replace = F)
keepCells <- sample(rownames(colData(sim.sce)), size = 1500, replace = F)

# Extract raw counts
counts <- as.matrix(sim.sce@assays@data@listData$counts)
cell.metadata <- as.data.frame(colData(sim.sce))
gene.metadata <- as.data.frame(rowData(sim.sce))

# Subset the counts
counts.reduced <- counts[keepGenes, keepCells]
cell.metadata.reduced <- cell.metadata[keepCells, ]
gene.metadata.reduced <- gene.metadata[keepGenes, ]

sim.sce <- SingleCellExperiment(list(counts = counts.reduced))
colData(sim.sce) <- DataFrame(cell.metadata.reduced)
rowData(sim.sce) <- DataFrame(gene.metadata.reduced)

# Extract counts
counts <- as.matrix(sim.sce@assays@data@listData$counts)

# Perform Quantile Normalization as per-tradeSeq paper
normCounts <- calcNormCounts(counts, cat = "FQNorm")

# Extract Cell_metadata
cell_metadata <- as.data.frame(colData(sim.sce))

# Extract Gene_metadata
gene_metadata <- as.data.frame(rowData(sim.sce))

# Prepare Input
pseudotime_table <- cell_metadata[, c("Cell", "Step", "Group")]
lineage_table <- cell_metadata[, c("Cell", "Step", "Group")]

# Add Pseudotime Info
pseudotime_table$Pseudotime1 <- pseudotime_table$Step
pseudotime_table$Pseudotime2 <- pseudotime_table$Step
pseudotime_table <- pseudotime_table[, c("Pseudotime1", "Pseudotime2")]

# Hard Assignmnet for Lineage
lineage_table$Lineage1 <- ifelse(lineage_table$Group == "Path1", 1, 0)
lineage_table$Lineage2 <- ifelse(lineage_table$Group == "Path2", 1, 0)
lineage_table <- lineage_table[, c("Lineage1", "Lineage2")]

mbm <- microbenchmark(
  "TradeSeq" = {
    # Fit GAM
    sce.tradeseq <- fitGAM(
      counts = normCounts,
      pseudotime = pseudotime_table,
      cellWeights = lineage_table,
      parallel = T,
      nknots = 4, verbose = FALSE
    )
    gc()

    # One of the test
    patternRes <- patternTest(sce.tradeseq)
    gc()
  },
  "ScMaSigPro" = {
    # Running scMaSigPro
    scmp.obj <- as_scmp(sim.sce, from = "sce")

    gc()

    # Compress
    scmp.obj <- squeeze(
      scmp.ob = scmp.obj,
      time.col = "Step",
      path.col = "Group",
      method = "Sturges",
      drop.fac = 0.6,
      verbose = T,
      cluster.count.by = "sum"
    )
    gc()

    # Make Design
    scmp.obj <- sc.make.design.matrix(scmp.obj,
      degree = 2,
      time.col = "binnedTime",
      path.col = "path"
    )
    gc()

    # Run p-vector
    scmp.obj <- sc.p.vector(
      scmpObj = scmp.obj, verbose = F, min.obs = 10,
      counts = T, theta = 1,
      offset = T, epsilon = 0.00001
    )
    gc()

    # Run-Step-2
    scmp.obj <- sc.T.fit(
      data = scmp.obj, verbose = F,
      step.method = "backward",
      family = scmp.obj@scPVector@family,
      offset = T
    )
    gc()
  },
  times = 2
)

data <- summary(mbm) %>% as.data.frame()

# Define the order of the bars
data$expr <- factor(data$expr, levels = c("TradeSeq", "ScMaSigPro"))

# Create the bar plot
compareBar <- ggplot(data, aes(x = expr, y = mean, fill = expr)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Execution Times",
    x = "Method",
    y = "Time (seconds)"
  ) +
  theme_minimal()

# Save
ggsave(
  plot = compareBar,
  filename = paste0(resPath, "CompareBarTime.png"),
  dpi = 600
)
