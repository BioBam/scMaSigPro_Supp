# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(maSigPro))
library(microbenchmark)
source("old_Scripts/calcNormCounts.R")

# Set Paths relative to project
inPath <- "supp/08_Computation_Time_Comparison_with_TradeSeq/data/input/"
outPath <- "supp/08_Computation_Time_Comparison_with_TradeSeq/data/"

# ReadData
sce.obj <- readRDS(paste0(
  inPath,
  "zi_mid_60_mid_0_shape_0.25.RDS"
))

# Extract normalized counts
counts <- as.matrix(sce.obj@assays@data@listData$counts)

# Perform Quantile Normalization as per-tradeSeq paper
normCounts <- calcNormCounts(counts, cat = "FQNorm")

# Extract Cell_metadata
cell_metadata <- as.data.frame(colData(sce.obj))

# Extract Gene_metadata
gene_metadata <- as.data.frame(rowData(sce.obj))

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

# Evaluate K
# icMat <- evaluateK(counts = normCounts,
#                    pseudotime = pseudotime_table,
#                    cellWeights = lineage_table,
#                    k = 3:6,
#                    nGenes = 200, verbose = T)


input <- readRDS(paste0(inPath, "zi_60_mid_0_shape_0.25.RDS"))
pbCounts <- input$pbCounts
pbMeta <- input$pbMeta

# Some Variables
poly.degree <- 2
min.gene <- 10
theta.val <- 1
ep <- 0.00001


mbm <- microbenchmark(
  "TrdeSeq" = {
    sce.tradeseq <- fitGAM(
      counts = normCounts,
      pseudotime = pseudotime_table,
      cellWeights = lineage_table,
      parallel = T,
      nknots = 4, verbose = FALSE
    )

    patternRes <- patternTest(sce.tradeseq)
    gc()
  },
  "ScMaSigPro" = {
    design <- make.design.matrix(
      edesign = as.data.frame(pbMeta), # bulkMeta is edesign
      degree = 2,
      time.col = 1, repl.col = 2
    )

    gc <- capture.output(p.vector.fit <- p.vector(
      data = pbCounts, design = design, Q = 0.05,
      MT.adjust = "BH", min.obs = 6, counts = TRUE,
      theta = theta.val, epsilon = ep
    ))

    # Runnig T Step
    gc <- capture.output(tstep.fit <- T.fit(p.vector.fit,
      step.method = "backward",
      family = p.vector.fit$family,
      epsilon = ep
    ))
    gc()
  },
  times = 2
)
mbm


library(ggplot2)

# Create a data frame with the execution times
data <- data.frame(
  expr = c("TrdeSeq", "ScMaSigPro"),
  min = c(143.57062, 57.67871),
  lq = c(143.57062, 57.67871),
  mean = c(144.6410, 58.3368),
  median = c(144.6410, 58.3368),
  uq = c(145.71129, 58.99488),
  max = c(145.71129, 58.99488)
)

# Define the order of the bars
data$expr <- factor(data$expr, levels = c("TrdeSeq", "ScMaSigPro"))

# Create the bar plot
ggplot(data, aes(x = expr, y = mean, fill = expr)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Execution Times",
    x = "Method",
    y = "Time (seconds)"
  ) +
  theme_minimal()
