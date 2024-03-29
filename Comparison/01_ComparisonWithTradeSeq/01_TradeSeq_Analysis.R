# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2023

set.seed(007)

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))

# Set paths
inPath <- "/supp_data/ComparisonWithTradeSeq/simulated/sce/"
resPath <- "/supp_data/ComparisonWithTradeSeq/output/"
dir.create(resPath, showWarnings = F)
helpScriptsDir <- "R_Scripts/helper_function/"

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))

# ReadData
load(paste0(inPath, "testTradeSeq.RData"))

# Extract raw counts
counts <- as.matrix(sim.sce@assays@data@listData$counts)

# Perform Quantile Normalization as per-tradeSeq paper
normCounts <- FQnorm(counts)

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

# Evaluate K
# Choosing lowest AIC i.e. 5
# icMat <- evaluateK(counts = normCounts,
#                    pseudotime = pseudotime_table,
#                    cellWeights = lineage_table,
#                    k = 3:15,
#                    nGenes = 200, verbose = T)

# Fit GAM
sce.tradeseq <- fitGAM(
  counts = normCounts,
  pseudotime = pseudotime_table,
  cellWeights = lineage_table,
  parallel = F,
  nknots = 5,
  verbose = FALSE
)

# Save Fitted GAM
save(sce.tradeseq, file = paste0(resPath, "fitGam_TS_Results.RData"))

# Run Different Test
patternRes <- patternTest(sce.tradeseq)
# earlyRes <- earlyDETest(sce.tradeseq)
diffEndRes <- diffEndTest(sce.tradeseq)

# Save All Objects as list
additionalTest <- list(
  patternRes = patternRes,
  # earlyRes = earlyRes,
  diffEndRes = diffEndRes
)
save(additionalTest,
  file = paste0(resPath, "TS_AdditionalTest_ZI_60.RData")
)

# Extract Data
patternResCobra <- patternRes[, "pvalue", drop = F]
# earlyResCobra <- earlyRes[, "pvalue", drop = F]
diffEndResCobra <- diffEndRes[, "pvalue", drop = F]
patternResCobra$pvalue <- as.numeric(patternResCobra$pvalue)
# earlyResCobra$pvalue <- as.numeric(earlyResCobra$pvalue)
diffEndResCobra$pvalue <- as.numeric(diffEndResCobra$pvalue)

# Set Column Names
colnames(patternResCobra) <- c("TS_pattern")
# colnames(earlyResCobra) <- c("TS_early")
colnames(diffEndResCobra) <- c("TS_diffEnd")

# Order Data
patternResCobra <- patternResCobra[mixedorder(rownames(patternResCobra)), , drop = F]
# earlyResCobra <- earlyResCobra[mixedorder(rownames(earlyResCobra)), , drop = F]
diffEndResCobra <- diffEndResCobra[mixedorder(rownames(diffEndResCobra)), , drop = F]

# Create DF
TradeSeq_Clean <- cbind(
  patternResCobra,
  # earlyResCobra,
  diffEndResCobra
)

# Save Dataframe
save(TradeSeq_Clean,
  file = paste0(resPath, "TradeSeq_CobraInput.RData")
)
