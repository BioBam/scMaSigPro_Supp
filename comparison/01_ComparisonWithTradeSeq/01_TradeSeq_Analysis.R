# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2023

set.seed(007)

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))

# Set paths
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
rdsPath <- paste0(base_string, "comparison/sim/")
outPath <- paste0(base_string, "comparison/out/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Create Directory if does not exist
dir.create(figPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_hd, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_lr, showWarnings = FALSE, recursive = TRUE)
dir.create(tabPath, showWarnings = FALSE, recursive = TRUE)
dir.create(rdsPath, showWarnings = FALSE, recursive = TRUE)
dir.create(outPath, showWarnings = FALSE, recursive = TRUE)

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))

# ReadData
load(paste0(rdsPath, "testTradeSeq.RData"))

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
  parallel = T,
  nknots = 5,
  verbose = FALSE
)

# Save Fitted GAM
save(sce.tradeseq, file = paste0(outPath, "fitGam_TS_Results.RData"))

# Run Different Test
patternRes <- patternTest(sce.tradeseq)
diffEndRes <- diffEndTest(sce.tradeseq)

# Save All Objects as list
additionalTest <- list(
  patternRes = patternRes,
  diffEndRes = diffEndRes
)
save(additionalTest,
  file = paste0(outPath, "TS_AdditionalTest_ZI_60.RData")
)

# Extract Data
patternResCobra <- patternRes[, "pvalue", drop = F]
diffEndResCobra <- diffEndRes[, "pvalue", drop = F]
patternResCobra$pvalue <- as.numeric(patternResCobra$pvalue)
diffEndResCobra$pvalue <- as.numeric(diffEndResCobra$pvalue)

# Set Column Names
colnames(patternResCobra) <- c("TS_pattern")
colnames(diffEndResCobra) <- c("TS_diffEnd")

# Order Data
patternResCobra <- patternResCobra[mixedorder(rownames(patternResCobra)), , drop = F]
diffEndResCobra <- diffEndResCobra[mixedorder(rownames(diffEndResCobra)), , drop = F]

# Create DF
TradeSeq_Clean <- cbind(
  patternResCobra,
  diffEndResCobra
)

# Save Dataframe
save(TradeSeq_Clean,
  file = paste0(outPath, "TradeSeq_CobraInput.RData")
)
