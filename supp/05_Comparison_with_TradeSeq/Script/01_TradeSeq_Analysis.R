# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))
source("old_Scripts/calcNormCounts.R")

# Set Paths relative to project
inPath <- "supp/05_Comparison_with_TradeSeq/data/input/"
outPath <- "supp/05_Comparison_with_TradeSeq/data/output/"

# ReadData
sce.obj <- readRDS(paste0(inPath,
                          "sceObjects/zi_mid_60_mid_0_shape_0.25.RDS"))

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
lineage_table$Lineage1 <- ifelse(lineage_table$Group == "Path1", 1,0)
lineage_table$Lineage2 <- ifelse(lineage_table$Group == "Path2", 1,0)
lineage_table <- lineage_table[, c("Lineage1", "Lineage2")]

# Evaluate K
# icMat <- evaluateK(counts = normCounts, 
#                    pseudotime = pseudotime_table,
#                    cellWeights = lineage_table,
#                    k = 3:6, 
#                    nGenes = 200, verbose = T)

# Fit GAM
sce.tradeseq <- fitGAM(counts = normCounts, 
                       pseudotime = pseudotime_table,
                       cellWeights = lineage_table,
                       parallel = T,
                       nknots = 4, verbose = FALSE)

# Save Fitted GAM
outPath2 <- paste0(outPath, "tradSeqResults")
dir.create(
    outPath2,
    showWarnings = F, recursive = T
)
saveRDS(sce.tradeseq, paste0(outPath2,"/fitGam_ZI_60.RDS"))

# Run Different Test
patternRes <- patternTest(sce.tradeseq)
earlyRes <- earlyDETest(sce.tradeseq)
diffEndRes <- diffEndTest(sce.tradeseq)
saveRDS(list(patternRes = patternRes,
             associationRes = earlyRes,
             diffEndRes = diffEndRes),
        paste0(outPath2, "/AdditionalTest_ZI_60.RDS"))

# Extract Data
patternResCobra <- patternRes[, "pvalue",drop = F]
earlyResCobra <- earlyRes[, "pvalue",drop = F]
diffEndResCobra <- diffEndRes[, "pvalue",drop = F]
patternResCobra$pvalue <- as.numeric(patternResCobra$pvalue)
earlyResCobra$pvalue <- as.numeric(earlyResCobra$pvalue)
diffEndResCobra$pvalue <- as.numeric(diffEndResCobra$pvalue)

# Set Column Names
colnames(patternResCobra) <- c("TS_pattern")
colnames(earlyResCobra)<- c("TS_early")
colnames(diffEndResCobra)<- c("TS_diffEnd")

# Order Data
patternResCobra <- patternResCobra[mixedorder(rownames(patternResCobra)), , drop = F]
earlyResCobra <- earlyResCobra[mixedorder(rownames(earlyResCobra)), , drop = F]
diffEndResCobra <- diffEndResCobra[mixedorder(rownames(diffEndResCobra)), , drop = F]

# Create DF
TradeSeq_Clean <- cbind(patternResCobra, earlyResCobra, diffEndResCobra)

# Save Dataframe
saveRDS(TradeSeq_Clean, paste0(outPath2, "/TradeSeq_CobraInput_ZI_60.RDS"))