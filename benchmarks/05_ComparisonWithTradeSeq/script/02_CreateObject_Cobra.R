# Title: Creating Cobra Input for ScMaSigPro
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(gtools))

# Set paths
dirPath <- "benchmarks/05_ComparisonWithTradeSeq/data/input/"
resPath <- "benchmarks/05_ComparisonWithTradeSeq/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load result of 60% inflation
load(paste0(dirPath, "scmp.obj.sparsity.60.RData"))

# Get sol
sol <- showSol(scmpObj = scmp.obj, view = F, return = T)

# Select the column with R2 and P-value
sol <- sol[, c(1, 2)]

# Reset column names
colnames(sol) <- c("p_value", "rsq")

# Set NA p-value to 1
sol$p_value[is.na(sol$p_value)] <- 1

# Get genes with r2 > 0.6
sol.sel <- sol[sol$rsq >= 0.6, c(1, 2), drop = F]

# Load tradeSeq table
load(paste0(resPath, "TradeSeq_CobraInput_ZI_60.RData"))

# get the genes Not selected
undetected <- rownames(TradeSeq_Clean)[!(rownames(TradeSeq_Clean) %in% rownames(sol.sel))]

# P_value ==1 and Rsq ==0
undetected <- data.frame(
  row.names = undetected,
  "p_value" = rep(1, length(undetected)),
  "rsq" = rep(0, length(undetected))
)

# join
sol <- rbind(undetected, sol.sel)

# Reorder
sol <- sol[mixedorder(rownames(sol)), , drop = F]

# Select the column and rename
colnames(sol) <- c("scMSP_0.6", "rsq")

# Join with TradeSeq Data
cobra.dataset <- cbind(TradeSeq_Clean, sol[, 1, drop = F])

# Save Dataframe
save(cobra.dataset,
  file = paste0(resPath, "CobraInputObject.RData")
)
