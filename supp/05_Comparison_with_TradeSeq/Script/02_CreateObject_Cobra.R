# Title: Creating Cobra Input for ScMaSigPro
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(maSigPro))
suppressPackageStartupMessages(library(gtools))

# Set Paths relative to project
inPath <- "supp/05_Comparison_with_TradeSeq/data/output/"
outPath <- "supp/05_Comparison_with_TradeSeq/data/output/"

# Load Tstep
tstep <- readRDS(paste0(inPath, "MaSigPro/tstep/zi_60_mid_0_shape_0.25.RDS"))

# Get sol
sol <- as.data.frame(tstep$sol)

# Select the column with R2 and P-value
sol <- sol[,c(1,2)]

# Reset column names
colnames(sol) <- c("p_value", "rsq")

# Set NA p-value to 1
sol$p_value[is.na(sol$p_value)] <- 1

# Get genes with r2 > 0.6
sol.sel <- sol[sol$rsq >= 0.6, c(1,2), drop = F]

# Load tradeSeq table
ts.cobra <- readRDS(paste0(outPath,"tradSeqResults/TradeSeq_CobraInput_ZI_60.RDS"))

# get the genes Not selected 
undetected <- rownames(ts.cobra)[!(rownames(ts.cobra) %in% rownames(sol.sel))]

# P_value ==1 and Rsq ==0
undetected <- data.frame(row.names = undetected,
                         "p_value" = rep(1, length(undetected)),
                         "rsq" = rep(0, length(undetected)))

# join
sol <- rbind(undetected, sol.sel)

# Reorder
sol <- sol[mixedorder(rownames(sol)), , drop = F]

# Select the column and rename
colnames(sol) <- c("scMSP_0.6", "rsq")

# Join with TradeSeq Data
cobra.dataset <- cbind(ts.cobra, sol[,1, drop = F])

# Save Dataframe
saveRDS(cobra.dataset, paste0(outPath, "CobraInputObject.RDS"))
