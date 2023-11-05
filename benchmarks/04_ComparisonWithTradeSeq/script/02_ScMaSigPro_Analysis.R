# Title: Creating Cobra Input for ScMaSigPro
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(gtools))

# Set paths
dirPath <- "benchmarks/04_ComparisonWithTradeSeq/data/input/sce/"
resPath <- "benchmarks/04_ComparisonWithTradeSeq/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load result of 60% inflation
#load(paste0(dirPath, "Test_TradeSeq.RData"))

# Convert
scmp.obj <- as_scmp(sim.sce, from = "sce",
                    align_pseudotime = F,
                    additional_params = list(
                        labels_exist = TRUE,
                        existing_pseudotime_colname = "Step",
                        existing_path_colname = "Group"), verbose = F)

# Compress
scmp.obj <- squeeze(
    scmpObject = scmp.obj,
    bin_method = "Sturges",
    drop.fac = 0.7,
    verbose = F,
    cluster_count_by = "sum",
    split_bins = T,
    prune_bins = F,
    drop_trails = T,
    fill_gaps = F
)
sc.plot.bins.tile(scmp.obj)

# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
  poly_degree = 2
)

# Run p-vector
scmp.obj <- sc.p.vector(
  scmpObj = scmp.obj, verbose = T, min.obs = 1,
  counts = T, theta = 1,parallel = T,MT.adjust = "fdr",
  offset = T
)

# Run-Step-2
scmp.obj <- scMaSigPro::sc.T.fit(
  data = scmp.obj, verbose = T,
  step.method = "backward",parallel = T,
  family = scmp.obj@scPVector@family,
  offset = T
)

# Get sol
sol <- showSol(scmpObj = scmp.obj, view = F, return = T)

# Select the column with R2 and P-value
sol <- sol[, c(1, 2)]

# Reset column names
colnames(sol) <- c("p_value", "rsq")

# Set NA p-value to 1
sol$p_value[is.na(sol$p_value)] <- 1

# Get genes with r2 > 0.6
sol.sel <- sol[sol$rsq >= 0.7, c(1, 2), drop = F]

# Load tradeSeq table
load(paste0(resPath, "TradeSeq_CobraInput.RData"))

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
colnames(sol) <- c("scmp_0.6", "rsq")

# Join with TradeSeq Data
cobra.dataset <- cbind(TradeSeq_Clean, sol[, 1, drop = F])

# Save Dataframe
save(cobra.dataset,
  file = paste0(resPath, "CobraInputObject.RData")
)
