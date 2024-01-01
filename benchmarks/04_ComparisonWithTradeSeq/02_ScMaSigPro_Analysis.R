# Title: Creating Cobra Input for ScMaSigPro
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(gtools))

# Set paths
dirPath <- "/supp_data/benchmarks/04_ComparisonWithTradeSeq/simulated/sce/"
resPath <- "/supp_data/benchmarks/04_ComparisonWithTradeSeq/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load result of 60% inflation
load(paste0(dirPath, "testTradeSeq.RData"))

# Convert
scmp.obj <- as_scmp(sim.sce,
  from = "sce",
  align_pseudotime = F,
  additional_params = list(
    labels_exist = TRUE,
    exist_ptime_col = "Step",
    exist_path_col = "Group"
  ), verbose = F
)

# Compress
scmp.obj <- sc.squeeze(
  scmpObj = scmp.obj,
  bin_method = "Sturges",
  drop_fac = 0.6,
  verbose = F,
  aggregate = "sum",
  split_bins = F,
  prune_bins = F,
  drop_trails = F,
  fill_gaps = F
)
plotBinTile(scmp.obj)

# Make Design
scmp.obj <- sc.set.poly(scmp.obj,
  poly_degree = 2
)

# Run p-vector
scmp.obj <- sc.p.vector(
  scmpObj = scmp.obj, verbose = F,
  min_na = 1,
  parallel = T,
  offset = T,
  log_offset = TRUE,
  max_it = 1000
)

# Run-Step-2
scmp.obj <- sc.t.fit(
  scmpObj = scmp.obj, verbose = T,
  selection_method = "backward",
  parallel = T,
  nvar_correction = F
)

# Get sol
sol <- showSol(scmpObj = scmp.obj, view = F, return = T, includeInflu = T)

# Select the column with R2 and P-value
sol <- sol[, c(1, 2)]

# Reset column names
colnames(sol) <- c("p_value", "rsq")

# Set NA p-value to 1
sol$p_value[is.na(sol$p_value)] <- 1

# Select by pvalue
sol <- sol[sol$p_value <= 0.05, ]

# Get genes with r2 > 0.6
sol.sel <- sol[sol$rsq >= 0.6, c(1, 2), drop = F]

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
