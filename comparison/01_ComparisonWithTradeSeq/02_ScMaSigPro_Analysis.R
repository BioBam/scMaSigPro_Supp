# Title: Creating Cobra Input for ScMaSigPro
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(gtools))

# Set paths
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
rdsPath <- paste0(base_string, "comparison/sim/")
figPath <- paste0(base_string, "figures/")
outPath <- paste0(base_string, "comparison/out/")
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

# ReadData
load(paste0(rdsPath, "testTradeSeq.RData"))

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
load(paste0(outPath, "TradeSeq_CobraInput.RData"))

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
  file = paste0(outPath, "CobraInputObject.RData")
)
