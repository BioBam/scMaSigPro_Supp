# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2024

trash <- gc()
trash <- NULL

# Create Pid
pid_ts_8cpu <- 101
pid_scmp_8cpu <- 101
pid_ts_1cpu <- 101
pid_scmp_1cpu <- 101

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(viridis))

# Set paths
base_string <- "../../../scMaSigPro_supp_data/"
base_string_2 <- "../../"
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Load Base data
paramEstimates <- readRDS(paste0(base_string, "benchmarks/00_Parameter_Estimation/output/setty_et_al_d1_splatEstimates.RDS"))

# Create Directory if does not exist
dir.create(figPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_hd, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_lr, showWarnings = FALSE, recursive = TRUE)
dir.create(tabPath, showWarnings = FALSE, recursive = TRUE)

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))
source(paste0(helpScriptsDir, "add_gene_anno().R"))
source(paste0(helpScriptsDir, "calc_bin_size.R"))

cat(paste("\nInitiate Simulation:", timestamp(suffix = "", prefix = "", quiet = TRUE)))

# Create Base parameters/ Same for All groups
params.groups <- newSplatParams(
  batch.rmEffect = TRUE, # No Batch affect
  batchCells = 6000, # Number of Cells
  nGenes = 2000, # Number of Genes
  seed = 2022, # Set seed
  mean.rate = paramEstimates@mean.rate,
  mean.shape = paramEstimates@mean.shape,
  lib.scale = paramEstimates@lib.scale,
  lib.loc = paramEstimates@lib.loc,
  bcv.common = paramEstimates@bcv.common,
  bcv.df = paramEstimates@bcv.df,
  dropout.type = "experiment",
  group.prob = c(0.6, 0.4),
  path.from = c(0, 0),
  de.prob = 0.3,
  de.facLoc = 1,
  path.nonlinearProb = 0.3,
  path.sigmaFac = 0.5,
  out.facLoc = paramEstimates@out.facLoc,
  dropout.mid = paramEstimates@dropout.mid,
  out.facScale = paramEstimates@out.facScale,
  out.prob = paramEstimates@out.prob,
  path.skew = c(0.4, 0.6),
  dropout.shape = -0.5,
  path.nSteps = c(3000, 3000)
)

# Simulate Object
sim.sce <- splatSimulate(
  params = params.groups,
  method = "paths",
  verbose = F
)

cat(paste("\nSimulation completed:", timestamp(suffix = "", prefix = "", quiet = TRUE)))

trash <- gc()
trash <- NULL

# Add gene Info
gene.info <- add_gene_anno(sim.sce = sim.sce)
gene.info <- gene.info[mixedsort(gene.info$gene_short_name), ]

# Update the SCE Simulated Object
rowData(sim.sce) <- DataFrame(gene.info)

# Extract counts
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

cat(paste("\nClean-up for tradeSeq completed:", timestamp(suffix = "", prefix = "", quiet = TRUE)))

# Running scMaSigPro
scmp.obj <- as_scmp(sim.sce,
  from = "sce",
  align_pseudotime = T,
  additional_params = list(
    labels_exist = TRUE,
    exist_ptime_col = "Step",
    exist_path_col = "Group"
  ), verbose = F
)

# Squeeze
scmp.obj <- sc.squeeze(
  scmpObj = scmp.obj,
  drop_fac = 40,
  verbose = F,
  aggregate = "sum",
  split_bins = F,
  prune_bins = F,
  drop_trails = F,
  fill_gaps = F
)

# Make Design
scmp.obj <- sc.set.poly(scmp.obj,
  poly_degree = 3
)

# Cleanup
sim.sce <- NULL
counts <- NULL
rm(counts)
rm(sim.sce)
trash <- gc()
trash <- NULL

cat(paste("\nClean-up for scMaSigPro completed:", timestamp(suffix = "", prefix = "", quiet = TRUE)))

mbm <- microbenchmark(
  "TradeSeq_1_CPU" = {
    trash <- gc()
    trash <- NULL
    sce.tradeseq <- fitGAM(counts = normCounts, pseudotime = pseudotime_table, cellWeights = lineage_table, parallel = F, nknots = 3, verbose = FALSE)
    patternRes <- patternTest(sce.tradeseq)
    cat(paste("\nFinished TradeSeq with 1 CPU(", pid_ts_1cpu, "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
    pid_ts_1cpu <- pid_ts_1cpu + 1
    trash <- gc()
    trash <- NULL
  },
  "TradeSeq_8_CPU" = {
    trash <- gc()
    trash <- NULL
    sce.tradeseq <- fitGAM(counts = normCounts, pseudotime = pseudotime_table, cellWeights = lineage_table, parallel = T, nknots = 3, verbose = FALSE)
    patternRes <- patternTest(sce.tradeseq)
    cat(paste("\nFinished TradeSeq with 8 CPU(", pid_ts_8cpu, "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
    pid_ts_8cpu <- pid_ts_8cpu + 1
    trash <- gc()
    trash <- NULL
  },
  "ScMaSigPro_1_CPU" = {
    trash <- gc()
    trash <- NULL
    scmp.obj <- sc.p.vector(scmpObj = scmp.obj, verbose = FALSE, min_na = 1, parallel = F, offset = T, max_it = 1000)
    scmp.obj <- sc.t.fit(scmpObj = scmp.obj, verbose = FALSE, selection_method = "backward", parallel = F, offset = T)
    cat(paste("\nFinished ScMaSigPro with 1 CPU(", pid_scmp_1cpu, "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
    pid_scmp_1cpu <- pid_scmp_1cpu + 1
    trash <- gc()
    trash <- NULL
  },
  "ScMaSigPro_8_CPU" = {
    trash <- gc()
    trash <- NULL
    cat(paste("\nComputing ScMaSigPro with 8 CPU(", pid_scmp_8cpu, "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
    scmp.obj <- sc.p.vector(scmpObj = scmp.obj, verbose = FALSE, min_na = 1, parallel = T, n_cores = 8, offset = T, max_it = 1000)
    scmp.obj <- sc.t.fit(scmpObj = scmp.obj, verbose = FALSE, selection_method = "backward", parallel = TRUE, n_cores = 8, offset = T)
    cat(paste("\nFinished ScMaSigPro with 8 CPU(", pid_scmp_8cpu, "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
    pid_scmp_8cpu <- pid_scmp_8cpu + 1
    trash <- gc()
    trash <- NULL
  },
  times = 10
)

cat(paste("\nBenchmark completed:", timestamp(suffix = "", prefix = "", quiet = TRUE)))

# Extract Data for time
data <- mbm %>% as.data.frame()

# Convert to minutes
data$time <- data$time / 1e9 / 60

# Add function
data$expr <- as.factor(data$expr)

# Save the plot
write.table(
  x = data,
  file = paste0(tabPath, "3k_Cells_TS_Time_Profiling.txt"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

# Creating a violin plot for time
compare_violin <- ggplot(data, aes(x = expr, y = time, fill = expr, color = expr)) +
  geom_violin(trim = TRUE) +
  geom_jitter(width = 0.3, height = 0.3) +
  labs(title = "Processing Times", subtitle = "Evaluated 10 times", x = "", y = "Time (minutes)") +
  theme_minimal(base_size = 12) +
  scale_fill_viridis(discrete = TRUE) +
  coord_flip() +
  scale_color_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "top") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 0.5))

# Save the plot
ggsave(
  plot = compare_violin,
  filename = paste0(figPath_hd, "04_B_3k_ts_time.png"),
  dpi = 600, width = 10, height = 10
)
ggsave(
  plot = compare_violin,
  filename = paste0(figPath_lr, "04_B_3k_ts_time.png"),
  dpi = 150, width = 12, height = 10
)

cat(paste("\nScript completed:", timestamp(suffix = "", prefix = "", quiet = TRUE)))
