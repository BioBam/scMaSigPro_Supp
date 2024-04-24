# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(ggplot2))

# Set paths
dirPath <- "../scMaSigPro_supp_data/ComparisonWithTradeSeq/simulated/sce/"
resPath <- "../scMaSigPro_supp_data/ComparisonWithTradeSeq/output/"
helpScriptsDir <- "R_Scripts/helper_function/"
# dirPath <- "../../../scMaSigPro_supp_data/ComparisonWithTradeSeq/simulated/sce/"
# resPath <- "../../../scMaSigPro_supp_data/ComparisonWithTradeSeq/output/"
# helpScriptsDir <- "../../R_Scripts/helper_function/"

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))

# ReadData
load(paste0(dirPath, "testTradeSeq.RData"))

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
  drop_fac = 20,
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
rm(sim.sce)
gc()

# Initialize a data frame to store memory usage data
memory_data <- data.frame(expr = character(), MemoryBefore = numeric(), MemoryAfter = numeric(), MemoryUsed = numeric(), stringsAsFactors = FALSE)

# Define a helper function to capture memory usage and execute a code block
record_memory_and_run <- function(test_name, code_block) {
  gcinfo(FALSE)
  memory_before <- pryr::mem_used()
  code_block() # Execute the passed function containing the code
  memory_after <- pryr::mem_used()
  gcinfo(FALSE)
  memory_used <- memory_after - memory_before
  memory_data <<- rbind(memory_data, data.frame(expr = test_name, MemoryBefore = memory_before, MemoryAfter = memory_after, MemoryUsed = memory_used))
}

mbm <- microbenchmark(
  "TradeSeq_1_CPU" = record_memory_and_run("TradeSeq_1_CPU", function() {
    cat("\nComputing TradeSeq with 1 CPU\n")
    sce.tradeseq <- fitGAM(counts = normCounts, pseudotime = pseudotime_table, cellWeights = lineage_table, parallel = F, nknots = 3, verbose = FALSE)
    patternRes <- patternTest(sce.tradeseq)
    cat("\nFinished TradeSeq with 1 CPU\n")
  }),
  "TradeSeq_8_CPU" = record_memory_and_run("TradeSeq_8_CPU", function() {
    cat("\nComputing TradeSeq with 8 CPU\n")
    sce.tradeseq <- fitGAM(counts = normCounts, pseudotime = pseudotime_table, cellWeights = lineage_table, parallel = T, nknots = 3, verbose = FALSE)
    patternRes <- patternTest(sce.tradeseq)
    cat("\nFinished TradeSeq with 8 CPU\n")
  }),
  "ScMaSigPro_1_CPU" = record_memory_and_run("ScMaSigPro_1_CPU", function() {
    cat("\nComputing ScMaSigPro with 1 CPU\n")
    scmp.obj <- sc.p.vector(scmpObj = scmp.obj, verbose = FALSE, min_na = 1, parallel = F, offset = T, max_it = 1000)
    scmp.obj <- sc.t.fit(scmpObj = scmp.obj, verbose = FALSE, selection_method = "backward", parallel = F, offset = T)
    cat("\nFinished ScMaSigPro with 1 CPU\n")
  }),
  "ScMaSigPro_8_CPU" = record_memory_and_run("ScMaSigPro_8_CPU", function() {
    cat("\nComputing ScMaSigPro with 8 CPU\n")
    scmp.obj <- sc.p.vector(scmpObj = scmp.obj, verbose = FALSE, min_na = 1, parallel = T, n_cores = 8, offset = T, max_it = 1000)
    scmp.obj <- sc.t.fit(scmpObj = scmp.obj, verbose = FALSE, selection_method = "backward", parallel = FALSE, n_cores = 8, offset = T)
    cat("\nFinished ScMaSigPro with 8 CPU\n")
  }),
  times = 100
)

# Extract Data for time
data <- mbm %>% as.data.frame()

# Extract memory
memory_data <- memory_data %>%
  mutate(MemoryUsedMb = MemoryUsed / 2^20) %>%
  select(expr, MemoryUsedMb)

# Combine data and memory
data <- cbind(data, memory_data)
data <- data[, -3]

# Convert to minutes
data$time <- data$time / 1e9 / 60

# Creating a violin plot for time
compare_violin <- ggplot(data, aes(x = expr, y = time, fill = expr)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.1) +
  labs(title = "Violin Plot of Processing Times", subtitle = "Evaluated 20 times", x = "", y = "Time (minutes)") +
  theme_minimal(base_size = 16) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")

# Save the plot
ggsave(
  plot = compare_violin,
  filename = paste0(resPath, "04_A_tradeSeq_time_violin.png"),
  dpi = 300, width = 9
)
# Save the plot
saveRDS(
  data,
  file = paste0(resPath, "profiling.rds")
)

cat("\nSaved Computing Memory Usage\n")
#
# # Benchmark time
# mbm <- microbenchmark(
#   "TradeSeq_1_CPU" = {
#     # Fit GAM
#       cat("\nComputing TradeSeq with 1 CPU\n")
#     sce.tradeseq <- fitGAM(
#       counts = normCounts,
#       pseudotime = pseudotime_table,
#       cellWeights = lineage_table,
#       parallel = F,
#       nknots = 3, verbose = FALSE
#     )
#     gc()
#     # One of the test
#     patternRes <- patternTest(sce.tradeseq)
#     cat("\nFinished TradeSeq with 1 CPU\n")
#     gc()
#   },
#   "TradeSeq_8_CPU" = {
#       cat("\nComputing TradeSeq with 8 CPU\n")
#     # Fit GAM
#     sce.tradeseq <- fitGAM(
#       counts = normCounts,
#       pseudotime = pseudotime_table,
#       cellWeights = lineage_table,
#       parallel = T,
#       nknots = 3, verbose = FALSE
#     )
#     gc()
#
#     # One of the test
#     patternRes <- patternTest(sce.tradeseq)
#     cat("\nFinished TradeSeq with 8 CPU\n")
#     gc()
#   },
#   "ScMaSigPro_1_CPU" = {
#       cat("\nComputing ScMaSigPro with 1 CPU\n")
#     # Run p-vector
#     scmp.obj <- sc.p.vector(
#       scmpObj = scmp.obj, verbose = FALSE,
#       min_na = 1,
#       parallel = F,
#       offset = T,
#       max_it = 1000
#     )
#     gc()
#
#     # Run-Step-2
#     scmp.obj <- sc.t.fit(
#       scmpObj = scmp.obj, verbose = FALSE,
#       selection_method = "backward", parallel = F,
#       offset = T
#     )
#     gc()
#     cat("\nFinished ScMaSigPro with 1 CPU\n")
#   },
#   "ScMaSigPro_8_CPU" = {
#       cat("\nComputing ScMaSigPro with 8 CPU\n")
#     # Run p-vector
#     scmp.obj <- sc.p.vector(
#       scmpObj = scmp.obj, verbose = FALSE,
#       min_na = 1,
#       parallel = T,
#       n_cores = 8,
#       offset = T,
#       max_it = 1000
#     )
#     gc()
#
#     # Run-Step-2
#     scmp.obj <- sc.t.fit(
#       scmpObj = scmp.obj, verbose = F,
#       selection_method = "backward", parallel = FALSE,
#       n_cores = 8,
#       offset = T
#     )
#     gc()
#     cat("\nFinished ScMaSigPro with 8 CPU\n")
#   },
#   times = 20
# )
#
# # Extract Data
# data <- mbm %>% as.data.frame()
#
# # Convert to minutes
# data$time <- data$time / 1e9 / 60
#
# # Creating a violin plot
# compare_violin <- ggplot(data, aes(x = expr, y = time, fill = expr)) +
#   geom_violin(trim = FALSE) +
#   geom_jitter(width = 0.1) +
#   labs(title = "Violin Plot of Processing Times", subtitle = "Evaluated 20 times", x = "", y = "Time (minutes)") +
#   theme_minimal(base_size = 16) +
#   scale_fill_viridis(
#     alpha = 0.6,
#     discrete = TRUE,
#     breaks = c("TradeSeq_1_CPU", "ScMaSigPro_1_CPU", "TradeSeq_8_CPU", "ScMaSigPro_8_CPU")
#   ) +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
#
#
# # Save
# ggsave(
#   plot = compare_violin,
#   filename = paste0(resPath, "04_A_tradeSeq_time_violin.png"),
#   dpi = 300, width = 9
# )
