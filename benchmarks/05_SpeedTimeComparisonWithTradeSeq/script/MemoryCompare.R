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

# Set paths
dirPath <- "benchmarks/05_SpeedTimeComparisonWithTradeSeq/data/input/"
resPath <- "benchmarks/05_SpeedTimeComparisonWithTradeSeq/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))

# ReadData
load(paste0(dirPath, "sparsity_60.RData"))

# Readuce Dataset
keepGenes <- sample(rownames(rowData(sim.sce)), size = 2500, replace = F)
keepCells <- sample(rownames(colData(sim.sce)), size = 1500, replace = F)

# Extract raw counts
counts <- as.matrix(sim.sce@assays@data@listData$counts)
cell.metadata <- as.data.frame(colData(sim.sce))
gene.metadata <- as.data.frame(rowData(sim.sce))

# Subset the counts
counts.reduced <- counts[keepGenes, keepCells]
cell.metadata.reduced <- cell.metadata[keepCells, ]
gene.metadata.reduced <- gene.metadata[keepGenes, ]

sim.sce <- SingleCellExperiment(list(counts = counts.reduced))
colData(sim.sce) <- DataFrame(cell.metadata.reduced)
rowData(sim.sce) <- DataFrame(gene.metadata.reduced)

# Extract counts
counts <- as.matrix(sim.sce@assays@data@listData$counts)

# Perform Quantile Normalization as per-tradeSeq paper
normCounts <- calcNormCounts(counts, cat = "FQNorm")

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
scmp.obj <-as_scmp(sim.sce, from = "sce",
                   additional_params = list(
                       existing_pseudotime_colname = "Step",
                       existing_path_colname = "Group",
                       overwrite_labels = T), verbose = F)
# Squeeze
scmp.obj <- squeeze(
    scmpObject = scmp.obj,
    bin_method = "Sturges",
    drop.fac = 0.6,
    verbose = F,
    cluster_count_by = "sum"
)

# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
                                  poly_degree = 2
)

# Function to measure memory consumption for a block of code
measure_memory_diff <- function(expr) {
    start_mem <- mem_used()
    eval(expr)
    end_mem <- mem_used()
    return(end_mem - start_mem)
}

# Measure memory consumption for TradeSeq_1_CPU
TradeSeq_1_CPU <- measure_memory_diff(quote({
    sce.tradeseq <- fitGAM(
        counts = normCounts,
        pseudotime = pseudotime_table,
        cellWeights = lineage_table,
        parallel = F,
        nknots = 4, verbose = FALSE
    )
    gc()
    patternRes <- patternTest(sce.tradeseq)
    gc()
}))

# Measure memory consumption for TradeSeq_8_CPU
TradeSeq_8_CPU <- measure_memory_diff(quote({
    sce.tradeseq <- fitGAM(
        counts = normCounts,
        pseudotime = pseudotime_table,
        cellWeights = lineage_table,
        parallel = T,
        nknots = 4, verbose = FALSE
    )
    gc()
    patternRes <- patternTest(sce.tradeseq)
    gc()
}))

# Measure memory consumption for ScMaSigPro_1_CPU
ScMaSigPro_1_CPU <- measure_memory_diff(quote({
    scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj, verbose = F, min.obs = 5,
        counts = T, theta = 10, parallel = F,
        offset = T
    )
    gc()
    scmp.obj <- sc.T.fit(
        data = scmp.obj, verbose = F,
        step.method = "backward",
        family = scmp.obj@scPVector@family,
        offset = T, parallel = F
    )
    gc()
}))

# Measure memory consumption for ScMaSigPro_8_CPU
ScMaSigPro_8_CPU <- measure_memory_diff(quote({
    scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj, verbose = T, min.obs = 5,
        counts = T, theta = 10, parallel = T,
        offset = T
    )
    gc()
    scmp.obj <- sc.T.fit(
        data = scmp.obj, verbose = F,
        step.method = "backward",
        family = scmp.obj@scPVector@family,
        offset = T, parallel = T
    )
    gc()
}))

# Create a dataframe for plotting
mem_data <- data.frame(
    Method = c("TradeSeq_1_CPU", "TradeSeq_8_CPU", "ScMaSigPro_1_CPU", "ScMaSigPro_8_CPU"),
    Memory = c(TradeSeq_1_CPU, TradeSeq_8_CPU, ScMaSigPro_1_CPU, ScMaSigPro_8_CPU)
)

# Convert bytes to MB
mem_data$Memory_MB <- mem_data$Memory / (1024^2)

# Check the maximum memory in the data to decide if we need GB scale
max_mem_gb <- max(mem_data$Memory_MB) / 1024

# If the maximum memory is greater than 1 GB, convert to GB
if(max_mem_gb > 1) {
    mem_data$Memory_MB <- mem_data$Memory_MB / 1024
    y_label <- "Memory Used (GB)"
} else {
    y_label <- "Memory Used (MB)"
}

# Plotting
memBar <- ggplot(mem_data, aes(x = Method, y = Memory_MB, fill = Method)) +
    geom_bar(stat = "identity") +
    labs(title = "Memory Consumption Comparison",
         y = y_label,
         x = "Method") + coord_flip()+
    scale_fill_viridis(discrete = TRUE, name = "Custom Legend Title",
                       breaks = c("TradeSeq_1_CPU", "ScMaSigPro_1_CPU", "TradeSeq_8_CPU", "ScMaSigPro_8_CPU"),
                       labels = c("Custom Label 1", "Custom Label 2", "Custom Label 3", "Custom Label 4")) +
    theme_minimal(base_size = 20) + 
    theme(legend.position = "none", legend.justification = "left", legend.box.just = "left")


# Save
ggsave(
    plot = memBar,
    filename = paste0(resPath, "CompareBarMemory.png"),
    dpi = 800, width = 10
)
