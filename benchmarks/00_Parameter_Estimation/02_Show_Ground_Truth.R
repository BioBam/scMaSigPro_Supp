# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(coop))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))

# Set paths
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
rdsPath <- paste0(base_string, "benchmarks/00_Parameter_Estimation/input/")
imgPath <- paste0(base_string, "benchmarks/00_Parameter_Estimation/img/")
outPath <- paste0(base_string, "benchmarks/00_Parameter_Estimation/output/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Create Directory if does not exist
dir.create(figPath, showWarnings = FALSE, recursive = TRUE)
dir.create(imgPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_hd, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_lr, showWarnings = FALSE, recursive = TRUE)
dir.create(tabPath, showWarnings = FALSE, recursive = TRUE)
dir.create(rdsPath, showWarnings = FALSE, recursive = TRUE)
dir.create(outPath, showWarnings = FALSE, recursive = TRUE)

# Load Custom Function
source(paste0(helpScriptsDir, "plot_simulations().R"))
source(paste0(helpScriptsDir, "add_gene_anno().R"))
source(paste0(helpScriptsDir, "plot_loess.R"))

# Load
paramEstimates <- readRDS(paste0(base_string, "benchmarks/00_Parameter_Estimation/output/setty_et_al_d1_splatEstimates.RDS"))

# Create Base parameters/ Same for All groups
params.groups <- newSplatParams(
  batch.rmEffect = TRUE, # No Batch affect
  batchCells = 3000, # Number of Cells
  nGenes = 2000, # Number of Genes
  seed = 2022, # Set seed
  mean.rate = paramEstimates@mean.rate,
  mean.shape = paramEstimates@mean.shape,
  lib.scale = paramEstimates@lib.scale,
  lib.loc = paramEstimates@lib.loc,
  bcv.common = paramEstimates@bcv.common,
  bcv.df = paramEstimates@bcv.df,
  dropout.type = "experiment",
  group.prob = c(0.5, 0.5),
  path.from = c(0, 0),
  de.prob = 0.3,
  de.facLoc = 1,
  out.facLoc = paramEstimates@out.facLoc,
  dropout.mid = paramEstimates@dropout.mid,
  out.facScale = paramEstimates@out.facScale,
  out.prob = paramEstimates@out.prob,
  path.skew = c(0.5, 0.5),
  path.nSteps = c(1500, 1500),
  dropout.shape = paramEstimates@ dropout.shape
)

# Simulate the dataset
sim.sce <- splatSimulate(params.groups,
  method = "paths",
  verbose = F,
)

# Annotate Genes
rowData(sim.sce) <- DataFrame(add_gene_anno(sim.sce = sim.sce))

# Estimate
trueSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
simulatedSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
totSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)

# Plot Simulation and ground truth
sob.tc <- CreateSeuratObject(
  counts = sim.sce@assays@data@listData$TrueCounts,
  meta.data = as.data.frame(sim.sce@colData)
)
sob.sim <- CreateSeuratObject(
  counts = sim.sce@assays@data@listData$counts,
  meta.data = as.data.frame(sim.sce@colData)
)
sob.tc <- NormalizeData(sob.tc,
  normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = F
)
sob.sim <- NormalizeData(sob.sim,
  normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = F
)
sob.tc <- FindVariableFeatures(sob.tc,
  selection.method = "vst", nfeatures = 2000,
  verbose = F
)
sob.sim <- FindVariableFeatures(sob.sim,
  selection.method = "vst", nfeatures = 2000,
  verbose = F
)
sob.tc <- ScaleData(sob.tc, verbose = F)
sob.sim <- ScaleData(sob.sim, verbose = F)
sob.tc <- RunPCA(sob.tc, features = VariableFeatures(object = sob.tc), verbose = F)
sob.sim <- RunPCA(sob.sim, features = VariableFeatures(object = sob.sim), verbose = F)
sob.tc <- RunUMAP(sob.tc, dims = 1:10, verbose = F)
sob.sim <- RunUMAP(sob.sim, dims = 1:10, verbose = F)

# Create Plotting frame for PHATE
plt.data.tc <- data.frame(
  UMAP_1 = sob.tc@reductions[["umap"]]@cell.embeddings[, 1],
  UMAP_2 = sob.tc@reductions[["umap"]]@cell.embeddings[, 2],
  Simulated_Steps = sim.sce@colData$Step,
  Path = sim.sce@colData$Group
)
plt.data.sim <- data.frame(
  UMAP_1 = sob.sim@reductions[["umap"]]@cell.embeddings[, 1],
  UMAP_2 = sob.sim@reductions[["umap"]]@cell.embeddings[, 2],
  Simulated_Steps = sim.sce@colData$Step,
  Path = sim.sce@colData$Group
)

# Plot UMAP dimensions
plt.sim <- ggplot(plt.data.sim) +
  geom_point(
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = Simulated_Steps,
      shape = Path
    ),
    size = 1.5
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_color_viridis(option = "C") +
  ggtitle(
    paste("Base Simulation"),
    subtitle = paste(
      "True Sparsity:", trueSparsity,
      "| Simulated:", simulatedSparsity,
      "| Total:", totSparsity
    )
  )


similar.change.de <- plot_loess_fit(
  sce_obj = sim.sce, "Gene659", dfreedom = 1, log = T,
  plt_subtitle = "DE: Similar Change", point.alpha = 0.2,
  point.size = 1, line.size = 1, line.alpha = 1, ci = F,
) + scale_x_continuous(breaks = (seq(0, 1500, 300))) + labs(
  color = "Branching Path", shape = "Branching Path"
)

opposite.change.de <- plot_loess_fit(
  sce_obj = sim.sce, "Gene798", dfreedom = 1, log = T,
  plt_subtitle = "DE: Opposite Change", point.alpha = 0.2,
  point.size = 1, line.size = 1, line.alpha = 1, ci = F,
) + scale_x_continuous(breaks = (seq(0, 1500, 300))) + labs(
  color = "Branching Path", shape = "Branching Path"
)

one.change.de <- plot_loess_fit(
  sce_obj = sim.sce, "Gene309", dfreedom = 1, log = T,
  plt_subtitle = "DE: Change in One", point.alpha = 0.2,
  point.size = 1, line.size = 1, line.alpha = 1, ci = F,
) + scale_x_continuous(breaks = (seq(0, 1500, 300))) + labs(
  color = "Branching Path", shape = "Branching Path"
)
no.change <- plot_loess_fit(
  sce_obj = sim.sce, "Gene336", dfreedom = 1, log = T,
  plt_subtitle = "Not-Differentially Expressed", point.alpha = 0.2,
  point.size = 1, line.size = 1, line.alpha = 1, ci = F,
) + scale_x_continuous(breaks = (seq(0, 1500, 300))) + labs(
  color = "Branching Path", shape = "Branching Path"
)

gt.true <- ggarrange(similar.change.de, opposite.change.de, one.change.de, no.change,
  labels = c("A.", "B.", "C.", "D."), common.legend = T, legend = "bottom"
)

# Save
saveRDS(gt.true,
  file = paste0(imgPath, "MainFigure2A.rds")
)

gt.true <- ggarrange(similar.change.de, opposite.change.de, one.change.de, no.change,
  labels = c("B.", "C.", "D.", "E."), common.legend = T, legend = "bottom"
)
ground.truth <- ggarrange(plt.sim, gt.true, labels = c("A.", ""))
ground.truth

ggsave(
  plot = ground.truth, filename = paste0(figPath_hd, "supp_fig_2_Sim_and_Ground_Truth.png"),
  dpi = 600, width = 12, height = 5
)
ggsave(
  plot = ground.truth, filename = paste0(figPath_lr, "supp_fig_2_Sim_and_Ground_Truth.png"),
  dpi = 150, width = 12, height = 5
)
