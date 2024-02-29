# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(ggpubr))

# Load Custom Function
source("R_Scripts/helper_function/plot_simulations().R")
source("R_Scripts/helper_function/add_gene_anno().R")
source("R_Scripts/helper_function/plot_loess.R")

# Simulate a general splatter dataset
params.groups <- newSplatParams(
  batch.rmEffect = TRUE, # No Batch affect
  batchCells = 3000, # Number of Cells
  nGenes = 5000, # Number of Genes
  seed = 2022, # Set seed
  mean.rate = 0.3, mean.shape = 5, lib.scale = 0.2,
  lib.loc = 12, dropout.mid = 0,
  dropout.type = "experiment",
  group.prob = c(0.5, 0.5), path.from = c(0, 0),
  de.prob = 0.3, de.facLoc = 1, path.nonlinearProb = 0,
  path.sigmaFac = 0,
  dropout.shape = -0.4,
  path.nSteps = c(1500, 1500)
)

# Simulate the dataset
sim.sce <- splatSimulate(params.groups,
  method = "paths",
  verbose = F,
  dropout.shape = 0
)

# Drop cells
drop_cells <- sim.sce@colData[(sim.sce@colData$Step >= 600 & sim.sce@colData$Step <= 1000 & sim.sce@colData$Group == "Path1" |
  sim.sce@colData$Step >= 200 & sim.sce@colData$Step <= 500 & sim.sce@colData$Group == "Path1" |
  sim.sce@colData$Step >= 400 & sim.sce@colData$Step <= 700 & sim.sce@colData$Group == "Path2" |
  sim.sce@colData$Step >= 1100 & sim.sce@colData$Step <= 1400 & sim.sce@colData$Group == "Path2"), ][["Cell"]]
sim.sce <- sim.sce[, !(colnames(sim.sce) %in% drop_cells)]

# Plot PCA of the simulated dataset
step.plot <- plot_simulations(sim.sce,
  assay_type = "TrueCounts",
  plot3d = F, plot2d = T, frame = 2,
  title.2d = "Trajectory with missing intermediate cell states"
)
# Annotate Genes
rowData(sim.sce) <- DataFrame(add_gene_anno(sim.sce = sim.sce))

# Plot standard gene
opposite.change.de <- plot_loess_fit(
  sce_obj = sim.sce, "Gene691", dfreedom = 2, log = T,
  span = 0.1,
  plt_subtitle = "Differentially Expressed"
) + scale_x_continuous(breaks = (seq(1, 1500, 100)))
opposite.change.de

# Load ScMaSigPro
library(scMaSigPro)

# Create ScMaSigPro Object
scmp.obj <- as_scmp(sim.sce,
  from = "sce",
  additional_params = list(
    labels_exist = TRUE,
    existing_pseudotime_colname = "Step",
    existing_path_colname = "Group"
  ), verbose = F
)

# Compress
scmp.obj.gaps <- squeeze(
  scmpObject = scmp.obj,
  bin_method = "Sturges",
  drop.fac = 1,
  assay_name = "counts",
  verbose = F,
  cluster_count_by = "sum",
  split_bins = F,
  prune_bins = F,
  drop_trails = F,
  fill_gaps = F
)

# Compress
scmp.obj.no.gaps <- squeeze(
  scmpObject = scmp.obj,
  bin_method = "Sturges",
  drop.fac = 1,
  assay_name = "counts",
  verbose = F,
  cluster_count_by = "sum",
  split_bins = F,
  prune_bins = F,
  drop_trails = F,
  fill_gaps = T
)


# Compress
scmp.obj.no.gaps.tail <- squeeze(
  scmpObject = scmp.obj,
  bin_method = "Sturges",
  drop.fac = 1,
  assay_name = "counts",
  verbose = F,
  cluster_count_by = "sum",
  split_bins = F,
  prune_bins = F,
  drop_trails = T,
  fill_gaps = T
)

gap.tile <- sc.plot.bins.tile(scmp.obj.gaps)
no.gap.tile <- sc.plot.bins.tile(scmp.obj.no.gaps)
no.gap.tail.tile <- sc.plot.bins.tile(scmp.obj.no.gaps.tail)


# Plot standard gene
compressed.opposite.change.de.no.gaps.tail <- plot_loess_fit(
  sce_obj = scmp.obj.no.gaps.tail@compress.sce, "Gene691", dfreedom = 1, log = T,
  plt_subtitle = "Differentially Expressed",
  assay_name = "bulk.counts",
  time_col = scmp.obj.no.gaps.tail@addParams@bin_pseudotime_colname,
  path_col = scmp.obj.no.gaps.tail@addParams@path_colname
)

compare_bin_plot <- ggarrange(step.plot, opposite.change.de, compressed.opposite.change.de.no.gaps.tail,
  gap.tile, no.gap.tile, no.gap.tail.tile,
  labels = c("A.", "B.", "C,", "D.", "E.", "F.")
)

compare_bin_plot

ggsave(
  plot = compare_bin_plot, filename = paste0(
    "/supp_data/Figures/SuppData/supp_fig_4_NoGap_NoTails_compression.png"
  ),
  dpi = 600, width = 20, height = 8,
)
