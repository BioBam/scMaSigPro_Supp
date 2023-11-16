# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(ggpubr))

# Load Custom Function
source("R_Scripts/helper_function/plot_simulations().R")
source("R_Scripts/helper_function/add_gene_anno().R")

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

# Plot PCA of the simulated dataset
step.plot <- plot_simulations(sim.sce,
                              assay_type = "TrueCounts",
                              plot3d = F, plot2d = T, frame = 2,
                              title.2d = "A. Bifurcating Trajectory"
)
group.plot <- plot_simulations(sim.sce,
                               assay_type = "TrueCounts",
                               plot3d = F, plot2d = T, frame = 1,
                               title.2d = "B. Bifurcating Trajectory"
)

general.plot <- ggarrange(step.plot, group.plot)

ggsave(
    plot = general.plot, filename = paste0(
        "Figures/SuppData/supp_fig_2_general_bifurcation_trajectory.png"
    ),
    dpi = 600, width = 9
)

# Annotate Genes
rowData(sim.sce) <- DataFrame(add_gene_anno(sim.sce = sim.sce))

# Plot standard gene
opposite.change.de <- plot_loess_fit(
    sce_obj = sim.sce, "Gene691", dfreedom = 1, log = T,
    plt_subtitle = "Differentially Expressed"
)

# Load ScMaSigPro
library(scMaSigPro)

# Create ScMaSigPro Object
scmp.obj <- as_scmp(sim.sce, from = "sce",
                    additional_params = list(
                        labels_exist = TRUE,
                        existing_pseudotime_colname = "Step",
                        existing_path_colname = "Group"), verbose = F)

# Compress
scmp.obj <- squeeze(
    scmpObject = scmp.obj,
    bin_method = "Sturges",
    drop.fac = 1,
    verbose = F,
    cluster_count_by = "sum",
    split_bins = F,
    prune_bins = F,
    drop_trails = T,
    fill_gaps = T
)

# Plot standard gene
compressed.opposite.change.de <- plot_loess_fit(
    sce_obj = scmp.obj@compress.sce, "Gene691", dfreedom = 1, log = T,
    plt_subtitle = "Differentially Expressed", 
    assay_name = "bulk.counts",
    time_col = scmp.obj@addParams@bin_pseudotime_colname,
    path_col = scmp.obj@addParams@path_colname
)

trend_bulk_compare <- ggarrange(opposite.change.de,
          compressed.opposite.change.de,
          labels = c("A.", "B."))

rowData(sim.sce) <- DataFrame(add_gene_anno(sim.sce = sim.sce))

ggsave(
    plot = trend_bulk_compare, filename = paste0(
        "Figures/SuppData/supp_fig_3_trend_compression.png"
    ),
    dpi = 600, width = 9
)

# 
# similar.change.de <- plot_loess_fit(
#     sce_obj = sim.sce, "Gene14", dfreedom = 1, log = T,
#     plt_subtitle = "Differentially Expressed"
# )
# opposite.change.de <- plot_loess_fit(
#     sce_obj = sim.sce, "Gene691", dfreedom = 1, log = T,
#     plt_subtitle = "Differentially Expressed"
# )
# one.change.de <- plot_loess_fit(
#     sce_obj = sim.sce, "Gene67", dfreedom = 1, log = T,
#     plt_subtitle = "Differentially Expressed"
# )
# no.change <- plot_loess_fit(
#     sce_obj = sim.sce, "Gene1", dfreedom = 1, log = T,
#     plt_subtitle = "Not-Differentially Expressed"
# )
# 
# gt.true <- ggarrange(similar.change.de, opposite.change.de, one.change.de, no.change,
#                      labels = c("A.", "B.", "C.", "D.")
# )
# 
# ggsave(
#     plot = gt.true, filename = paste0(
#         imgDir,
#         "Fig2_What_is_differentially_Expressed.png"
#     ),
#     dpi = 600, height = 8, width = 9
# )
