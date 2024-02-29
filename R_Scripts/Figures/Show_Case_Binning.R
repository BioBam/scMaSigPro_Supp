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
        "/supp_data/Figures/SuppData/supp_fig_2_general_bifurcation_trajectory.png"
    ),
    dpi = 600, width = 9
)

# Annotate Genes
rowData(sim.sce) <- DataFrame(add_gene_anno(sim.sce = sim.sce))

# Plot standard gene
opposite.change.de <- plot_loess_fit(
    sce_obj = sim.sce, "Gene691", dfreedom = 1, log = T,
    plt_subtitle = "Differentially Expressed"
)+ scale_x_continuous(breaks = seq(0,max(sim.sce@colData$Step), 100))

# Load ScMaSigPro
library(scMaSigPro)

# Create ScMaSigPro Object
scmp.obj <- as_scmp(sim.sce,
                    from = "sce",
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
    drop_fac = 1,
    assay_name = "counts",
    verbose = F,
    aggregate = "sum",
    split_bins = F,
    prune_bins = F,
    drop_trails = T,
    fill_gaps = T
)

# Plot standard gene
compressed.opposite.change.de <- plot_loess_fit(
    sce_obj = scmp.obj@Dense, "Gene691", dfreedom = 1, log = T,
    plt_subtitle = "Differentially Expressed",
    assay_name = "bulk.counts",
    time_col = scmp.obj@Parameters@bin_ptime_col,
    path_col = scmp.obj@Parameters@path_col
)

trend_bulk_compare <- ggarrange(opposite.change.de,
                                compressed.opposite.change.de,
                                labels = c("A.", "B.")
)

rowData(sim.sce) <- DataFrame(add_gene_anno(sim.sce = sim.sce))

ggsave(
    plot = trend_bulk_compare, filename = paste0(
        "/supp_data/Figures/SuppData/supp_fig_3_trend_compression.png"
    ),
    dpi = 600, width = 9
)
