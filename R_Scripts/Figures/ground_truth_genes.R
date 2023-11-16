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

# Annotate Genes
rowData(sim.sce) <- DataFrame(add_gene_anno(sim.sce = sim.sce))

similar.change.de <- plot_loess_fit(
    sce_obj = sim.sce, "Gene14", dfreedom = 1, log = T,
    plt_subtitle = "Differentially Expressed"
)+ scale_x_continuous(breaks=(seq(1, 1500, 100)))
opposite.change.de <- plot_loess_fit(
    sce_obj = sim.sce, "Gene691", dfreedom = 1, log = T,
    plt_subtitle = "Differentially Expressed"
)+ scale_x_continuous(breaks=(seq(1, 1500, 100)))
one.change.de <- plot_loess_fit(
    sce_obj = sim.sce, "Gene67", dfreedom = 1, log = T,
    plt_subtitle = "Differentially Expressed"
)+ scale_x_continuous(breaks=(seq(1, 1500, 100)))
no.change <- plot_loess_fit(
    sce_obj = sim.sce, "Gene1", dfreedom = 1, log = T,
    plt_subtitle = "Not-Differentially Expressed"
)+ scale_x_continuous(breaks=(seq(1, 1500, 100)))
# 
gt.true <- ggarrange(similar.change.de, opposite.change.de, one.change.de, no.change,
                     labels = c("A.", "B.", "C.", "D.")
)

ggsave(
    plot = gt.true, filename = paste0(
        "Figures/SuppData/supp_fig_5_Gene_Ground _truth.png"
    ),
    dpi = 600, width = 9
)

# 