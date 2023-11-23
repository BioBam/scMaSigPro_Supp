##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: scMaSigPro Application ########################
##########################################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))

# Load scMaSigPro
suppressPackageStartupMessages(library(scMaSigPro))

# Monocl3 3 object
cds <- readRDS("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/Monocle3_Processed_Donor3.RDS")

# Convert the ScMaSigPro Object
scmp.obj <- as_scmp(cds, from = "cds",
                    align_pseudotime = F,
                    annotation_colname = "cell_type")

# Compress
scmp.obj <- squeeze(scmpObject = scmp.obj,
                    split_bins = F,
                    prune_bins = F,
                    drop_trails = F)

a <- sc.plot.bins.bar(scmp.obj)
b <- sc.plot.bins.tile(scmp.obj)
a +b


# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
                                  poly_degree = 3,
)

scmp.obj@distribution <- MASS::negative.binomial(10)

# Run p-vector
scmp.obj <- sc.p.vector(
    parallel = T,
    scmpObj = scmp.obj, verbose = T,
    max_it = 10000,
    globalTheta = T,
    logOffset = T,
    useInverseWeights = F,
    logWeights = F,
    offset = T,
    min.obs = 1)

# Run-Step-2
scmp.obj <- sc.T.fit(
    scmpObj = scmp.obj, verbose = T,
    step.method = "backward"
)

# Save the object for further analysis
save(scmp.obj, file = paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/scMaSigPro_Donor2_Monocle3.RData"))
