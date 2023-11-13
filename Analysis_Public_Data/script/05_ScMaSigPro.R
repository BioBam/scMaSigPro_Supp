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
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))

# Load scMaSigPro
#install.packages("../scMaSigPro_1.0.1.tar.gz")
suppressPackageStartupMessages(library(scMaSigPro))

# Load CDS object
load("Analysis_Public_Data/data/rep3/rep3_processed.RData")

# Monocl3 3 object
cds

# Convert the ScMaSigPro Object
scmp.obj <- as_scmp(cds, from = "cds",
                    annotation_colname = "predicted.celltype.l2")

# Compress
scmp.obj <- squeeze(scmpObject = scmp.obj,
                    
                    split_bins = F,
                    prune_bins = F,
                    drop_trails = F)

sc.plot.bins.bar(scmp.obj)
sc.plot.bins.tile(scmp.obj)
sc.fraction.bin(scmp.obj)

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
save(scmp.obj, file = paste0("Analysis_Public_Data/data/rep3/rep3_scMaSigPro_Results.RData"))
