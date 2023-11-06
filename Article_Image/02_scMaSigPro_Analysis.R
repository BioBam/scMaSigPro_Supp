# Title: Run ScMaSigPro on simulated datasets with different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "Article_Image/data/sce/"
dir.create("Article_Image/data/output/", showWarnings = F)
helpScriptsDir <- "R_Scripts/helper_function/"

# Load data
load(file = paste0(dirPath, "ArticleFigure.RData"))

# Create Object
scmp.obj <- as_scmp(sim.sce, from = "sce",
                    additional_params = list(
                        labels_exist = TRUE,
                        existing_pseudotime_colname = "Step",
                        existing_path_colname = "Group"),
                    verbose = F)
            
            
# Compress
scmp.obj <- squeeze(
    scmpObject = scmp.obj,
    bin_method = "Doane",
    drop.fac = 0.7, verbose = F,
    cluster_count_by = "sum",
    split_bins = T,
    prune_bins = F,
    drop_trails = T,
    fill_gaps = F
    )

sc.plot.bins.tile(scmp.obj)
sc.plot.bins.bar(scmp.obj)

# Make Design
scmp.obj <- sc.make.design.matrix(
    scmpObject = scmp.obj,
    poly_degree = 1)

# Run p-vector
scmp.obj <- sc.p.vector(
    scmpObj = scmp.obj, verbose = F, min.obs = 1,
    offset = T, parallel = T,
    family = MASS::negative.binomial(theta = 0.5, link = "log")#logit
    )

# Run-Step-2
scmp.obj <- sc.T.fit(
    parallel = T,
    scmpObj = scmp.obj, verbose = F,
    step.method = "backward",
    offset = T
    )
            
# Save Object
save(scmp.obj, file = paste0("Article_Image/data/output/scmp.obj.article.image.RData"))
