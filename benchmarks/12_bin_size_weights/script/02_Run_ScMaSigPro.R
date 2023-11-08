# Title: Run ScMaSigPro on simulated datasets with different Skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/12_bin_size_weights/data/input/sce/"
dir.create("benchmarks/12_bin_size_weights/data/output/", showWarnings = F)
helpScriptsDir <- "R_Scripts/helper_function/"

# Load names of files
dataSets <- list.files(paste0(dirPath))

# Set-up a for loop
for (i in c("TRUE", "FALSE")) {
    poly.degree <- 2
    min.gene <- 6
    theta.val <- 1
    ep <- 0.00001
    
    cat(paste("\nUsing bin size as weights:", i))
    
    # Load Data
    load(file = paste0(dirPath, dataSets))
    
    tryCatch(
        expr = {
            # Convert
            scmp.obj <- as_scmp(sim.sce, from = "sce",
                                additional_params = list(
                                    labels_exist = TRUE,
                                    existing_pseudotime_colname = "Step",
                                    existing_path_colname = "Group"), verbose = F)

            # Split bin false
            scmp.obj <- squeeze(
                scmpObject = scmp.obj,
                bin_method = "Sturges",
                drop.fac = 0.7,
                verbose = F,
                cluster_count_by = "sum",
                split_bins = FALSE,
                prune_bins = F,
                drop_trails = F,
                fill_gaps = F
            )
            
            # Make Design
            scmp.obj <- sc.make.design.matrix(scmp.obj, poly_degree = poly.degree)
            
            # Run p-vector
            scmp.obj <- sc.p.vector(
                scmpObj = scmp.obj, verbose = F, min.obs = 1,
                offset = T, parallel = T, useWeights = as.logical(i)
            )
            
            # Run-Step-2
            scmp.obj <- sc.T.fit(
                parallel = T,
                scmpObj = scmp.obj, verbose = F,
                step.method = "backward",
                offset = T
            )

            # Save Object
            save(scmp.obj, file = paste0("benchmarks/12_bin_size_weights/data/output/scmp.obj.bin.weights.", i, ".RData"))
            
            # Validate
            cat(paste("\nCompleted for", i))
        },
        error = function(e) {
            cat(paste("\nFailed for", i, "because", e$message))
        }
    )
}
