# Title: Run ScMaSigPro on simulated datasets with different lengths of Pseudotime
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/03_Differentiation_Speed/data/simulated/sce/"
dir.create("benchmarks/03_Differentiation_Speed/data/output/", showWarnings = F)
helpScriptsDir <- "R_Scripts/helper_function/"

# Load names of files
dataSets <- list.files(paste0(dirPath))
names(dataSets) <- str_remove(
  str_remove(dataSets, pattern = "Arm_"),
  ".RData"
)

# Set-up a for loop
for (i in names(dataSets)) {
    poly.degree <- 2
    min.gene <- 6
    theta.val <- 1
    ep <- 0.00001
    
    cat(paste("\nRunning for Arm:", i))
    
    #stop("Expected Stop")
    
    # Load Data
    load(file = paste0(dirPath, dataSets[i]))
    
    tryCatch(
        expr = {
            # Convert
            scmp.obj <- as_scmp(sim.sce, from = "sce",
                                align_pseudotime = F,
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
                fill_gaps = F
            )
            sc.plot.bins.tile(scmp.obj)
            
            # Make Design
            scmp.obj <- sc.make.design.matrix(scmp.obj,
                                              poly_degree = poly.degree)
            
            # Run p-vector
            scmp.obj <- sc.p.vector(
                scmpObj = scmp.obj, verbose = F, min.obs = 1,
                offset = T, parallel = T
            )
            
            # Run-Step-2
            scmp.obj <- sc.T.fit(
                parallel = T,
                scmpObj = scmp.obj, verbose = F,
                step.method = "backward",
                offset = T
            )
            
            # Save Object
            save(scmp.obj, file = paste0("benchmarks/03_Differentiation_Speed/data/output/scmp.obj.Arm.", i, ".RData"))
            
            # Validate
            cat(paste("\nCompleted for", i))
        },
        error = function(e) {
            cat(paste("\nFailed for", i, "because", e$message))
        }
    )
}
