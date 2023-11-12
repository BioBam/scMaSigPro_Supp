# Title: Run ScMaSigPro on simulated datasets with different lengths of Pseudotime
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/03_Different_Length/data/simulated/sce/"
dir.create("benchmarks/03_Different_Length/data/output/", showWarnings = F)
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
    drop_fac <- 1
    gTheta<- FALSE
    maxit<- 100
    split.bins = F
    drop_trails  = F
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
            
            if(i %in% c("800_and_2200", "600_and_2400")){
                drop_trails = T
                split.bins = T
            }
            # Compress
            scmp.obj <- squeeze(
                scmpObject = scmp.obj,
                bin_method = "Sturges",
                drop.fac = drop_fac,
                verbose = F,
                cluster_count_by = "sum",
                split_bins = split.bins,
                prune_bins = F,
                drop_trails = drop_trails,
                fill_gaps = F
            )
            # Make Design
            scmp.obj <- sc.make.design.matrix(scmp.obj,
                                              poly_degree = poly.degree)
            
            # Run p-vector
            scmp.obj <- sc.p.vector(
                scmpObj = scmp.obj, verbose = F, min.obs = 1, parallel = T,
                offset = T, logOffset = T,
                useWeights = T, logWeights = T, useInverseWeights = F,
                max_it = maxit, 
                globalTheta = gTheta
            )
            
            # Run-Step-2
            scmp.obj <- sc.T.fit(
                parallel = T,
                scmpObj = scmp.obj, verbose = F,
                step.method = "backward"
            )
            
            # Save Object
            save(scmp.obj, file = paste0("benchmarks/03_Different_Length/data/output/scmp.obj.Arm.", i, ".RData"))
            
            # Validate
            cat(paste("\nCompleted for", i))
        },
        error = function(e) {
            cat(paste("\nFailed for", i, "because", e$message))
        }
    )
}
