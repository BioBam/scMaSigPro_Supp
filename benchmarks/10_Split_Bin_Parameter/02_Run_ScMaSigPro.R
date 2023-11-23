# Title: Run ScMaSigPro on simulated datasets with different Skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/10_Split_bin_Parameter/simulated/sce/"
outPath <- "/supp_data/benchmarks/10_Split_bin_Parameter/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
    str_split_i(dataSets, pattern = "_", i = 2),
    ".RData"
)

# Split_bin_parameterization
spl_bin <- list("SplitBins"= TRUE,
                "noSplitBins" = FALSE)

# Set-up a for loop
for (i in names(dataSets)) {
    
    maxit<- 100
    
    # Get skew value
    skew_value <- dataSets[[i]]
    
    # Load Data
    load(file = paste0(inPath, dataSets[i]))
    
    # Run for splits
    for (j in names(spl_bin)){
        split_bin <- spl_bin[[j]]
        
    cat(paste("\nRunning for skew:", i, "with", j))
    
    tryCatch(
        expr = {
            # Convert
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
                split_bins = split_bin,
                prune_bins = F,
                drop_trails = F,
                fill_gaps = F
            )
            
            # Make Design
            scmp.obj <- sc.make.design.matrix(scmp.obj, poly_degree = 2)
            
            # Run p-vector
            scmp.obj <- sc.p.vector(
                scmpObj = scmp.obj, verbose = F, min.obs = 1, parallel = T,
                offset = T, logOffset = T,
                useWeights = T, logWeights = F, useInverseWeights = F,
                max_it = maxit, 
                globalTheta = FALSE
            )
            
            # Run-Step-2
            scmp.obj <- sc.T.fit(
                parallel = T,
                scmpObj = scmp.obj, verbose = F,
                step.method = "backward"
            )
            
            # Save Object
            save(scmp.obj, file = paste0(outPath, "scmp.obj.skew.", i, ".", j, ".RData"))
            
            # Validate
            cat(paste("\nCompleted for", i, "with", j))
        
        },
        error = function(e) {
            cat(paste("\nFailed for", i, "because", e$message))
        }
    )
    }
}
