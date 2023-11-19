# Title: Run ScMaSigPro on simulated datasets with different Binning
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/07_BinningMethod/simulated/sce/"
outPath <- "/supp_data/benchmarks/07_BinningMethod/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
    str_split_i(dataSets, pattern = "_", i = 2),
    ".RData"
)
# Load Data
load(file = paste0(inPath, dataSets["50"]))

meth.list <- list("Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane", "Scott.Normal")
names(meth.list) <- meth.list

#meth.list <- meth.list[names(meth.list) %in% c("Freedman.Diaconis", "Scott.Normal")]

# Set-up a for loop
for (i in names(meth.list)) {
    poly.degree <- 2
    drop_fac <- 1
    gTheta<- FALSE
    maxit<- 100
    fam <- MASS::negative.binomial(10)
    split.bins = F

    cat(paste("\nRunning for sparsity:", i))
    
    if(i %in% c("Freedman.Diaconis", "Scott.Normal")){
        drop_fac <- 0.3
    }
    
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
                bin_method = i,
                drop.fac = drop_fac,
                verbose = F,
                cluster_count_by = "sum",
                split_bins = split.bins,
                prune_bins = F,
                drop_trails = F,
                fill_gaps = F
            )
            
            # Make Design
            scmp.obj <- sc.make.design.matrix(scmp.obj,
                                              poly_degree = poly.degree)
            
            # Run p-vector
            scmp.obj <- sc.p.vector(
                scmpObj = scmp.obj, verbose = F, min.obs = 1, parallel = T,
                offset = T, 
                logOffset = F,
                useWeights = T, 
                logWeights = F, 
                useInverseWeights = F,
                max_it = maxit, 
                family = fam
            )
            
            # Run-Step-2
            scmp.obj <- sc.T.fit(
                parallel = T,
                scmpObj = scmp.obj, verbose = F,
                step.method = "backward"
            )
            
            # Save Object
            save(scmp.obj, file = paste0(outPath, "scmp.obj.method.", i, ".RData"))
            
            # Validate
            cat(paste("\nCompleted for", i))
        },
        error = function(e) {
            cat(paste("\nFailed for", i, "because", e$message))
        }
    )
}
