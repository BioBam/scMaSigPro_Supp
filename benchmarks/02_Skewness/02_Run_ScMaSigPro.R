# Title: Run ScMaSigPro on simulated datasets with different Skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/02_Skewness/simulated/sce/"
outPath <- "/supp_data/benchmarks/02_Skewness/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
    str_split_i(dataSets, pattern = "_", i = 2),
    ".RData"
)

# Set-up a for loop
for (i in names(dataSets)) {
    poly.degree <- 2
    drop_fac <- 1
    gTheta<- FALSE
    maxit<- 100
    split.bins = F
  
  if (i > 0.3 | i < 0.7){
      split.bins <- TRUE
  }

  cat(paste("\nRunning for skew:", i))

  # stop("Expected Stop")

  # Load Data
  load(file = paste0(inPath, dataSets[i]))

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
            offset = T, logOffset = T,
            useWeights = T, logWeights = F, useInverseWeights = F,
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
      save(scmp.obj, file = paste0(outPath, "scmp.obj.skew.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}
