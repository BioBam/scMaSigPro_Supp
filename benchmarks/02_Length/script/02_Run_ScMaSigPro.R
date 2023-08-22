# Title: Run ScMaSigPro on simulated datasets with different lengths of Pseudotime
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/02_Length/data/simulated/sce/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load names of files
dataSets <- list.files(paste0(dirPath))
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = "_", i = 2),
  ".RData"
)

# Set-up a for loop
for (i in names(dataSets)) {
    
    poly.degree <- 2
    min.gene <- 6
    theta.val <- 1
    ep <- 0.00001
    
  cat(paste("\nRunning for LenEq:", i))

  # stop("Expected Stop")

  # Load Data
  load(file = paste0(dirPath, dataSets[i]))

  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as_scmp(sim.sce, from = "sce")

      # Compress
      scmp.obj <- squeeze(
        scmp.ob = scmp.obj,
        time.col = "Step",
        path.col = "Group",
        method = "Sturges",
        drop.fac = 0.6,
        verbose = T,
        cluster.count.by = "sum"
      )

      # Make Design
      scmp.obj <- sc.make.design.matrix(scmp.obj,
        degree = poly.degree,
        time.col = "binnedTime",
        path.col = "path"
      )

      # Run p-vector
      scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj, verbose = F, min.obs = min.gene,
        counts = T, theta = theta.val,
        offset = T, epsilon = ep
      )

      # Run-Step-2
      scmp.obj <- sc.T.fit(
        data = scmp.obj, verbose = F,
        step.method = "backward",
        family = scmp.obj@scPVector@family,
        offset = T
      )

      # Save Object
      save(scmp.obj, file = paste0("benchmarks/02_Length/data/output/scmp.obj.LenEq.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i))
    }
  )
}
