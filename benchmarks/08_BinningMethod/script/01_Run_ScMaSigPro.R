# Title: Run ScMaSigPro on simulated datasets with different lengths of Pseudotime
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "benchmarks/09_BinningMethod/data/input/"
outPath <- "benchmarks/09_BinningMethod/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load names of files
load(paste0(inPath, "sparsity_60.RData"))

# Define Avaible distributions
binning.methods <- list(
  # Sturges = "Sturges",
  # Freedman.Diaconis = "Freedman.Diaconis",
  # Sqrt = "Sqrt",
  # Rice = "Rice",
  # Doane = "Doane",
  Scott.Normal = "Scott.Normal"
)

# Set-up a for loop
for (i in names(binning.methods)) {
  poly.degree <- 2
  min.gene <- 6
  ep <- 0.00001

  cat(paste("\nRunning for Method:", i))

  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as_scmp(sim.sce, from = "sce")

      # Compress
      scmp.obj <- squeeze(
        scmp.ob = scmp.obj,
        time.col = "Step",
        path.col = "Group",
        method = i,
        drop.fac = 0.3,
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
        counts = T,
        offset = T, epsilon = ep,
        family = MASS::negative.binomial(10)
      )

      # Run-Step-2
      scmp.obj <- sc.T.fit(
        data = scmp.obj, verbose = F,
        step.method = "backward",
        family = scmp.obj@scPVector@family,
        offset = T
      )

      # Save Object
      save(scmp.obj, file = paste0(outPath, "scmp.obj.Method.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i))
    }
  )
}
