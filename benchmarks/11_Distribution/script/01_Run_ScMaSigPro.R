# Title: Run ScMaSigPro on simulated datasets with different lengths of Pseudotime
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "benchmarks/08_Distribution/data/input/"
outPath <- "benchmarks/08_Distribution/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load names of files
load(paste0(inPath, "sparsity_60.RData"))

# Define Avaible distributions
avail.dist <- list(
  "Gaussian" = gaussian(link = "identity"),
  "Poisson" = poisson(link = "log"),
  "Quasi" = quasi(link = "identity", variance = "constant"),
  "Quasipoisson" = quasipoisson(link = "log"),
  "Negative_Binomial_theta_1" = MASS::negative.binomial(1)
)

# Set-up a for loop
for (i in names(avail.dist)) {
  poly.degree <- 2
  min.gene <- 6
  ep <- 0.00001

  cat(paste("\nRunning for Family:", i))

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
        counts = F,
        offset = T, epsilon = ep,
        family = avail.dist[[i]]
      )

      # Run-Step-2
      scmp.obj <- sc.T.fit(
        data = scmp.obj, verbose = F,
        step.method = "backward",
        family = scmp.obj@scPVector@family,
        offset = T
      )

      # Save Object
      save(scmp.obj, file = paste0(outPath, "scmp.obj.Fam.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i))
    }
  )
}
