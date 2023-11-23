# Title: Run ScMaSigPro on simulated datasets with different Skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/06_Normalization/simulated/sce/"
outPath <- "/supp_data/benchmarks/06_Normalization/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
  str_remove(dataSets, pattern = "type_"),
  ".RData"
)

# Set-up a for loop
for (i in names(dataSets)) {
  poly.degree <- 2
  min.gene <- 6
  theta.val <- 10
  ep <- 0.00001

  if (i %in% c("clr", "logN", "rc", "sct")) {
    fam <- gaussian()
  } else {
    fam <- MASS::negative.binomial(10)
  }

  cat(paste("\nRunning for Type:", i))

  # Load Data
  load(file = paste0(inPath, dataSets[i]))

  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as_scmp(sce.obj,
        from = "sce",
        additional_params = list(
          labels_exist = TRUE,
          existing_pseudotime_colname = "Step",
          existing_path_colname = "Group"
        ), verbose = F
      )

      # Compress
      scmp.obj <- squeeze(
        scmpObject = scmp.obj,
        bin_method = "Sturges",
        drop.fac = 1,
        verbose = F,
        cluster_count_by = "sum",
        split_bins = F,
        prune_bins = F,
        drop_trails = F,
        fill_gaps = F
      )

      # Make Design
      scmp.obj <- sc.make.design.matrix(scmp.obj,
        poly_degree = poly.degree
      )

      if (i != "offset") {
        # Run p-vector
        scmp.obj <- sc.p.vector(
          scmpObj = scmp.obj, verbose = F, min.obs = 1,
          offset = F, parallel = T,
          family = fam
        )

        # Run-Step-2
        scmp.obj <- sc.T.fit(
          parallel = T,
          scmpObj = scmp.obj, verbose = F,
          step.method = "backward",
          offset = F
        )
      } else if (i == "offset") {
        # Run p-vector
        scmp.obj <- sc.p.vector(
          scmpObj = scmp.obj, verbose = F, min.obs = 1,
          family = fam,
          offset = T, parallel = T
        )

        # Run-Step-2
        scmp.obj <- sc.T.fit(
          parallel = T,
          scmpObj = scmp.obj, verbose = F,
          step.method = "backward",
          offset = T
        )
      }

      # Save Object
      save(scmp.obj, file = paste0(outPath, "scmp.obj.norm.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}
