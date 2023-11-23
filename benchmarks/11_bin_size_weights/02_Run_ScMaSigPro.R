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

# Generate Combination of Parameters
# Create an empty list to store the combinations
param_list <- list()

# Generate all possible combinations of logical values
param_combinations <- expand.grid(
  offset = c(TRUE, FALSE),
  useWeights = c(TRUE, FALSE),
  useInverseWeights = c(TRUE, FALSE),
  useBinWeightAsOffset = c(TRUE, FALSE),
  logOffset = c(TRUE, FALSE)
)

# Loop through each combination and add it to the list
for (i in 1:nrow(param_combinations)) {
  combination <- param_combinations[i, ]

  # Create a name for the combination
  combination_name <- paste0(
    "offset_", ifelse(combination$offset, "T", "F"),
    "_UseWeights_", ifelse(combination$useWeights, "T", "F"),
    "_UseInverseWeights_", ifelse(combination$useInverseWeights, "T", "F"),
    "_UseBinWeightAsOffset_", ifelse(combination$useBinWeightAsOffset, "T", "F"),
    "_logOffset_", ifelse(combination$logOffset, "T", "F")
  )

  # Create a named vector with logical values for true parameters
  param_vector <- sapply(combination, function(value) value)

  # Add the named vector to the list
  param_list[[combination_name]] <- param_vector
}

# Print the named list
param_list

# Set-up a for loop
for (i in names(param_list)) {
  poly.degree <- 2
  min.gene <- 6
  theta.val <- 1
  ep <- 0.00001

  cat(paste("\nRunning for:", i))

  # Load Data
  load(file = paste0(dirPath, dataSets))

  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as_scmp(sim.sce,
        from = "sce",
        additional_params = list(
          labels_exist = TRUE,
          existing_pseudotime_colname = "Step",
          existing_path_colname = "Group"
        ), verbose = F
      )

      # Split bin false
      scmp.obj <- squeeze(
        scmpObject = scmp.obj,
        bin_method = "Sturges",
        drop.fac = 1,
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
        parallel = T,
        offset = as.logical(param_list[[i]][["offset"]]),
        useInverseWeights = as.logical(param_list[[i]][["useInverseWeights"]]),
        useBinWeightAsOffset = as.logical(param_list[[i]][["useBinWeightAsOffset"]]),
        useWeights = as.logical(param_list[[i]][["useWeights"]]),
        logOffset = as.logical(param_list[[i]][["logOffset"]])
      )

      # Run-Step-2
      scmp.obj <- sc.T.fit(
        parallel = T,
        scmpObj = scmp.obj, verbose = F,
        step.method = "backward",
        offset = T
      )

      # Save Object
      save(scmp.obj, file = paste0("benchmarks/12_bin_size_weights/data/output/scmp.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}
