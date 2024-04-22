# Title: Run ScMaSigPro on simulated datasets with different levels of skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
base_string <- "/supp_data/benchmarks/"
inputObjectPath <- paste0(base_string, "03_Different_Length/sim/")
outputObjectPath <- paste0(base_string, "03_Different_Length/out/")
tabPath <- paste0(base_string, "03_Different_Length/tab/")
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Create Missing Directory
dir.create(outputObjectPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inputObjectPath))
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = "_", i = 2),
  ".RData"
)

# Zero-Inflation.evaluation
eval.list <- list()

# Set-up a for loop
for (i in names(dataSets)) {
  # Set variables
  poly.degree <- 3
  drop_fac <- 1
  maxit <- 100
  fam <- MASS::negative.binomial(10)

  cat(paste("\nRunning for Length:", i))

  # Load Data
  load(file = paste0(inputObjectPath, dataSets[i]))

  # Evaluate
  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as_scmp(sim.sce,
        from = "sce",
        align_pseudotime = F,
        additional_params = list(
          labels_exist = TRUE,
          exist_ptime_col = "Step",
          exist_path_col = "Group"
        ), verbose = F
      )

      # Compress
      scmp.obj <- sc.squeeze(
        scmpObj = scmp.obj,
        bin_method = "Sturges",
        drop_fac = drop_fac,
        verbose = F,
        aggregate = "sum",
        split_bins = FALSE,
        prune_bins = F,
        drop_trails = F,
        fill_gaps = F
      )

      # Make Design
      scmp.obj <- sc.set.poly(scmp.obj, poly_degree = poly.degree)

      # Run p-vector
      scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj,
        verbose = F,
        min_na = 1,
        parallel = T,
        offset = T,
        log_offset = T,
        max_it = maxit,
        family = fam
      )

      # Run-Step-2
      scmp.obj <- sc.t.fit(
        parallel = T,
        scmpObj = scmp.obj,
        verbose = F,
        selection_method = "backward"
      )

      # Save Object
      save(scmp.obj, file = paste0(outputObjectPath, "scmp.obj.len.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
      # Evaluate
      row_data <- as.data.frame(
        rowData(scmp.obj@Sparse)
      )[, c("gene_short_name", "status")]

      # Set binary labels
      gene.change <- rep(1, length(rownames(row_data[row_data$status != "No_Change", ])))
      gene.no.change <- rep(0, length(rownames(row_data[row_data$status == "No_Change", ])))

      # Add names
      names(gene.change) <- rownames(row_data[row_data$status != "No_Change", ])
      names(gene.no.change) <- rownames(row_data[row_data$status == "No_Change", ])

      # Ground truth
      groundTruth <- c(gene.change, gene.no.change)

      # Get Performance
      performance.measure <- as.data.frame(get_performance_ROCR(
        scmpObj = scmp.obj,
        groundTruth = groundTruth,
        r2_sequence = seq(0.00, 0.95, 0.05),
        include_influ = TRUE
      ))

      # Add to list
      performance.measure[["parameter"]] <- "len"
      performance.measure[["parameter.value"]] <- paste(i, collapse = "_")
      eval.list[[i]] <- performance.measure
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}

# Combine
evaluation.frame <- do.call(rbind, eval.list)

# Write
write.table(evaluation.frame, paste0(tabPath, "03_Length_Performance.Table.tsv"),
  sep = "\t", row.names = F, quote = F
)
