# Title: Run ScMaSigPro on Normalization
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(MASS))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/06_Offset/simulated/sce/"
outPath <- "/supp_data/benchmarks/06_Offset/output/"
outPath2 <- "Tables/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
dataSets <- dataSets[dataSets %in% grep(x = dataSets, pattern = "norm", value = T)]
names(dataSets) <- str_remove(str_remove(
  dataSets,
  ".RData"
), "norm_")

# Zero-Inflation.evaluation
eval.list <- list()

# Set-up a for loop
for (i in names(dataSets)) {
  # Set variables
  poly.degree <- 2
  drop_fac <- 1
  maxit <- 100

  cat(paste("\nRunning for sparsity:", i))

  # Change parameters according to the normalization strategy
  if (i == "SCT") {
    cat(paste("Performing Seurat's SCT"))
    use_Offset <- FALSE
    fam <- gaussian()
  } else if (i == "libSize") {
    cat(paste("Performing Seurat's libSize"))
    use_Offset <- FALSE
    fam <- gaussian()
  } else if (i == "logLibSize") {
    cat(paste("Performing Seurat's logLibSize"))
    use_Offset <- FALSE
    fam <- gaussian()
  } else if (i == "rawCounts") {
    cat(paste("Using Raw Counts"))
    use_Offset <- FALSE
    fam <- negative.binomial(20)
  } else if (i == "FQ") {
    cat(paste("Performing FQ"))
    use_Offset <- FALSE
    fam <- MASS::negative.binomial(20)
  } else if (i == "rawCounts_Offset") {
    cat(paste("rawCounts_Offset"))
    use_Offset <- TRUE
    fam <- MASS::negative.binomial(20)
  }

  # Load Data
  load(file = paste0(inPath, dataSets[i]))

  # Evaluate
  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as.scmp(sim_sce,
        from = "sce",
        align_pseudotime = F,
        additional_params = list(
          labels_exist = TRUE,
          existing_pseudotime_colname = "Step",
          existing_path_colname = "Group"
        ), verbose = F
      )

      # Compress
      scmp.obj <- sc.squeeze(
        scmpObject = scmp.obj,
        bin_method = "Sturges",
        drop_fac = drop_fac,
        verbose = F,
        cluster_count_by = "sum",
        split_bins = F,
        prune_bins = F,
        drop_trails = F,
        fill_gaps = F
      )

      # Make Design
      scmp.obj <- sc.set.poly(scmp.obj, poly_degree = poly.degree)

      # Run p-vector
      scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj, verbose = F,
        min.na = 1, parallel = F,
        offset = use_Offset,
        logOffset = F,
        max_it = maxit,
        family = fam
      )

      # Run-Step-2
      scmp.obj <- sc.t.fit(
        parallel = T,
        scmpObj = scmp.obj, verbose = F,
        step.method = "backward"
      )

      # Save Object
      save(scmp.obj, file = paste0(outPath, "scmp.obj.norm.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
      # Evaluate
      row_data <- as.data.frame(
        rowData(scmp.obj@sparse)
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
      performance.measure[["parameter"]] <- "Norm"
      performance.measure[["parameter.value"]] <- i
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
write.table(evaluation.frame, paste0(outPath2, "04_Normalization.Table.tsv"),
  sep = "\t", row.names = F, quote = F
)
