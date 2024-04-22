# Title: Run ScMaSigPro on simulated datasets with different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/01_Sparsity/output/"
outPath <- "/supp_data/benchmarks/01_Sparsity/output/"
outPath2 <- "/supp_data/Tables/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = ".zi.", i = 2),
  ".RData"
)

# Zero-Inflation.evaluation
eval.list <- list()

# Set-up a for loop
for (i in names(dataSets)) {
  # Load data
  scmp_ob <- load(paste0(inPath, "scmp.obj.zi.", i, ".RData"))
  scmpObj <- scmp.obj
  scmp.obj <- NULL

  # Get the sol frame
  solData <- showSol(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE)

  # Extract the counts
  bulkCounts <- as.data.frame(as.matrix(scmpObj@Dense@assays@data$bulk.counts))

  # Create Prediction frame
  predictionData <- data.frame(features = rownames(bulkCounts))

  # Add R2
  predictionData[predictionData[["features"]] %in% rownames(solData), "RSQ"] <- solData[["R-squared"]]

  # Add P-Value
  predictionData[predictionData[["features"]] %in% rownames(solData), "pValue"] <- solData[["p-value"]]

  # Set NA to 0 and 1
  predictionData[is.na(predictionData[["RSQ"]]), "RSQ"] <- 0
  predictionData[is.na(predictionData[["pValue"]]), "pValue"] <- 1

  # Add Remaining columns
  for (i in c(4:ncol(solData))) {
    predictionData[predictionData[["features"]] %in% rownames(solData), colnames(solData[, i, drop = FALSE])] <- solData[[colnames(solData[, i, drop = FALSE])]]
    predictionData[is.na(predictionData[[colnames(solData[, i, drop = FALSE])]]), colnames(solData[, i, drop = FALSE])] <- 1
  }

  # Get the predictions per R2 Value
  sigPredictionData <- predictionData[predictionData[["pValue"]] <= 0.05, , drop = FALSE]

  # Get the values based on RSQ
  sigPredictionData <- sigPredictionData[sigPredictionData[["RSQ"]] >= 0.6, , drop = FALSE]

  # Get column for each pattern
  gene_3_poly_cols <- grep(
    x = colnames(sigPredictionData), ignore.case = TRUE,
    pattern = paste0(scmpObj@Parameters@bin_ptime_col, "3"),
    value = TRUE
  )
  back_gene_3_poly_cols <- colnames(sigPredictionData)[!colnames(sigPredictionData) %in% c(gene_3_poly_cols, "features", "RSQ", "pValue", "p.valor_Path2vsPath1")]
  gene_2_poly_cols <- grep(
    x = colnames(sigPredictionData), ignore.case = TRUE,
    pattern = paste0(scmpObj@Parameters@bin_ptime_col, "2"),
    value = TRUE
  )
  back_gene_2_poly_cols <- colnames(sigPredictionData)[!colnames(sigPredictionData) %in% c(gene_2_poly_cols, "features", "RSQ", "pValue", "p.valor_Path2vsPath1")]
  gene_1_poly_cols <- grep(
    x = colnames(sigPredictionData), ignore.case = TRUE,
    pattern = paste0(scmpObj@Parameters@bin_ptime_col, "x"),
    value = TRUE
  )
  back_gene_1_poly_cols <- colnames(sigPredictionData)[!colnames(sigPredictionData) %in% c(gene_1_poly_cols, "features", "RSQ", "pValue", "p.valor_Path2vsPath1")]

  # Apply conditions
  gene_3_poly <- sigPredictionData[
    rowSums(sapply(sigPredictionData[gene_3_poly_cols], function(x) x <= 0.05)) == length(gene_3_poly_cols) &
      rowSums(sapply(sigPredictionData[back_gene_3_poly_cols], function(x) x >= 0.05)) == length(back_gene_3_poly_cols),
  ]



  plotTrend(
    scmpObj = scmpObj,
    feature_id = gene_3_poly[15, 1],
    significant = FALSE
  )





  gene_3_poly <- sigPredictionData[rowSums(sapply(sigPredictionData[gene_3_poly_cols], function(x) x <= 0.05)) == length(gene_3_poly_cols), ]

  colnames(gene_3_poly) <- colnames(sigPredictionData)
  for (i in grep(
    x = colnames(sigPredictionData), ignore.case = TRUE,
    pattern = paste0(scmpObj@Parameters@bin_ptime_col, "3"), value = TRUE
  )) {
    gene_3_poly <- rbind(gene_3_poly, sigPredictionData[sigPredictionData[[i]] <= 0.5, , drop = FALSE])
  }
  gene_3_poly <- gene_3_poly[-1, , drop = FALSE]
  grep(
    x = colnames(sigPredictionData), ignore.case = TRUE,
    pattern = "3", value = TRUE
  )


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
        scmpObj = scmp.obj, verbose = F,
        min_na = 1, parallel = T,
        offset = TRUE,
        log_offset = TRUE,
        max_it = maxit,
        family = fam
      )

      # Run-Step-2
      scmp.obj <- sc.t.fit(
        parallel = T,
        scmpObj = scmp.obj, verbose = F,
        selection_method = "backward"
      )

      # Save Object
      save(scmp.obj, file = paste0(outPath, "scmp.obj.zi.", i, ".RData"))

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
      performance.measure[["parameter"]] <- "ZI"
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
write.table(evaluation.frame, paste0(outPath2, "01_ZI_Performance.Table.tsv"),
  sep = "\t", row.names = F, quote = F
)
