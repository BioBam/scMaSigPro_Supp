# FUnction to calculate performance measures
get_performance_ROCR <- function(r2_sequence, groundTruth, scmpObj, include_influ) {
  # Get the sol frame
  solData <- showSol(scmpObj, view = FALSE, return = TRUE, includeInflu = include_influ)

  # Extract the counts
  bulkCounts <- as.data.frame(as.matrix(scmpObj@dense@assays@data$bulk.counts))

  # Create Prediction frame
  predictionData <- data.frame(features = rownames(bulkCounts))

  # Add R2
  predictionData[predictionData[["features"]] %in% rownames(solData), "RSQ"] <- solData[["R-squared"]]

  # Add P-Value
  predictionData[predictionData[["features"]] %in% rownames(solData), "pValue"] <- solData[["p-value"]]

  # Set NA to 0 and 1
  predictionData[is.na(predictionData[["RSQ"]]), "RSQ"] <- 0
  predictionData[is.na(predictionData[["pValue"]]), "pValue"] <- 1

  # Generate perdictions one-by-one
  results <- lapply(r2_sequence, function(i) {
    # Get the predictions per R2 Value
    sigPredictionData <- predictionData[predictionData[["pValue"]] <= 0.05, ]

    # Get the values based on RSQ
    sigPredictionData <- sigPredictionData[sigPredictionData[["RSQ"]] >= i, ]

    # Generate prediction vector
    selected_genes <- sigPredictionData[, "features"]

    # Set predictions to 1
    selected <- rep(1, length(selected_genes))

    # Set names
    names(selected) <- selected_genes

    # Set the unselected genes to false
    unselected_genes <- predictionData[!(predictionData[["features"]] %in% names(selected_genes)), "features"]

    # Set unselected as 0
    unselected <- rep(0, length(unselected_genes))

    # Set names and drop selected
    names(unselected) <- unselected_genes
    unselected <- unselected[!(names(unselected) %in% names(selected))]

    # Predictions
    predictions <- c(selected, unselected)

    # Align according to ground truth
    predictions <- predictions[names(groundTruth)]

    # Calculate binary metrics
    results <- calculate_metrics_binary(
      predictions = predictions,
      groundTruth = groundTruth,
      RSQ = i
    )

    # Return
    return(results)
  })

  # Generate Performance
  performance <- do.call(rbind, results)

  return(performance)
}
