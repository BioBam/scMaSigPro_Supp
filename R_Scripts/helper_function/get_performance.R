get_performance <- function(r2_sequence, gene_no_change,
                            gene_change, scmp_obj) {
  # Create Dataframe to plot Performance
  df.performance <- data.frame(
    VARIABLE = 0, TP = 0, FP = 0, TN = 0, FN = 0,
    TPR = 0, FPR = 0, FNR = 0, TNR = 0,
    INFLU = 0, INFLU_FP = 0, INFLU_FN = 0, ACCURACY = 0,
    PRECISION = 0
  )

  all_genes <- c(gene_no_change, gene_change)

  # Run For all Values of R square
  for (j in r2_sequence) {
    # Get Sig genes
    capture.output(
      scmp_obj <- scMaSigPro::sc.get.siggenes(
        scmpObj = scmp_obj,
        rsq = j,
        vars = "groups"
      )
    )

    sigs <- list(
      sig.genes = scmp_obj@sig.genes@sig.genes,
      summary = scmp_obj@sig.genes@summary
    )


    if (is.null(sigs$sig.genes)) {
      classifications <- list(
        detected.genes = "No Significant Genes found",
        undetected.gene = all_genes
      )
    } else {
      # Path2 and Path3vsPath2
      shared <- intersect(sigs$summary[, 2], sigs$summary[, 1])

      # Path2
      grp1 <- setdiff(sigs$summary[, 1], sigs$summary[, 2])

      # Path3vsPath2
      grp2 <- setdiff(sigs$summary[, 2], sigs$summary[, 1])

      # All detected
      detected <- c(shared, grp2, grp1)

      # Remove empty element
      detected <- detected[!(detected == " ")]

      # Get the undetected genes
      undetected <- all_genes[!(all_genes %in% detected)]

      classifications <- list(detected.genes = detected, undetected.gene = undetected)
    }

    # True Positive
    tp <- intersect(gene_change, classifications$detected.genes)

    # True Negative
    tn <- intersect(gene_no_change, classifications$undetected.gene)

    # False Positive
    fp <- intersect(gene_no_change, classifications$detected.genes)

    # False Negative
    fn <- intersect(gene_change, classifications$undetected.gene)

    # Influential Genes
    influ <- colnames(scmp_obj@scTFit@influ.info)

    # Influential + False Positive
    fp_influential <- intersect(fp, colnames(scmp_obj@scTFit@influ.info))
    fn_influential <- intersect(fn, colnames(scmp_obj@scTFit@influ.info))

    # Sensitivity, True Positive Rate / Recall
    senstivity <- length(tp) / (length(tp) + length(fn))

    # False Negative Rate
    fnr <- length(fn) / (length(tp) + length(fn))

    # Specificity / True Negative Rate
    specificity <- length(tn) / (length(tn) + length(fp))

    # False Positive Rate
    fpr <- 1 - specificity

    # Accuracy
    accuracy <- (length(tp) + length(tn)) / (length(tp) + length(tn) + length(fp) + length(fn))

    # Calculate Precison
    precision <- length(tp) / (length(tp) + length(fp))

    # Make the vector of results
    res.perfom <- data.frame(
      VARIABLE = j,
      TP = length(tp),
      FP = length(fp),
      TN = length(tn),
      FN = length(fn),
      TPR = round(senstivity, 3),
      FPR = round(fpr, 3),
      FNR = round(fnr, 3),
      TNR = round(specificity, 3),
      INFLU = length(influ),
      INFLU_FP = length(fp_influential),
      INFLU_FN = length(fn_influential),
      ACCURACY = round(accuracy, 3),
      PRECISION = round(precision, 3)
    )

    # Add vector to the frame
    df.performance <- rbind(df.performance, res.perfom)
    gc()
  }

  # Remove the first Row
  df.performance <- df.performance[-1, ]

  return(df.performance)
}
