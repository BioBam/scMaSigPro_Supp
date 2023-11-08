# Get ROCR
suppressPackageStartupMessages(library(ROCR))

# Calculate binary
calculate_metrics_binary <- function(predictions, groundTruth, RSQ) {
    # Ensure that the ground truth and predictions are aligned
    labels <- unname(groundTruth)
    predictions <- unname(predictions[names(groundTruth)])
    
    # Calculate the confusion matrix elements
    TP <- sum(predictions == 1 & labels == 1)
    TN <- sum(predictions == 0 & labels == 0)
    FP <- sum(predictions == 1 & labels == 0)
    FN <- sum(predictions == 0 & labels == 1)
    
    # Calculate performance metrics
    accuracy <- (TP + TN) / (TP + FP + TN + FN)
    precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    recall <- TP / (TP + FN)  # Recall is also known as sensitivity
    specificity <- TN / (TN + FP)
    f1_score <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
    TPR <- recall  # True Positive Rate is the same as recall
    FPR <- FP / (FP + TN)  # False Positive Rate
    
    # Compile the results
    results <- c(
        TP = round(TP),
        FP = round(FP),
        TN = round(TN),
        FN = round(FN),
        TPR = round(TPR,4),
        FPR = round(FPR,4),
        Accuracy = round(accuracy,4),
        Precision = round(precision,4),
        Recall = round(recall,4),
        Specificity = round(specificity,4),
        F1_Score = round(f1_score,4),
        RSQ = RSQ
    )
    
    return(results)
}
