# Evaluation with iCobra
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SingleCellExperiment))

# Set Path
dirPath <- "benchmarks/04_ComparisonWithTradeSeq/data/input/sce/"
resPath <- "benchmarks/04_ComparisonWithTradeSeq/data/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Read TradeSeq Data
load(paste0(resPath, "CobraInputObject.RData"))

# Convert all values to numeric
cobraInput <- as.data.frame(apply(cobra.dataset, 2, FUN = function(x) {
    as.numeric(x)
}))
rownames(cobraInput) <- rownames(cobra.dataset)

# Read Ground truth
load(paste0(dirPath, "Test_TradeSeq.RData"))

# Extract Gene counts
row_data <- as.data.frame(rowData(sim.sce))

# Extract status
gt <- row_data[, "status", drop = F]

# Add status
gt$status <- ifelse(gt$status == "No_Change", 0, 1)

# Reformat
gtInput <- apply(gt, 2, FUN = function(x) {
    as.integer(x)
})
rownames(gtInput) <- rownames(gt)
gtInput <- as.data.frame(gtInput)
library(ROCR)
library(ggplot2)

generate_metrics <- function(data, pval_column, ground_truth){
    # Ensure the data is sorted the same way as ground_truth
    data <- data[order(row.names(data)), ]
    data$status <- ground_truth$status
    
    data <- data[!is.na(data[[pval_column]]), ]
    
    # Convert p-values to binary classification
    thresholds <- seq(0, 1, by = 0.05)
    recall_vals <- list()
    precision_vals <- list()
    FPR_vals <- list()
    TPR_vals <- list()
    
    for (i in 1:length(thresholds)) {
        threshold <- thresholds[i]
        data$prediction <- ifelse(data[[pval_column]] < threshold, 1, 0)
        
        # Calculate metrics
        # We use 1-pvalue as predictions, since higher values of prediction 
        # should indicate more confidence in a positive classification.
        pred <- prediction(1 - data[[pval_column]], data$status)
        
        perf_prec <- performance(pred, measure = "prec", x.measure = "rec")
        
        perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
        
        recall_vals[i] <- perf_prec@x.values
        precision_vals[i] <- perf_prec@y.values
        FPR_vals[i] <- perf_prec@x.values[[1]]
        TPR_vals[i] <- perf_prec@y.values[[1]]
    }
    
    
    View(precision_vals)
    stop()
    
    # Plot Precision-Recall
    pr_df <- data.frame(Recall = recall_vals, Precision = precision_vals)
    
    View(pr_df)
    stop()
    
    pr_plot <- ggplot(pr_df, aes(x = Recall, y = Precision)) +
        geom_line() +
        labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision")
    
    print(pr_plot)
    
    # Plot ROC
    roc_df <- data.frame(FPR = FPR_vals, TPR = TPR_vals)
    
}
# Assuming your data is stored in a dataframe named df and ground truth in gt
pval_columns <- c("TS_pattern", "TS_early", "TS_diffEnd", "scmp_0.6")
plots <- lapply(pval_columns, generate_metrics, data = as.data.frame(cobraInput), ground_truth = gt)



