# Evaluation with iCobra
suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pROC))

# Set paths
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
rdsPath <- paste0(base_string, "comparison/sim/")
figPath <- paste0(base_string, "figures/")
outPath <- paste0(base_string, "comparison/out/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Create Directory if does not exist
dir.create(figPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_hd, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_lr, showWarnings = FALSE, recursive = TRUE)
dir.create(tabPath, showWarnings = FALSE, recursive = TRUE)
dir.create(rdsPath, showWarnings = FALSE, recursive = TRUE)
dir.create(outPath, showWarnings = FALSE, recursive = TRUE)

# Read TradeSeq Data
load(paste0(outPath, "CobraInputObject.RData"))

# Convert all values to numeric
cobraInput <- as.data.frame(apply(cobra.dataset, 2, FUN = function(x) {
  as.numeric(x)
}))
rownames(cobraInput) <- rownames(cobra.dataset)

# Read Ground truth
load(paste0(rdsPath, "testTradeSeq.RData"))

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

# Create Cobra data
cob_data <- COBRAData(
  pval = cobraInput,
  truth = gtInput
)

# COBRAapp(cob_data)

# Assuming cobraInput contains continuous probabilities
predictions <- cbind(cobraInput, gtInput)

# Binarize Predictions
predictions$scmp_0.6 <- ifelse(predictions$scmp_0.6 <= 0.05, 1, 0)
predictions$TS_diffEnd <- ifelse(predictions$TS_diffEnd <= 0.05, 1, 0)
predictions$TS_pattern <- ifelse(predictions$TS_pattern <= 0.05, 1, 0)

# AUC
auc_scMaSigPro <- auc(predictions$status, predictions$scmp_0.6)
auc_tradeSeq_diffEnd <- auc(predictions$status, predictions$TS_diffEnd)
auc_tradeSeq_pattern <- auc(predictions$status, predictions$TS_pattern)

# Calculate Adjusted p-value
cobradata_custom <- calculate_adjp(cob_data)

# Calculate Performance
cobraperf <- calculate_performance(cobradata_custom,
  binary_truth = "status",
  splv = "none",
  maxsplit = 3
)
# Calculate for plots
cobraplot <- prepare_data_for_plot(cobraperf,
  colorscheme = "Dark2",
  facetted = TRUE
)

ROC <- plot_roc(cobraplot, title = "ROC")
ROC <- ROC + theme_classic(base_size = 14) +
  scale_color_manual(
    labels = c(
      paste("scMaSigPro |", round(auc_scMaSigPro, 3)),
      paste("tradeSeq diffEnd() |", round(auc_tradeSeq_diffEnd, 3)),
      paste("tradeSeq pattern() |", round(auc_tradeSeq_pattern, 3))
    ),
    values = c(
      colorConesa(7)[1],
      colorConesa(7)[3],
      # colorConesa(7)[5],
      colorConesa(7)[7]
    ),
  ) +
  labs(
    title = "ROC Curve: Comparison with TradeSeq",
    subtitle = "R-Square threshold: 0.6",
    y = "True Positive Rate (Sensitivity)",
    x = "False Positive Rate (1-Specificity)",
    color = "Method | AUC"
  ) + theme(
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
    # legend.position = "bottom"
    legend.position = c(0.5, 0.5), legend.justification = c("left", "top")
  ) +
  geom_vline(xintercept = 0.01, colour = "lightgrey", linetype = "dotted") + # Highlighted the x-intercept of 0.01
  geom_vline(xintercept = 0.05, colour = "lightgrey", linetype = "dotted") +
  geom_vline(xintercept = 0.1, colour = "lightgrey", linetype = "dotted") +
  scale_y_continuous(breaks = unique(c(seq(0.8, 1, 0.05), seq(0.5, 1, 0.1))), limits = c(0.5, 1)) +
  scale_x_continuous(breaks = unique(c(c(0.05, 0.01, 0.1), seq(0.2, 1, 0.1))), limits = c(0, 1))

print(ROC)

saveRDS(ROC, file = paste0(outPath, "MainFigure_2C.rds"))
#
# # Save Truth table with results
# df <- cbind(gt, cobraInput)
# saveRDS(df, file = paste0(outPath, "Prediction.RDS"))
#
# plot_overlap(cobraplot)
#
#
# # Define Cutom plot for the iCOBRA
# data <- cbind(cobraInput, gtInput)
#
# # Load the libraries
# library(pROC)
# library(ggplot2)
#
# # Calculate the ROC curve and AUC for each p-value column
# roc_TS_pattern <- roc(data$status, data$TS_pattern)
# roc_TS_diffEnd <- roc(data$status, data$TS_diffEnd)
# roc_scmp_0.6 <- roc(data$status, data$scmp_0.6)
#
# # Get AUC values
# auc_TS_pattern <- auc(roc_TS_pattern)
# auc_TS_diffEnd <- auc(roc_TS_diffEnd)
# auc_scmp_0.6 <- auc(roc_scmp_0.6)
#
# # Convert ROC data to data frames for ggplot
# roc_TS_pattern_df <- data.frame(
#   fpr = 1 - roc_TS_pattern$specificities,
#   tpr = roc_TS_pattern$sensitivities,
#   method = paste("TS_pattern (AUC =", round(auc_TS_pattern, 2), ")")
# )
#
# roc_TS_diffEnd_df <- data.frame(
#   fpr = 1 - roc_TS_diffEnd$specificities,
#   tpr = roc_TS_diffEnd$sensitivities,
#   method = paste("TS_diffEnd (AUC =", round(auc_TS_diffEnd, 2), ")")
# )
#
# roc_scmp_0.6_df <- data.frame(
#   fpr = 1 - roc_scmp_0.6$specificities,
#   tpr = roc_scmp_0.6$sensitivities,
#   method = paste("scmp_0.6 (AUC =", round(auc_scmp_0.6, 2), ")")
# )
#
# # Combine the data frames
# roc_data <- rbind(roc_TS_pattern_df, roc_TS_diffEnd_df, roc_scmp_0.6_df)
#
# # Plot the ROC curves with AUC values
# ggplot(roc_data, aes(x = fpr, y = tpr, color = method)) +
#   geom_line() +
#   geom_abline(linetype = "dashed") +
#   labs(
#     title = "ROC Curve with AUC",
#     x = "False Positive Rate",
#     y = "True Positive Rate"
#   ) +
#   theme_minimal() +
#   theme(legend.title = element_blank())
