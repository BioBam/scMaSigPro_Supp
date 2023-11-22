# Evaluation with iCobra
suppressPackageStartupMessages(library(iCOBRA))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(RColorConesa))

# Set Path
dirPath <- "/supp_data/benchmarks/04_ComparisonWithTradeSeq/simulated/sce/"
resPath <- "/supp_data/benchmarks/04_ComparisonWithTradeSeq/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Read TradeSeq Data
load(paste0(resPath, "CobraInputObject.RData"))

# Convert all values to numeric
cobraInput <- as.data.frame(apply(cobra.dataset, 2, FUN = function(x) {
  as.numeric(x)
}))
rownames(cobraInput) <- rownames(cobra.dataset)

# Read Ground truth
load(paste0(dirPath, "testTradeSeq.RData"))

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

#COBRAapp(cob_data)

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
ROC <- ROC+ theme_classic(base_size = 15) +
    scale_color_manual(
        labels = c("scMaSigPro",
                   "tradeSeq diffEnd()", 
                   #"TS earlyDETest()",
                   "tradeSeq pattern()"),
                   values = c(colorConesa(7)[1],
                              colorConesa(7)[3],
                              #colorConesa(7)[5],
                              colorConesa(7)[7])
                   ,
    )+
    labs(
        title = "ROC Curve: Comparison with TradeSeq",
        subtitle = "R-Square threshold: 0.6",
        y = "True Positive Rate (Sensitivity)",
        x = "False Positive Rate (1-Specificity)",
        color = "Method"
    ) +theme(
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
        #legend.position = "bottom"
        legend.position = c(0.5, 0.5), legend.justification = c("left", "top")) +
    geom_vline(xintercept = 0.01, colour = "lightgrey", linetype = "dotted") +  # Highlighted the x-intercept of 0.01
    geom_vline(xintercept = 0.05, colour = "lightgrey",linetype = "dotted") +
    geom_vline(xintercept = 0.1, colour = "lightgrey",linetype = "dotted")  +
    scale_y_continuous(breaks = unique(c(seq(0.8, 1, 0.05), seq(0, 1, 0.1))))+
scale_x_continuous(breaks = unique(c(c(0.05, 0.01), seq(0.1, 1, 0.1))))

print(ROC)

# ggsave(plot = ROC,
#        path = "Figures/MainArticle",
#        dpi = 1000,  filename = "Figure1_B.png",
#        width = 6, height = 5)
# saveRDS(ROC, file = "Figures/MainArticle/Figure1_B.RDS")

# TPRvsFDR <- plot_fdrtprcurve(cobraplot, title = "TPR vs FDR")
# TPR <- plot_tpr(cobraplot, title = "TPR")
# FPR <- plot_fpr(cobraplot, title = "FPR")
# 
# ggsave(filename = paste0(resPath, "ROC.png"), plot = ROC, dpi = 600, width = 7, height = 5)
# ggsave(filename = paste0(resPath, "TPRvsFDR.png"), plot = TPRvsFDR, dpi = 600, width = 7, height = 5)
# ggsave(filename = paste0(resPath, "TPR.png"), plot = TPR, dpi = 600, width = 7, height = 5)
# ggsave(filename = paste0(resPath, "FPR.png"), plot = FPR, dpi = 600, width = 7, height = 5)
