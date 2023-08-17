# Evaluation with iCobra
suppressPackageStartupMessages(library(iCOBRA))

# Set path
inPath <- "supp/05_Comparison_with_TradeSeq/data/output/"

# Read TradeSeq Data
CobraInput <- readRDS(paste0(inPath,"CobraInputObject.RDS"))

# Convert all values to numeric
Input <- as.data.frame(apply(CobraInput, 2, FUN = function(x){as.numeric(x)}))
rownames(Input) <- rownames(CobraInput) 

# Read Ground truth
sce.obj <- readRDS(paste0("supp/05_Comparison_with_TradeSeq/data/input/PseudobulkInput/zi_60_mid_0_shape_0.25.RDS"))
gt <- sce.obj$gt
gt <- as.data.frame(as.matrix(gt[mixedorder(rownames(gt)),, drop = F]))
gt <- gt[,4,drop=F]
colnames(gt) <- "status"
gtInput <- apply(gt, 2, FUN = function(x){as.integer(x)})
rownames(gtInput) <- rownames(gt)
gtInput <- as.data.frame(gtInput)

# Create Cobra data
cob_data <- COBRAData(pval = CobraInput,
                      truth = gtInput,
                      score = gt[,-1, drop = F])

# Calculate Adjusted p-value
cobradata_custom <- calculate_adjp(cob_data)

# Calculate Performance
cobraperf <- calculate_performance(cobradata_custom,
                                   binary_truth = "status",
                                   splv = "none",
                                   maxsplit = 2)
# Calculate for plots
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", 
                                   facetted = TRUE)

a <- plot_roc(cobraplot, title = "ROC")
b <- plot_fdrtprcurve(cobraplot, title = "TPR vs FDR")
c <- plot_tpr(cobraplot, title = "TPR")
d <- plot_fpr(cobraplot, title = "FPR")
ggpubr::ggarrange(a,b,c,d)
