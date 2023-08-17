# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(maSigPro))
suppressPackageStartupMessages(library(gtools))

# Load Prefix
parentFix <- "benchmarking_Offset/"
scriptFix <- paste0("old_Scripts/")
inPrefix <- paste0(parentFix, "data/input/PB_Input/")
outPrefix <- paste0(parentFix, "data/output/OffsetAnalysis/")

# Load Custim Functions
load <- function() {
    source(paste0(scriptFix, "plot_simulations().R"))
    source(paste0(scriptFix, "add_gene_anno().R"))
    source(paste0(scriptFix, "calc_bin_size.R"))
    source(paste0(scriptFix, "calcNormCounts.R"))
    source(paste0(scriptFix, "entropy_discretize.R"))
    source(paste0(scriptFix, "make_bulk_design.R"))
    source(paste0(scriptFix, "make_bulk_counts.R"))
    source(paste0(scriptFix, "create_range.R"))
    source(paste0("old_Scripts/classify_sig_genes.R"))
    source("Testing_And_Development/old_maSigPro/p.vector_test.R")
}
load()

# 42: zi_shape_-0.2_
# 52: zi_shape_0_
# 62: zi_shape_0.2_
# 87: zi_shape_1_

# Read data list
data.list = readRDS("benchmarking_Offset/data/input/PB_Input/zi_shape_1_RawCounts.RDS")

# Extract the design file
exp.design.file <- data.list[["cell_meta"]]

# Create List of all the counts
count_list <- list(
    sum = data.list[["count_table"]][["sum"]],
    mean = data.list[["count_table"]][["mean"]],
    median = data.list[["count_table"]][["median"]]
)

# Get the counts
ct <- count_list[["sum"]]

if (all(colnames(ct) == rownames(exp.design.file)) == F) { stop("MisMatch") }
        
# Create MaSigPro Design File
design <- make.design.matrix(
    edesign = as.data.frame(exp.design.file), # bulkMeta is edesign
    degree = poly.degree,
    time.col = 1, repl.col = 2
    )
                
# Run P-Vector
p.vector.fit <- p.vector_test(
    data = ct, design = design, Q = 0.05,
    MT.adjust = "BH", min.obs = min.gene, counts = TRUE,
    theta = theta.val, epsilon = ep
    )
                
# Runnig T Step
tstep.fit <- T.fit_test(p.vector.fit,
                        step.method = "backward",
                        family = p.vector.fit$family,
                        )
tstepObj <- tstep.fit

# Look for corresponding Ground truth
groundTruth <- read.table(file = "benchmarking_Offset/data/input/GroundTruth/zi_shape_1_GroundTruth.tsv",
                          header = T, sep = "\t", row.names = 1
)

# Checking For R-Square
gene.change <- rownames(groundTruth[groundTruth$cobraTruth == 1,])
gene.no.change <- rownames(groundTruth[groundTruth$cobraTruth == 0,])

# Create Dataframe to plot Performance
df.performance <- data.frame(VARIABLE = 0, TP = 0, FP = 0, TN = 0, FN = 0,
                             TPR = 0, FPR = 0, FNR = 0, TNR = 0, 
                             INFLU = 0, INFLU_FP = 0, ACCURACY = 0)

# Run For all Values of R square
for (j in seq(0.05, 0.95, 0.05)){
    
    # Classify Genes
    classifications <- classify_sig_genes(masigpro.tstep = tstepObj, 
                                          rsq.value = j, draw.venn = F,
                                          all_genes = rownames(groundTruth), 
                                          num_group = 2)
    
    # True Positive
    tp <- intersect(gene.change, classifications$detected.genes) 
    
    # True Negative
    tn <- intersect(gene.no.change, classifications$undetected.gene)
    
    # False Positive
    fp <- intersect(gene.no.change, classifications$detected.genes) 
    
    # False Negative
    fn <- intersect(gene.change, classifications$undetected.gene)
    
    # Influential Genes
    influ <- colnames(tstepObj$influ.info)
    
    # Influential + False Positive
    fp_influential <- intersect(fp,colnames(tstepObj$influ.info)) 
    
    # Sensitivity, True Positive Rate / Recall
    senstivity <- length(tp)/(length(tp) + length(fn))
    
    # False Negative Rate
    fnr <- length(fn)/(length(tp) + length(fn))
    
    # Specificity / True Negative Rate
    specificity <- length(tn)/(length(tn) + length(fp))
    
    # False Positive Rate
    fpr <- 1- specificity
    
    # Accuracy
    accuracy <- (length(tp)+length(tn))/(length(tp)+length(tn)+length(fp)+length(fn))
    
    # Make the vector of results
    res.perfom <- data.frame(VARIABLE = j,
                             TP = length(tp),
                             FP = length(fp),
                             TN = length(tn),
                             FN = length(fn),
                             TPR = round(senstivity,3),
                             FPR = round(fpr,3),
                             FNR = round(fnr,3),
                             TNR = round(specificity,3),
                             INFLU = length(influ),
                             INFLU_FP = length(fp_influential),
                             ACCURACY = round(accuracy, 3)
    )
    # Add vector to the frame
    df.performance <- rbind(df.performance, res.perfom)
    gc()
}

# Remove the first Row
df.performance <- df.performance[-1,]

measure_name <- paste0("benchmarking_Offset/data/output/OffsetAnalysis/Performance_82_Offset.tsv")
write.table(df.performance, measure_name, row.names = F, sep = "\t",
            quote = F)

df_82_off <- read.table("benchmarking_Offset/data/output/OffsetAnalysis/Performance_82_Offset.tsv", sep = "\t", header = T)
df_62_off <- read.table("benchmarking_Offset/data/output/OffsetAnalysis/Performance_62_Offset.tsv", sep = "\t", header = T)
df_52_off <- read.table("benchmarking_Offset/data/output/OffsetAnalysis/Performance_52_Offset.tsv", sep = "\t", header = T)
df_42_off <- read.table("benchmarking_Offset/data/output/OffsetAnalysis/Performance_42_Offset.tsv", sep = "\t", header = T)


# 42: zi_shape_-0.2_
# 52: zi_shape_0_
# 62: zi_shape_0.2_
# 87: zi_shape_1_

df_42 <- read.table("benchMarking_Normalization/data/output/evaluation/singleTable/zi_shape_-0.2_RawCounts_sum_Performance.tsv", sep = "\t", header = T)
df_52 <- read.table("benchMarking_Normalization/data/output/evaluation/singleTable/zi_shape_0_RawCounts_sum_Performance.tsv", sep = "\t", header = T)
df_62 <- read.table("benchMarking_Normalization/data/output/evaluation/singleTable/zi_shape_0.2_RawCounts_sum_Performance.tsv", sep = "\t", header = T)
df_82 <- read.table("benchMarking_Normalization/data/output/evaluation/singleTable/zi_shape_1_RawCounts_sum_Performance.tsv", sep = "\t", header = T)


df_42_off$Offset <- "Yes"
df_52_off$Offset <- "Yes"
df_62_off$Offset <- "Yes"
df_82_off$Offset <- "Yes"

df_42$Offset <- "No"
df_52$Offset <- "No"
df_62$Offset <- "No"
df_82$Offset <- "No"

df_42_off$ZI <- "42"
df_52_off$ZI <- "52"
df_62_off$ZI <- "62"
df_82_off$ZI <- "82"

df_42$ZI <- "42"
df_52$ZI <- "52"
df_62$ZI <- "62"
df_82$ZI <- "82"

df.performance <- rbind(df_42_off,
      df_52_off,
      df_62_off,
      df_82_off,
      df_42,
      df_52,
      df_62,
      df_82)

# ROC
a <- ggplot(df.performance, aes(x = FPR, y = TPR, linetype = Offset,
                                color = ZI))+ 
    geom_point()+
    geom_path(linewidth = 1.5, alpha = 0.6)+
    scale_x_continuous(breaks = seq(0,0.5, 0.05))+
    scale_y_continuous(breaks = seq(0.5,1, 0.05))+
    scale_color_brewer(palette = "Set1") +
    labs(title= "ROC-curve, Different Values of R-Square",
         x = "False Positive Rate (1-Specificity)",
         y = "True Positive Rate (Sensitivity)") + theme_classic()+
    theme(panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
          panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
          legend.position = "none",
          legend.key.size = unit(4, "cm"),
          legend.key.width = unit(2, "cm"),
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 14),
          axis.text = element_text(size = rel(1.5)),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
    )+ geom_vline(xintercept = 0.01, colour= "lightgrey") + geom_vline(xintercept = 0.05, colour= "grey")+
    geom_vline(xintercept = 0.1, colour= "darkgrey") + guides(color = guide_legend(key_width = unit(3, "cm"), key_height = unit(4, "cm")))


b<- ggplot(df.performance, aes(x = VARIABLE, y = ACCURACY,
                               linetype = Offset, color = ZI))+
    geom_point()+
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.05))+
    geom_path(linewidth = 1.5, alpha = 0.6)+
    scale_y_continuous(breaks =  seq(0.7, 1, 0.05))+
    scale_color_brewer(palette = "Set1") +
    labs(title= "Accuracy Against Changing R Square",
         subtitle = "Red.Dots: False Negatives",
         x = "Increasing Value for R-Square",
         y = "Accuracy") + theme_classic()+
    theme(panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
          panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
          axis.text = element_text(size = rel(1.5)),
          legend.key.size = unit(4, "cm"),
          legend.key.width = unit(2, "cm"),
          legend.position = "left",
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
    )

ggpubr::ggarrange(a,b, ncol=2)
