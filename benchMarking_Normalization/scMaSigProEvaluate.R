# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(maSigPro))
suppressPackageStartupMessages(library(gtools))

# Load Prefix
parentFix <- "benchMarking_Normalization/"
scriptFix <- paste0("old_Scripts/")
inPrefix <- paste0(parentFix, "data/output/scMaSigPro/")
outPrefix <- paste0(parentFix, "data/output/scMaSigPro/")

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
    source(paste0(scriptFix, "classify_sig_genes.R"))
}

# Load File names
all_file_names <- list.files(inPrefix)
all_file_names <- str_remove(all_file_names, ".RDS")
all_file_names <- str_split(all_file_names, "_")
names(all_file_names) <- list.files(inPrefix)
#all_file_names <- all_file_names[names(all_file_names) == "zi_shape_-0.2_TrueCounts_sum_tstep.RDS"]

# Some Variables
rsq <- 0.6

load()


## Run for Loop one-by-one
for (i in names(all_file_names)) {
    
    # Set Variables
    file_name <- i
    file_name_list <- all_file_names[[i]]
    
    # Read the file
    tstepObj <- readRDS(paste0(inPrefix, file_name))
    
    # Look for corresponding Ground truth
    gt_name <- paste0(parentFix, "data/output/GroundTruth/", paste(file_name_list[1],file_name_list[2], file_name_list[3], "GroundTruth.tsv",sep = "_"))
    groundTruth <- read.table(file = gt_name,
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
    
    # ROC
    roc <- ggplot(df.performance, aes(x = FPR, y = TPR))+
        geom_point(size = 2, alpha = 0.7, color = "red")+
        geom_path(size = 1, alpha = 0.7, color = "blue")+
        scale_x_continuous(breaks = seq(0,0.5, 0.05))+
        scale_y_continuous(breaks = seq(0.5,1, 0.05))+
        labs(title= "ROC-curve, Different Values of R-Square", 
             x = "False Positive Rate (1-Specificity)", 
             y = "True Positive Rate (Sensitivity)") + theme_classic()+
        geom_text(data = df.performance, aes(x = FPR, y = TPR, label = VARIABLE), 
                  nudge_y = 0.02, size = 4, nudge_x =0) +
        theme(panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
              panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
              axis.text = element_text(size = rel(1.5)),
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
        )+ geom_vline(xintercept = 0.01, colour= "green") + geom_vline(xintercept = 0.05, colour= "orange")  + geom_vline(xintercept = 0.1, colour= "red")
    
    # GGPlot2
    acc<- ggplot(df.performance, aes(x = VARIABLE, y = ACCURACY))+
        geom_point(size = 2, alpha = 0.7, color = "red")+
        scale_x_continuous(breaks = seq(0.1, 0.9, 0.05))+
        geom_path(size = 1, alpha = 0.7, color = "blue")+
        scale_y_continuous(breaks =  seq(0.7, 1, 0.05))+
        labs(title= "Accuracy Against Chnaging R Square",
             subtitle = "Red.Dots: False Negatives",
             x = "Increasing Value for R-Square", 
             y = "Accuracy") + theme_classic()+
        geom_text(data = df.performance, aes(x = VARIABLE, y = ACCURACY, label = FN), 
                  nudge_y= 0.01, nudge_x = 0.01, size = 3) +
        theme(panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
              panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
              axis.text = element_text(size = rel(1.5)),
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
        )
    
    roc_name <- paste0("benchMarking_Normalization/data/output/evaluation/singleROC/",
                       str_remove(i, "tstep.RDS"), "ROC.png")
    png(roc_name, width = 600, height = 600)
    print(roc)
    dev.off()
    
    acc_name <- roc_name <- paste0("benchMarking_Normalization/data/output/evaluation/singleACC/",
                                   str_remove(i, "tstep.RDS"), "ACC.png") 
    png(acc_name, width = 600, height = 600)
    print(acc)
    dev.off()
    
    measure_name <- roc_name <- paste0("benchMarking_Normalization/data/output/evaluation/singleTable/",
                                   str_remove(i, "tstep.RDS"), "Performance.tsv") 
    write.table(df.performance, measure_name, row.names = F, sep = "\t",
                quote = F)
}


# Load all the count tables in the 
ob <- readRDS("benchMarking_Normalization/data/output/SimObjectsSCE/zi_shape_-0.2.RDS")
clr <- as.matrix(read.table("benchMarking_Normalization/data/output/CLRNorm/zi_shape_-0.2_CLRNorm.tsv",
                  row.names = 1, header = T, sep = "\t"))
clr <- clr[, rownames(as.data.frame(colData(ob)))]
ob@assays@data$clr <- clr

fqnorm <- as.matrix(read.table("benchMarking_Normalization/data/output/FQNorm/zi_shape_-0.2_FQNorm.tsv",
                            row.names = 1, header = T, sep = "\t"))
fqnorm <- fqnorm[, rownames(as.data.frame(colData(ob)))]
ob@assays@data$fqnorm <- fqnorm

sct <- as.matrix(read.table("benchMarking_Normalization/data/output/SctNorm/zi_shape_-0.2_SctNorm.tsv",
                               row.names = 1, header = T, sep = "\t"))
sct <- sct[, rownames(as.data.frame(colData(ob)))]
ob@assays@data$sct <- sct

cpm <- as.matrix(read.table("benchMarking_Normalization/data/output/CPMNorm/zi_shape_-0.2_CPMNorm.tsv",
                            row.names = 1, header = T, sep = "\t"))
cpm <- cpm[, rownames(as.data.frame(colData(ob)))]
ob@assays@data$cpm <- cpm

logn <- as.matrix(read.table("benchMarking_Normalization/data/output/LogNorm/zi_shape_-0.2_LogNorm.tsv",
                            row.names = 1, header = T, sep = "\t"))
logn <- logn[, rownames(as.data.frame(colData(ob)))]
ob@assays@data$logn <- logn

rc <- as.matrix(read.table("benchMarking_Normalization/data/output/RcNorm/zi_shape_-0.2_RcNorm.tsv",
                             row.names = 1, header = T, sep = "\t"))
rc <- rc[, rownames(as.data.frame(colData(ob)))]
ob@assays@data$rc <- rc

plot_simulations(ob, assay_type = "TrueCounts", plot2d = T, plot3d = F, frame = 2)
plot_simulations(ob, assay_type = "counts", plot2d = T, plot3d = F, frame = 2 )
plot_simulations(ob, assay_type = "cpm", plot2d = T, plot3d = F, frame = 2)
plot_simulations(ob, assay_type = "clr", plot2d = T, plot3d = F, frame = 2)
plot_simulations(ob, assay_type = "logn", plot2d = T, plot3d = F, frame = 2)
plot_simulations(ob, assay_type = "rc", plot2d = T, plot3d = F, frame = 2)
plot_simulations(ob, assay_type = "fqnorm", plot2d = T, plot3d = F, frame = 2)
plot_simulations(ob, assay_type = "sct", plot2d = T, plot3d = F, frame = 2)

