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
inPrefix <- paste0(parentFix, "data/output/evaluation/singleTable/")
outPrefix <- paste0(parentFix, "data/output/ZI/")
dir.create(outPrefix,showWarnings = F)

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
all_file_names <- str_split(all_file_names, "_")
names(all_file_names) <- list.files(inPrefix)

# Some Variables
rsq <- 0.6

load()

performance.df <- data.frame()

## Run for Loop one-by-one
for (i in names(all_file_names)) {
    
    
    performance <- as.data.frame(read.table(paste0(inPrefix, i),
                                            header = T, sep = "\t"))
    
    performance$NORM_METHOD <- all_file_names[[i]][4]
    performance$BIN_METHOD <- all_file_names[[i]][5]
    performance$ZI_SHAPE <-  all_file_names[[i]][3]
    
    performance.df <- rbind(performance.df, performance)
    
}

# Checking Zero Inflation
zi.sum <- performance.df[performance.df$BIN_METHOD == "sum", ]
zi.sum.norm <- zi.sum[zi.sum$NORM_METHOD %in% c("RawCounts", "FQNorm"),]

zi.sum.norm[zi.sum.norm$ZI_SHAPE == -0.2, "INFLATION"] <- "42"
zi.sum.norm[zi.sum.norm$ZI_SHAPE == 0, "INFLATION"] <- "52"
zi.sum.norm[zi.sum.norm$ZI_SHAPE == 0.2, "INFLATION"] <- "62"
zi.sum.norm[zi.sum.norm$ZI_SHAPE == 1, "INFLATION"] <- "82"

ggplot(zi.sum.norm, aes(x = VARIABLE, y = ACCURACY,
                        linetype = NORM_METHOD,
                        color = INFLATION))+
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
          legend.position = "bottom",
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
    )

# ROC
ggplot(zi.sum.norm, aes(x = FPR, y = TPR, linetype = NORM_METHOD, color = INFLATION))+
    #geom_point()+
    geom_path(linewidth = 1.5, alpha = 0.6)+
    scale_x_continuous(breaks = seq(0,0.5, 0.05))+
    scale_y_continuous(breaks = seq(0.5,1, 0.05))+
    scale_color_brewer(palette = "Set1") +
    labs(title= "ROC-curve, Different Values of R-Square", 
         x = "False Positive Rate (1-Specificity)", 
         y = "True Positive Rate (Sensitivity)") + theme_classic()+
    theme(panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
          panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
          legend.position = "bottom",
          legend.key.size = unit(4, "cm"),
          legend.key.width = unit(2, "cm"), 
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 14),
          axis.text = element_text(size = rel(1.5)),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
    )+ geom_vline(xintercept = 0.01, colour= "lightgrey") + geom_vline(xintercept = 0.05, colour= "grey")+
    geom_vline(xintercept = 0.1, colour= "darkgrey") + guides(color = guide_legend(key_width = unit(3, "cm"), key_height = unit(4, "cm")))

