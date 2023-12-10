##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

set.seed(007)

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(scMaSigPro))

# Load object
object.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
    sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_cds.RDS"))
})


# Create ScMaSigPro Input
for (rep_i in c(rep_vec)){
    
    if (rep_i == "rep1") {
        path1_name <-"MPP_to_GMP"
        path2_name <-"MPP_to_CLP"
        root_pp = c("Y_26")
        path1_pp = c("Y_2", "Y_8", "Y_9", "Y_17", "Y_26", "Y_31", "Y_34", "Y_35", "Y_38", "Y_41", "Y_58", "Y_69", "Y_76", "Y_85", "Y_87", "Y_95", "Y_116", "Y_121", "Y_133", "Y_137", "Y_146", "Y_199", "Y_210")
        path2_pp = c("Y_6", "Y_14", "Y_26", "Y_47", "Y_52", "Y_54", "Y_55", "Y_78", "Y_107", "Y_117", "Y_134", "Y_147", "Y_148", "Y_178", "Y_184", "Y_203", "Y_218")
        individual <- "Donor-1"
        age <- "35"
        sex <- "Male"
    } else if (rep_i == "rep2") {
        path1_name <-"EMP_to_ProgEryth"
        path2_name <-"EMP_to_ProgMk"
        root_pp = c("Y_49")
        path1_pp <- c("Y_49", "Y_74", "Y_92", "Y_113", "Y_143", "Y_154", "Y_194")
        path2_pp =  c("Y_28", "Y_42", "Y_49", "Y_65", "Y_136")
        age <- "28"
        sex <- "Female"
    } else if (rep_i == "rep3") {
        path1_name <-"EMP_to_ProgMk"
        path2_name <-"EMP_to_ProgEryth"
        root_pp = c("Y_109")
        path1_pp = c("Y_6", "Y_27", "Y_30", "Y_56", "Y_72", "Y_82", "Y_86", "Y_90", "Y_104", "Y_108", "Y_109", "Y_113", "Y_129", "Y_139", "Y_142", "Y_162", "Y_176", "Y_200", "Y_204", "Y_224", "Y_248")
        path2_pp =  c("Y_1", "Y_44", "Y_48", "Y_58", "Y_76", "Y_109", "Y_114", "Y_147", "Y_148", "Y_170", "Y_183", "Y_186", "Y_203", "Y_206", "Y_234", "Y_236", "Y_241")
        individual <- "Donor-3"
        age <- "19"
        sex <- "Female"
    }
    
    ob <- m3_select_path(cdsObj = object.list[[rep_i]],
                         use_shiny = T, plot_purity = F,
                         annotation_col = "cell_type",
                         m3_pp = list(
                             root_pp = root_pp,
                             path1_pp = path1_pp,
                             path2_pp = path2_pp, 
                             path1_name =path1_name,
                             path2_name = path2_name)
    )
    
    saveRDS(
        ob,
        paste0(dirPath, rep_i, "/", "scMaSigPro_Input_.RDS")
    )
    
}
