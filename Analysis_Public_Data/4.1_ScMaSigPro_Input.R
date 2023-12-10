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
for (rep_i in c(rep_vec)[1]){
    
    if (rep_i == "rep1") {
        path1_name <-"EMP_EarlyErythrocyte"
        path2_name <-"EMP_ProgMk"
        root_pp = c("Y_39")
        path1_pp = c("Y_39", "Y_62", "Y_85", "Y_144", "Y_160", "Y_167", "Y_168")
        path2_pp = c("Y_2", "Y_3", "Y_39", "Y_43", "Y_97")
        individual <- "Donor-1"
        age <- "35"
        sex <- "Male"
    } else if (rep_i == "rep2") {
        path1_name <-"HSC_GMP"
        path2_name <-"HSC_CLP"
        root_pp = c("Y_58")
        path1_pp <- c("Y_5", "Y_14", "Y_19", "Y_20", "Y_23", "Y_40", "Y_43", "Y_49", "Y_51", "Y_53", "Y_58", "Y_65", "Y_86", "Y_89", "Y_97", "Y_101", "Y_114", "Y_127", "Y_147", "Y_156", "Y_164", "Y_176", "Y_183", "Y_185", "Y_190", "Y_195", "Y_205")
        path2_pp =  c("Y_4", "Y_11", "Y_22", "Y_58", "Y_77", "Y_107", "Y_117", "Y_121", "Y_168", "Y_171", "Y_192", "Y_203", "Y_207", "Y_216")
        age <- "28"
        sex <- "Female"
    } else if (rep_i == "rep3") {
        path1_name <-"MPP_ProB"
        path2_name <-"MPP_pmDC"
        root_pp = c("Y_20")
        path1_pp =  c("Y_8", "Y_13", "Y_20", "Y_34", "Y_35", "Y_37", "Y_50", "Y_71", "Y_76", "Y_89", "Y_95", "Y_113", "Y_115", "Y_119", "Y_127", "Y_130", "Y_133", "Y_153", "Y_160", "Y_198", "Y_202", "Y_205", "Y_218", "Y_223")
        path2_pp =  c("Y_4", "Y_10", "Y_18", "Y_20", "Y_29", "Y_46", "Y_47", "Y_48", "Y_51", "Y_54", "Y_61", "Y_80", "Y_84", "Y_107", "Y_116", "Y_141", "Y_150", "Y_163", "Y_169", "Y_175", "Y_181", "Y_187", "Y_188", "Y_217", "Y_219", "Y_233")
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
