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
rep_vec <- rep_vec[-3]
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
    #rep_i = "rep3"
    
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
        path1_name <-"CLP_pre_mDC"
        path2_name <-"CLP_pre_pDC"
        root_pp = c("Y_17")
        path1_pp <- c("Y_2", "Y_17", "Y_21", "Y_24", "Y_25", "Y_26", "Y_27", "Y_30", "Y_37", "Y_38", "Y_40")
        path2_pp = c("Y_4", "Y_8", "Y_9", "Y_10", "Y_11", "Y_13", "Y_17", "Y_20", "Y_28", "Y_29", "Y_31", "Y_41")
        age <- "28"
        sex <- "Female"
    } else if (rep_i == "rep3") {
        path1_name <-"HSC_LMPP"
        path2_name <-"HSC_EMP"
        root_pp = c("Y_62")
        path1_pp =  c('Y_4', 'Y_5', 'Y_6', 'Y_11', 'Y_18', 'Y_26', 'Y_29', 'Y_48', 'Y_53', 'Y_56', 'Y_62', 'Y_68', 'Y_79', 'Y_82', 'Y_84', 'Y_89', 'Y_91', 'Y_94', 'Y_98', 'Y_103', 'Y_106')
        path2_pp =  c('Y_2', 'Y_12', 'Y_14', 'Y_23', 'Y_25', 'Y_30', 'Y_32', 'Y_38', 'Y_42', 'Y_45', 'Y_55', 'Y_62', 'Y_65', 'Y_66', 'Y_67', 'Y_70', 'Y_72', 'Y_76', 'Y_78', 'Y_95', 'Y_96', 'Y_99', 'Y_109')
        individual <- "Donor-3"
        age <- "19"
        sex <- "Female"
    }
    
    ob <- m3_select_path(cdsObj = object.list[[rep_i]],
                         use_shiny = F, plot_purity = F,
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
