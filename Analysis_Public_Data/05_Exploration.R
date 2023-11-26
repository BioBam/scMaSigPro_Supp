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
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Load
scMaSigPro.list <- lapply(rep_vec, function(don) {
    ob <- readRDS(
        paste0("/supp_data/Analysis_Public_Data/scMaSigPro_Processed_", don, ".RDS")
    )
})

# Explore donor-1
scmp.obj <- scMaSigPro.list$rep1
