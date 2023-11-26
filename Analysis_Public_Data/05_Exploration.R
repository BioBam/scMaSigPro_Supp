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

# Extract Sol
sol <- showSol(scmp.obj, includeInflu = T)
    
# Change NA 
sol[is.na(sol$`p-value`), "p-value"] <- 1
sol[is.na(sol$`R-squared`), "R-squared"] <- 0
sol[,!(colnames(sol) %in% c("p-value", "R-squared"))][is.na(sol[,!(colnames(sol) %in% c("p-value", "R-squared"))])] <- 1

# Get genes with pvalue as significant
sol <- sol[sol$`p-value` <= 0.05, ]
sol <- sol[sol$`R-squared` >= 0.6, ]

# Remove the columns
sol <- sol[,!(colnames(sol) %in% c("p-value", "R-squared"))]

# Get the ones for which only beta changes
beta0.sol <- sol[sol$p.valor_beta0 <= 0.05, ]
beta0.sol.not <- sol[sol$p.valor_beta0 >= 0.05, ]

sc.PlotGroups(scmp.obj,
              sample(rownames(beta0.sol.not), 1),
              logs = T, logType = "log10")
