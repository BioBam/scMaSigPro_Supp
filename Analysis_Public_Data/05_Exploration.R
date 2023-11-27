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

# Select the Significant genes
scmp.obj <- sc.get.siggenes(scmpObj = scmp.obj,
                            rsq = 0.7,
                            Q = 0.05,
                            vars = "groups",
                            significant.intercept = "dummy",
                            term.Q = 0.05,
                            includeInflu = T)
    

sc.PlotGroups(scmp.obj, "GATA1", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "MPO", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "EPOR", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "ELANE", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "GATA2", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "IRF2", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "AZU1", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "CD34", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "ASXL1", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "SERPINB10", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "RNASE3", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "MS4A3", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "PRTN3", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "CTSG", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "RNASE2", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "NPW", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

# GMP
sc.PlotGroups(scmp.obj, "MYCT1", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "PBX1", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "IGSF10", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1, significant = F)

sc.PlotGroups(scmp.obj, "CYTL1", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1, significant = T)


sc.PlotGroups(scmp.obj, "HPGDS", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "AVP", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "NPR3", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)

sc.PlotGroups(scmp.obj, "CRHBP", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1)


sc.PlotGroups(scmp.obj, "CCL18", logType = "log", logs = T,
              smoothness = 1, pseudoCount = 1, significant = F)

