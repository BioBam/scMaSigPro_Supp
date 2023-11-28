##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

library(scMaSigPro)

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

# Get genes for Path2(EMP)

emp_up <- sc.get.features(scmpObj = scmp.obj,
                          query = "unique",
                          rsq = 0.7,
                          significant.intercept = "none",
                          vars = "groups",includeInflu = F,
                          union.ref = "Path2vsPath1",
                          union.target = "Path1",
                          union.ref.trend = "down",
                          union.target.trend = "up")
emp_up <- emp_up[order(emp_up)]
emp_up

GATA2 <- sc.PlotGroups(scmp.obj, "GATA2", logType = "log", logs = T,
              smoothness = 0.001, pseudoCount = 1, significant = F)

GATA2 <- GATA2 + theme_classic(base_size = 15) +theme( legend.box = "vertical",
               legend.direction = "vertical",
               legend.position = c(0.2, 0.2),
               axis.text.x = element_text(angle = 45, hjust = 1),
               panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
               panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),)
GATA2

saveRDS(GATA2, file = "Figures/MainArticle/MainArticle_FigureD.RDS")
