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
        paste0("/supp_data/Analysis_Public_Data/",don,"/","scMaSigPro_Processed_", don, ".RDS")
    )
})

# Explore donor-1
scmp.obj <- scMaSigPro.list$rep2

# Get genes for Path2(EMP)
emp_up <- sc.get.features(scmpObj = scmp.obj,
                          query = "union",
                          rsq = 0.7,
                          unique.group = "Path2vsPath1",
                          significant.intercept = "dummy",
                          vars = "groups",includeInflu = T,
                          union.ref = "Path1",
                          union.target = "Path2vsPath1",
                          union.ref.trend = "stable",
                          union.target.trend = "any")
emp_up <- emp_up[order(emp_up)]

az_emp <- c("MYCT1","CRHBP","NPR3","AVP","GATA2","HPGDS","CYTL1","CRYGD","IGSF10","PBX1")
az_gmp <- c("SERPINB10","RNASE3","MS4A3","PRTN3","ELANE","AZU1","CTSG","RNASE2","RETN","NPW")
emp <- emp_up[emp_up %in% az_emp]
gmp <- emp_up[emp_up %in% az_gmp]

print(paste("For EMP", length(emp), "out of" ,length(az_emp), "i.e.", paste(emp_up[emp_up %in% az_emp], collapse = ","),
            "|-------| For GMP", length(gmp), "out of" ,length(az_gmp)))

stop()




write.table(emp_up, "test.tsv", sep = "\t", quote = F, row.names = F)
ITGA2B
plotTrend(scmp.obj, "BANK1", logType = "log", logs = T,
                       smoothness = 0.001, pseudoCount = 1, significant = F)


# 
# GATA2 <- sc.PlotGroups(scmp.obj, "GATA2", logType = "log", logs = T,
#               smoothness = 0.001, pseudoCount = 1, significant = F)
# 
# GATA2 <- GATA2 + theme_classic(base_size = 15) +theme( legend.box = "vertical",
#                legend.direction = "vertical",
#                legend.position = c(0.2, 0.2),
#                axis.text.x = element_text(angle = 45, hjust = 1),
#                panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
#                panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),)
# GATA2
# 
# saveRDS(GATA2, file = "Figures/MainArticle/MainArticle_FigureD.RDS")
