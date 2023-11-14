##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: scMaSigPro Application ########################
##########################################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))

# Load scMaSigPro
suppressPackageStartupMessages(library(scMaSigPro))

# Load object
object.list <- lapply(c("don1", "don2", "don3"), function(don){
    readRDS(paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/scMaSigPro_Processed_",don, ".RDS"))
})
names(object.list) <- c("don1", "don2", "don3")

# Extract Gene with high r2
# Run scMaSigPro
scMaSigPro.list.gene <- lapply(object.list, function(don){
    scmp.obj <- sc.get.siggenes(scmpObj = don, rsq = 0.7, vars = "groups")
    sols <- showSol(scmp.obj, view = F, return = T, includeInflu = F)
    sols[is.na(sols$`p-value`), "p-value"] <- 0
    sols <- sols[sols$`p-value` <= 0.05, ]
    sols <- sols[sols$`R-squared` >= 0.7,]
    return(rownames(sols))
})

don1.plots <- list()
don2.plots <- list()
don3.plots <- list()

# Plot Gene Expression
for (i in scMaSigPro.list.gene$don1){
    don1.plots[[i]] <- sc.PlotGroups(scmpObj = object.list$don1,feature_id = i,
                                    logs = T,
                                    logType = "log")
}
for (i in scMaSigPro.list.gene$don2){
    don2.plots[[i]] <- sc.PlotGroups(scmpObj = object.list$don2,feature_id = i,
                                     logs = T,
                                     logType = "log")
}
for (i in scMaSigPro.list.gene$don3){
    don3.plots[[i]] <- sc.PlotGroups(scmpObj = object.list$don3,feature_id = i,
                                     logs = T,
                                     logType = "log")
}

# Donor-3
don1.plots$MPO # MPO Goes up in Path1- Monocyte and Path2 is Dendritic
don2.plots$AZU1.2 # AZU1 Goes up in path2 which is dendritic, path 1 is Myeloid Monocyte
don3.plots$SCT.1 # SCT goes up in path2 is dentric predc, path 1 is monocyte like

# Compettive Expression for Myeloid bias
#Cdk6 downregulation in myeloid progenitors releases Runx1 from Cdk6 inhibition, thereby allowing terminal differentiation.
don3.plots$CDK6 # Path2 goes down for densritic
don3.plots$RUNX2 # Goes up in dendritic

don3.plots$CDK6 / don3.plots$RUNX2

# 
don3.plots$LMO2 # Expression goes down during lineage commitment



mpo <- sc.PlotGroups(scmpObj = object.list$don1,feature_id = "MPO",
              logs = T,smoothness = 0.001,
              logType = "log")


mpo <- mpo + ggtitle("Myeloperoxidase (MPO) Expression",
              subtitle = "R-Square: 0.926") + 
    theme_classic(base_size = 15) +
    theme(legend.box = "vertical",
          legend.direction = "vertical",
          legend.title=element_text(size=13), 
          legend.text=element_text(size=9),
          panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
          panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
          #legend.position = "bottom"
          legend.position = c(0.8, 0.5), legend.justification = c("left", "top")
    )
    

ggsave(plot = mpo,
       filename = "Figure1_C.png",
       path = "Article_Image/",
       dpi = 1200, width = 5, height = 5)
saveRDS(
    mpo, "Article_Image/Figure1_C.RDS"
)
