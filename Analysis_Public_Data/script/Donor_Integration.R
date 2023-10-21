##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: Azimuth Annotation + Monocle3 Input Prepare ###
##########################################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(ggpubr))

# Prefix
prefixIn <- "Analysis_Public_Data/data"
prefixOut <- "Analysis_Public_Data/data"

# Get file names
rep_vec <- list.dirs(prefixIn, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "Setty_et_al_2019_Integrated_sob.h5seurat",  "Human_Cell_Atlas"))]
names(rep_vec) <- rep_vec

# Run lapply
sob.list <- lapply(rep_vec, function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
    # Step-1: Add Annotation for donors
    if (rep_i == "rep1") {
        individual <- "Donor_1"
        age <- "35"
        sex <- "Male"
    } else if (rep_i == "rep2") {
        individual <- "Donor_2"
        age <- "28"
        sex <- "Female"
    } else if (rep_i == "rep3") {
        individual <- "Donor_3"
        age <- "19"
        sex <- "Female"
    }
    
    # Step-2: Load the Seurat Object
    sob.azimuth <- LoadH5Seurat(
        paste0(paste(inPath, rep_i,sep = "/"), "/",rep_i, "_azimuth.h5seurat"),
        verbose = FALSE
    )
    
    # Add Columns to Metadata
    sob.azimuth@meta.data$individual <- individual
    sob.azimuth@meta.data$age <- age
    sob.azimuth@meta.data$sex <- sex
    
    # Add Cell types
    sob.azimuth@meta.data$azimuth_celltype <- sob.azimuth@meta.data$predicted.celltype.l2
    
    # return
    return(sob.azimuth)
})

# Select Integration Features
features <- SelectIntegrationFeatures(object.list = sob.list, nfeatures = 4000)

# Define Anchor set
anchors <- FindIntegrationAnchors(
    object.list = sob.list,
    anchor.features = features
)

# Integrate the datasets
settyDonors <- IntegrateData(anchorset = anchors)

# Set Deafult assay
DefaultAssay(settyDonors) <- "integrated"

# Run the standard workflow for visualization and clustering
settyDonors <- ScaleData(settyDonors, verbose = FALSE)
settyDonors <- RunPCA(settyDonors, npcs = 50, verbose = FALSE)
settyDonors <- RunUMAP(settyDonors, reduction = "pca", dims = 1:50)

# Plot UMAP
umap_A <- DimPlot(settyDonors, reduction = "umap", pt.size = 1, group.by = "orig.ident") +
    ggtitle("Setty et al, 2019 (Integrated Donors)") +
    scale_color_hue(l = 50) + theme(legend.position = "bottom")

# Plot UMAP
umap_Cell <- DimPlot(settyDonors, reduction = "umap", pt.size = 1, group.by = "azimuth_celltype") +
    ggtitle("Setty et al, 2019 (Integrated Donors)") +
    scale_color_hue(l = 50) + theme(legend.position = "bottom")

# Plot UMAP
umap_Sex <- DimPlot(settyDonors, reduction = "umap", pt.size = 1, group.by = "sex") +
    ggtitle("Setty et al, 2019 (Integrated Donors)") +
    scale_color_hue(l = 50) + theme(legend.position = "bottom")

# Plot UMAP
umap_Age <- DimPlot(settyDonors, reduction = "umap", pt.size = 1, group.by = "age") +
    ggtitle("Setty et al, 2019 (Integrated Donors)") +
    scale_color_hue(l = 50) + theme(legend.position = "bottom")

# Plot Integrated data
integrated_umap <- ggarrange(umap_A, umap_Age, umap_Cell, umap_Sex,
                             labels = c("A. Integrated Data", "B. Age", "C. Cell", "D. Sex"))


# Write Seurat H5
file_name <- paste(prefixOut, "Setty_et_al_2019_Integrated_sob", sep = "/")
SaveH5Seurat(
    object = settyDonors, filename = file_name,
    overwrite = T, verbose = FALSE
)

# Save
ggsave(integrated_umap,
       filename = paste(prefixOut, "Setty_et_al_2019_Integrated_UMAP.png", sep = "/"),
       dpi = 1400, limitsize = FALSE, width = 12, height = 12
)
