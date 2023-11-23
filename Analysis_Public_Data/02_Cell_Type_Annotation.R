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

# Prefix
prefixIn <- "Analysis_Public_Data/data"
prefixOut <- "Analysis_Public_Data/data"

# Get file names
rep_vec <- list.dirs(prefixIn, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "Setty_et_al_2019_Integrated_sob.h5seurat", "Human_Cell_Atlas"))]
names(rep_vec) <- rep_vec

# Run lapply
umap.list <- lapply(rep_vec, function(rep_i, inPath = prefixIn, outPath = prefixOut){
    
    # Step-1: Add Annotation for donors
    if (rep_i == "rep1") {
        individual <- "Donor-1"
        age <- "35"
        sex <- "Male"
    } else if (rep_i == "rep2") {
        individual <- "Donor-2"
        age <- "28"
        sex <- "Female"
    } else if (rep_i == "rep3") {
        individual <- "Donor-3"
        age <- "19"
        sex <- "Female"
    }
    
    # Step-2: Load the Seurat Object
    sob.prs <- LoadH5Seurat(
        paste0(paste(inPath, rep_i,sep = "/"), "/",rep_i, "_prs.h5seurat"),
        verbose = FALSE
    )
    
    # Step-3: Azimuth Annotations
    sob.bm <- RunAzimuth(
        query = sob.prs,
        reference = paste(inPath,
            "Azimuth_Human_BoneMarrow",sep = "/"
        )
    )
    
    # Plot the Annotations
    p <- DimPlot(sob.bm, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
        ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
        scale_color_hue(l = 50) + theme(legend.position = "bottom")
    
    ggsave(p,
           filename = paste(outPath, rep_i, paste0(rep_i, "_All_Anno_Azimuth.png"), sep = "/"),
           dpi = 1400, limitsize = FALSE, width = 8, height = 8
    )
    
    # Subset the Scores
    #sob.bm <- subset(sob.bm, subset = predicted.celltype.l2.score >= 0.6)

    # Plot the Annotations
    p <- DimPlot(sob.bm, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
        ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
        scale_color_hue(l = 50) + theme(legend.position = "bottom")

    # Subset the Major Cell States representing all lineage
    sob.bm <- subset(sob.bm,
                      subset = predicted.celltype.l2 %in% c(
                          "GMP", "LMPP","HSC",
                          "CLP", "EMP"
                      )
    )

    # Plot the Annotations
    umap <- DimPlot(sob.bm, reduction = "umap", pt.size = 1, group.by = "predicted.celltype.l2") +
        ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
        scale_color_hue(l = 50) + theme(legend.position = "bottom")

    # Save
    ggsave(umap,
           filename = paste(outPath, rep_i, paste0(rep_i, "_Final_Cells.png"), sep = "/"),
           dpi = 1400, limitsize = FALSE, width = 8, height = 8
    )
    
    # Save 
    file_name <- paste(outPath, rep_i, paste(rep_i, "azimuth", sep = "_"), sep = "/")
    SaveH5Seurat(
        object = sob.bm, filename = file_name,
        overwrite = T, verbose = FALSE
    )
})


# Selected cells
# A good representation of the major lineages and stages of differentiation. Here's a brief overview of what each cell type represents and how they fit into the lineage flow:
#     
# 1. **HSC (Hematopoietic Stem Cell)**: These are the multipotent stem cells that give rise to all blood cell types. They are at the top of the lineage hierarchy.
# 
# 2. **LMPP (Lymphoid-primed Multipotent Progenitor)**: These cells are derived from HSCs and have the potential to give rise to lymphoid and myeloid lineages but are more biased towards lymphoid cells.
# 
# 3. **GMP (Granulocyte-Macrophage Progenitor)**: These cells are also derived from HSCs and are biased towards forming granulocytes and macrophages.
# 
# 4. **CLP (Common Lymphoid Progenitor)**: These are progenitors for the lymphoid lineage, leading to cells like T cells, B cells, and NK cells.
# 
# 5. **EMP (Erythrocyte-Megakaryocyte Progenitor)**: These cells are specialized for generating erythrocytes (red blood cells) and megakaryocytes (which produce platelets).
# 
# 6. **Early Eryth (Early Erythrocytes)**: These are early-stage red blood cells still undergoing the process of maturation.
# 
# 7. **Late Eryth (Late Erythrocytes)**: These are more mature forms of red blood cells, closer to their final, circulating form.
#
# 8. **pre B: Precursor B cells, an early stage in B cell development.
#
# 9. **pre-mDC: Precursor myeloid Dendritic Cells.
#
# 10. **pre-pDC: Precursor plasmacytoid Dendritic Cells.
#
# 11. **pro B: Pro-B cells, a developmental stage in B cell maturation.
# 
# Creating a hematopoietic tree from the cell types provided:
#     
#     **Hematopoietic Stem Cell (HSC)**
#     - **Marker Genes**: CD34, CD133
# - **Source**: Blood and bone marrow.
# 
# - **Lymphoid Lineage**:
#     - **LMPP (Lymphoid-primed Multipotent Progenitor)**
#     
#     - **CLP (Common Lymphoid Progenitor)**
#     - **Marker Genes**: IL7R
# - **Source**: Bone marrow.
# 
# - **NK (Natural Killer cells)**
#     - **Marker Genes**: NCAM1 (CD56)
# - **Source**: Blood and bone marrow.
# 
# - **pro B**
#     - **Marker Genes**: CD19, VPREB1
# - **Source**: Bone marrow.
# 
# - **pre B**
#     - **Marker Genes**: CD79A, IGHM
# - **Source**: Bone marrow.
# 
# - **transitional B**
#     - **Marker Genes**: CD24, CD38
# - **Source**: Blood.
# 
# - **Naive B**
#     - **Marker Genes**: CD19, CD20
# - **Source**: Blood and lymphoid tissues.
# 
# - **CD4 Naive**
#     - **Marker Genes**: CD4, CCR7, SELL (CD62L)
# - **Source**: Blood and lymphoid tissues.
# 
# - **CD8 Memory**
#     - **Marker Genes**: CD8, CD45RO
# - **Source**: Blood and lymphoid tissues.
# 
# - **Myeloid Lineage**:
#     - **GMP (Granulocyte-Macrophage Progenitor)**
#     
#     - **Macrophage**
#     - **Marker Genes**: CD68, CD14
# - **Source**: Tissues and blood.
# 
# - **CD14 Mono (Monocytes)**
#     - **Marker Genes**: CD14, LY86
# - **Source**: Blood.
# 
# - **cDC2 (Conventional Dendritic Cell type 2)**
#     - **Marker Genes**: CD11c, CD1c
# - **Source**: Blood and tissues.
# 
# - **pDC (Plasmacytoid Dendritic Cells)**
#     - **Marker Genes**: CD123, BDCA-2
# - **Source**: Blood and lymphoid tissues.
# 
# - **pre-mDC (Precursor myeloid Dendritic Cells)**
#     - Likely similar markers to mature mDCs but in a precursor state.
# 
# - **EMP (Erythrocyte-Megakaryocyte Progenitor)**
#     
#     - **Platelet**
#     - **Marker Genes**: CD61 (ITGB3), CD41 (ITGA2B)
# - **Source**: Blood.
# 
# - **Early Eryth (Early Erythrocytes)**
#     - **Marker Genes**: GYPA, HBE1
# - **Source**: Bone marrow.
# 
# - **Late Eryth (Late Erythrocytes)**
#     - **Marker Genes**: HBB, HBA1
# - **Source**: Blood.
# 
# **Others (Not directly in the hematopoietic tree but related to hematopoiesis)**:
#     - **Stromal**
#     - **Marker Genes**: VIM, FAP
# - **Source**: Bone marrow and other tissues.
# 
# - **ASDC**: Not a standard abbreviation.
# - **BaEoMa**: Represents a combination lineage; markers would depend on the specific cell type in focus (Basophil, Eosinophil, or Mast cell).
# - **pre-pDC**: Likely precursor markers to mature pDCs.
# - **Prog Mk**: Progenitor Megakaryocytes, might have markers similar to mature megakaryocytes but in a precursor state.
