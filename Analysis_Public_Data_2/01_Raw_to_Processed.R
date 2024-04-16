#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Pre-Process ################
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))

# Prefix
dirPath <- "/home/priyansh/gitDockers/scMaSigPro_supp_data/Analysis_Public_Data_2"

# Load Custom loader for scanpy h5 (Form Human Developmental Atlas)
source("R_Scripts/helper_function/read_h5_hca_develop_to_seurat.R")

# Load Object
sob.raw <- read_h5_hca_develop_to_seurat(filename = "HTA08_v01_A06_Science_human_tcells.h5ad",
                                     filepath = dirPath,
                                     add_gene_annotations = TRUE,
                                     project_name = "t_cell_hca_developmental",
                                     assay_name = "RNA",
                                     min_cells = 100,
                                     min_features = 100)

# Calculate Percentage of Mitochondrial Reads
sob.raw[["percent.mt"]] <- PercentageFeatureSet(sob.raw, pattern = "^MT-")

# Plot QC
violins.raw.nFeature_RNA <- VlnPlot(sob.raw, features = "nFeature_RNA") + theme(legend.position = "none") + xlab("") +
    ggtitle(paste("HTA08_V01_A06_SCIENCE_HUMAN_TCELLS"), subtitle = "nFeature_RNA")
violins.raw.nCount_RNA <- VlnPlot(sob.raw, features = "nCount_RNA") + theme(legend.position = "none") + xlab("") +
    ggtitle(paste("HTA08_V01_A06_SCIENCE_HUMAN_TCELLS"), subtitle = "nCount_RNA")
violins.raw.percent.mt <- VlnPlot(sob.raw, features = "percent.mt") + theme(legend.position = "none") + xlab("") +
    ggtitle(paste("HTA08_V01_A06_SCIENCE_HUMAN_TCELLS"), subtitle = "percent.mt")
violins.raw <- ggarrange(violins.raw.nFeature_RNA, violins.raw.nCount_RNA,
                         violins.raw.percent.mt, labels = c("A.", "B.", "C."),
                         ncol = 3, nrow = 1)

# As data is already processed and log transforme we will not perform any filtering or Normalization

# We will remove cells without any cell-type annotataions
keep_cells <- sob.raw@meta.data[!is.na(sob.raw@meta.data$cell_types),"barcode"] 

# Subset
sob_prs <- subset(sob.raw, barcode %in% keep_cells)
sob.raw <-NULL

# Identify Variable Feature (Atleast 12000)
sob_prs <- FindVariableFeatures(
    sob_prs, 
    selection.method = "vst", 
    nfeatures = 12000)

# Scale Data
sob_prs <- ScaleData(sob_prs)

# Compute PCA with Variable Features
sob_prs <- RunPCA(sob_prs, features = VariableFeatures(object = sob_prs),
                  npcs = 200, seed.use = 123)

# Save Seurat Object as an H5
SaveH5Seurat(
    sob_prs, filename = "input_seurat_object", filepath = dirPath,
    overwrite = TRUE, verbose = TRUE)
