# Data Downloaded from https://www.ebi.ac.uk/gxa/sc/experiments/E-HCAD-6/downloads
# Date: Nov 14, 2023

# Set paths
dataInput = "Analysis_Public_Data/data/SingleCellExperimentAtlas/E-HCAD-6-quantification-raw-files/"

# Load libs
library(Seurat)
library(tidyverse)
library(biomaRt)

# Load Raw Counts
rawMat = ReadMtx(mtx = paste0(dataInput, "E-HCAD-6.aggregated_filtered_counts.mtx"),
                 cells = paste0(dataInput, "E-HCAD-6.aggregated_filtered_counts.mtx_cols"),
                 features = paste0(dataInput, "E-HCAD-6.aggregated_filtered_counts.mtx_rows"))

# Extract cells
cell_meta_temp = data.frame(
    barcode = str_split_i(colnames(rawMat), pattern = "-", i = 2),
    group = str_split_i(colnames(rawMat), pattern = "-", i = 1)
)

# Load Metainformation
cell_metadata = read.table(
    "Analysis_Public_Data/data/SingleCellExperimentAtlas/HaematopoieticProfiling-10x_cell_type_2020-03-12.csv",
    sep = ",", header = T)

# Merge
cell_metadata = merge(cell_metadata, cell_meta_temp, by = "barcode")

# Extract gene information
gene_metadata = data.frame(
    ensembl_gene_id  = rownames(rawMat),
    row.names = rownames(rawMat)
)

# Set ensg
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# Get symbols
ensembl_to_gensymbol <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'), 
                              filters = 'ensembl_gene_id', 
                              values = gene_metadata$ensembl_gene_id, 
                              mart = ensembl, uniqueRows = F)

# Convert Missing to ENGIDs
ensembl_to_gensymbol[ensembl_to_gensymbol$hgnc_symbol == "","hgnc_symbol"] = ensembl_to_gensymbol[ensembl_to_gensymbol$hgnc_symbol == "","ensembl_gene_id"]

# Set columns
colnames(ensembl_to_gensymbol) <- c("gene_short_name", "ensembl_gene_id")

# Mapp annotated IDs
gene_metadata <- merge(gene_metadata, ensembl_to_gensymbol, by = "ensembl_gene_id")

# Check for duplicates
length(gene_metadata$ensembl_gene_id) - length(unique(gene_metadata$ensembl_gene_id))

# Drop 3 duplicates
gene_metadata <- gene_metadata[!duplicated(gene_metadata$ensembl_gene_id),]

# Function to add suffix to duplicates
add_suffix <- function(x) {
    return(paste0(x, ".", ave(seq_along(x), x, FUN = seq_along)))
}

# Apply the function only to duplicated ensembl_gene_id values
gene_metadata$gene_short_name_unique <- ifelse(
    duplicated(gene_metadata$gene_short_name) | duplicated(gene_metadata$gene_short_name, fromLast = TRUE), 
                                    add_suffix(gene_metadata$gene_short_name), 
    gene_metadata$gene_short_name)


# Set rownames
rownames(gene_metadata) <- gene_metadata$gene_short_name_unique

# Remove mismatch from the counts
rawMat = rawMat[gene_metadata$ensembl_gene_id,]

# reset rownames
rownames(rawMat) = gene_metadata$gene_short_name_unique

# Abbreviate Cell Names
cell_metadata[cell_metadata$annotated_cell_identity.text == "dendtritic progenitor", "cell_type"] <- "pDend" 
cell_metadata[cell_metadata$annotated_cell_identity.text == "hematopoietic stem cell", "cell_type"] <- "HSC" 
cell_metadata[cell_metadata$annotated_cell_identity.text == "erythroid progenitor", "cell_type"] <- "pEryth" 
cell_metadata[cell_metadata$annotated_cell_identity.text == "myeloid progenitor", "cell_type"] <- "pMyel" 
cell_metadata[cell_metadata$annotated_cell_identity.text == "monocyte progenitor", "cell_type"] <- "pMono" 
cell_metadata[cell_metadata$annotated_cell_identity.text == "megakaryocyte progenitor", "cell_type"] <- "pMega" 
cell_metadata[cell_metadata$annotated_cell_identity.text == "common lymphoid progenitor", "cell_type"] <- "CLP" 
cell_metadata[cell_metadata$annotated_cell_identity.text == "hematopoietic multipotent progenitor", "cell_type"] <- "HMP" 

# Create group-barcode
cell_metadata$groupBarcode <- paste(cell_metadata$group, cell_metadata$barcode, sep = "-")
cell_metadata <- cell_metadata[!duplicated(cell_metadata$groupBarcode), ]
rownames(cell_metadata) <- cell_metadata$groupBarcode

# Subset counts and Create seurat object
rawMat <- rawMat[, colnames(rawMat) %in% rownames(cell_metadata)]

# Create Seurat Object
sob <- CreateSeuratObject(counts = rawMat, min.cells = 3, min.features = 200)
# Add metainfo
sob@meta.data <- cell_metadata[colnames(sob),]

# Basic QC
sob <- NormalizeData(sob)
sob <- FindVariableFeatures(sob)
sob <- ScaleData(sob)
sob <- RunPCA(sob, features = VariableFeatures(object = sob))
sob <- FindNeighbors(sob, dims = 1:10)
sob <- FindClusters(sob, resolution = 0.5)
sob <- RunUMAP(sob, dims = 1:10)
DimPlot(sob, reduction = "umap", group.by = "annotated_cell_identity.text")

# Seprate Groups
sob1 <- subset(sob, subset = cell_suspension.biomaterial_core.biomaterial_id %in% grep("P1", x = sob@meta.data$cell_suspension.biomaterial_core.biomaterial_id, value = T))
sob2 <- subset(sob, subset = cell_suspension.biomaterial_core.biomaterial_id %in% grep("P2", x = sob@meta.data$cell_suspension.biomaterial_core.biomaterial_id, value = T))
sob3 <- subset(sob, subset = cell_suspension.biomaterial_core.biomaterial_id %in% grep("P3", x = sob@meta.data$cell_suspension.biomaterial_core.biomaterial_id, value = T))

# Sumsample and Recompute 
sob.list <- parallel::mclapply(list(sob1, sob2, sob3), function(s){
    s <- subset(s, annotated_cell_identity.text %in% c("hematopoietic multipotent progenitor",
                                                       "myeloid progenitor",
                                                       "dendtritic progenitor",
                                                       "monocyte progenitor"))
    s <- NormalizeData(s)
    s <- FindVariableFeatures(s)
    s <- ScaleData(s)
    s <- RunPCA(s, features = VariableFeatures(object = s))
    s <- FindNeighbors(s, dims = 1:10)
    s <- FindClusters(s, resolution = 0.5)
    s <- RunUMAP(s, dims = 1:10)
}, mc.cores = 24)

names(sob.list) <- c("don1", "don2", "don3")

# Extract and Save Data
lapply(names(sob.list), function(i){
    
    j <- sob.list[[i]]
    print(j)
    
    # Extract Cell Metadata
    cell.meta <- j@meta.data
    
    # Counts
    counts <- j@assays$RNA@counts
    
    # SaveR
    saveRDS(cell.meta, paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/", i,".meta.RDS"))
    saveRDS(counts, paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/", i,".counts.RDS"))
    return(NULL)
})
