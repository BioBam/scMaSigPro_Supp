# Data Downloaded from https://www.ebi.ac.uk/gxa/sc/experiments/E-HCAD-6/downloads

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

# Separate per donor
cellMetaDataD1 <- cell_metadata[cell_metadata$cell_suspension.biomaterial_core.biomaterial_id %in% 
                                    grep("P1",
                                        x = cell_metadata$cell_suspension.biomaterial_core.biomaterial_id,
                                        ignore.case = T, value = T),]
cellMetaDataD2 <- cell_metadata[cell_metadata$cell_suspension.biomaterial_core.biomaterial_id %in% 
                                    grep("P2",
                                         x = cell_metadata$cell_suspension.biomaterial_core.biomaterial_id,
                                         ignore.case = T, value = T),]
cellMetaDataD3 <- cell_metadata[cell_metadata$cell_suspension.biomaterial_core.biomaterial_id %in% 
                                    grep("P3",
                                         x = cell_metadata$cell_suspension.biomaterial_core.biomaterial_id,
                                         ignore.case = T, value = T),]

# Separate counts
rawMatD1 = rawMat[, paste(cellMetaDataD1$group, cellMetaDataD1$barcode, sep = "-")]
rawMatD2 = rawMat[, paste(cellMetaDataD2$group, cellMetaDataD2$barcode, sep = "-")]
rawMatD3 = rawMat[, paste(cellMetaDataD3$group, cellMetaDataD3$barcode, sep = "-")]

# Save the counts and day
saveRDS(rawMatD1, file = "Analysis_Public_Data/data/SingleCellExperimentAtlas/Rep1RawCounts.RDS")
saveRDS(rawMatD2, file = "Analysis_Public_Data/data/SingleCellExperimentAtlas/Rep2RawCounts.RDS")
saveRDS(rawMatD3, file = "Analysis_Public_Data/data/SingleCellExperimentAtlas/Rep3RawCounts.RDS")

# Save Metadata
saveRDS(gene_metadata, file = "Analysis_Public_Data/data/SingleCellExperimentAtlas/Gene_MetaData.RDS")

# Save Cell Metadata
saveRDS(cellMetaDataD1,
        file = "Analysis_Public_Data/data/SingleCellExperimentAtlas/Cell_MetaData_D1.RDS")
saveRDS(cellMetaDataD2,
        file = "Analysis_Public_Data/data/SingleCellExperimentAtlas/Cell_MetaData_D2.RDS")
saveRDS(cellMetaDataD3,
        file = "Analysis_Public_Data/data/SingleCellExperimentAtlas/Cell_MetaData_D3.RDS")