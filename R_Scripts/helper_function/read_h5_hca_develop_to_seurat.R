# Read data from human cell atlas developmental H5 and Return as a seurat object

read_h5_hca_develop_to_seurat <- function(filename = "HTA08_v01_A06_Science_human_tcells.h5ad",
                                          filepath = "/home/priyansh/gitDockers/scMaSigPro_supp_data/Analysis_Public_Data_2",
                                          project_name = "HCA_Developmental",
                                          assay_name = "RNA",
                                          min_cells = 10,
                                          min_features = 10,
                                          add_gene_annotations = FALSE) {
  # Reuire Library
  suppressPackageStartupMessages(require(Seurat))
  suppressPackageStartupMessages(require(hdf5r))
  suppressPackageStartupMessages(require(rhdf5))
  suppressPackageStartupMessages(require(tidyverse))
  suppressPackageStartupMessages(require(biomaRt))
  suppressPackageStartupMessages(require(Matrix))
  #
  filename <- "HTA08_v01_A06_Science_human_tcells.h5ad"
  filepath <- "/supp_data/Analysis_Public_Data_2"
  project_name <- "HCA_Developmental_tcells"
  assay_name <- "RNA"
  min_cells <- 1000
  min_features <- 500

  # Create filename
  file_path <- paste(filepath, filename, sep = "/")

  # Load frame
  slot_frame <- data.frame(h5ls(file = file_path))

  # Load counts
  data <- h5read(file_path, paste0("/X/data"))
  indices <- h5read(file_path, paste0("/X/indices"))
  indptr <- h5read(file_path, paste0("/X/indptr"))

  # Load barcodes
  cell_metadata_encoded <- as.data.frame(h5read(file_path, paste0("/obs")))
  barcodes <- cell_metadata_encoded$index

  # Load feature Ids
  feature_ids <- as.data.frame(h5read(file_path, paste0("/var")))
  colnames(feature_ids) <- "gene_short_name"
  feature_ids[["feature_ids_symbol"]] <- make.unique(feature_ids[["gene_short_name"]])
  feature_ids$gene_short_name <- str_remove(feature_ids$gene_short_name, "\\..*")

  # Create Count Matrix
  count_matrix <- sparseMatrix(
    i = indices + 1, p = indptr, x = as.vector(data),
    repr = "C", dimnames = list(feature_ids$feature_ids_symbol, barcodes)
  )

  # Create Cell Metadata
  # Load the encoding
  cell_metadata_labels <- h5read(file_path, paste0("/uns"))

  # Create Cell Metadata
  cell_metadata <- data.frame(matrix(ncol = length(cell_metadata_labels), nrow = length(barcodes)))
  colnames(cell_metadata) <- str_remove(names(cell_metadata_labels), "_categories")

  # Loop over labels
  for (label in names(cell_metadata_labels)) {
    # get all values for the label
    label_values <- cell_metadata_labels[[label]]

    # Create factors with encoded metadata
    encodel_label_string <- str_remove(label, "_categories")
    encoded_labels <- factor(cell_metadata_encoded[[encodel_label_string]], levels = 1:length(label_values), labels = label_values)

    cell_metadata[[encodel_label_string]] <- encoded_labels
  }

  # Add Barcode Labels
  rownames(cell_metadata) <- barcodes

  # Set as columns
  cell_metadata <- cell_metadata %>%
    rownames_to_column(var = "barcode")

  # Replace any space in column names by "_" and make lower case
  colnames(cell_metadata) <- tolower(str_replace_all(colnames(cell_metadata), " ", "_"))

  # Convert all the columns except the "barcode" column as factors
  cell_metadata <- cell_metadata %>%
    mutate_if(is.character, as.factor)
  rownames(cell_metadata) <- cell_metadata$barcode

  # Load embeddings
  embeddings <- h5read(file_path, paste0("/obsm"),
    compoundAsDataFrame = FALSE
  )
  embeddings_2d <- embeddings$X_umap
  embeddings_2d <- t(embeddings_2d)
  embeddings_2d <- as.data.frame(embeddings_2d)
  rownames(embeddings_2d) <- barcodes

  # If Annotations are Requested
  if (add_gene_annotations) {
    # Connect to the Ensembl database
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    # Query the Ensembl database to get the desired information
    attributes <- c(
      "ensembl_gene_id",
      "external_gene_name",
      "description",
      "entrezgene_id",
      "gene_biotype"
    )

    # Create Query Without the Version
    query_vector <- unique(feature_ids$gene_short_name)

    # Call
    results <- getBM(
      attributes = attributes,
      filters = "hgnc_symbol",
      values = query_vector,
      mart = ensembl
    )

    # Merge the results back with the original data frame
    merged_data <- merge(feature_ids, results, by.x = "gene_short_name", by.y = "external_gene_name", all.x = TRUE)

    # Remove Duplicates if any
    merged_data <- merged_data %>% distinct()

    # Update
    feature_ids <- merged_data

    # Subset countmatrix by gene metadata ids
    count_matrix <- count_matrix[merged_data$feature_ids_symbol, ]
  } else {
    rownames(feature_ids) <- feature_ids$feature_ids_symbol
  }


  # Create Seurat Object
  assay_ob <- CreateAssayObject(count_matrix)
  rownames(feature_ids) <- as.vector(rownames(assay_ob))
  sob <- CreateSeuratObject(assay_ob,
    project = project_name,
    assay = assay_name,
    min.cells = min_cells, min.features = min_features,
    meta.data = cell_metadata
  )

  # Add Gene Level metatdata
  sob@assays$RNA@meta.features <- feature_ids

  # Return
  return(sob)
}
