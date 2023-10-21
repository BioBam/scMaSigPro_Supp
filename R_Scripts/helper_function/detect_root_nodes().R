find_root_pp <- function(cds.object, cell = "HSC", cell_col = "predicted.celltype.l2", top_nodes = 5){
    
    # Retrieve the cell metadata from the given cds object
    cell_meta <- as.data.frame(colData(cds.object))
    
    # Subset the metadata to only include the column specified by 'cell_col'
    cell_meta <- cell_meta[, cell_col, drop = F]
    
    # Add a 'barcode' column using the rownames of cell_meta
    cell_meta$barcode <- rownames(cell_meta)
    
    # Extract the principal graph from the cds object
    ppgraph <- principal_graph_aux(cds.object)
    
    # Get the closest vertex data from the principal graph
    vertex <- as.data.frame(ppgraph@listData$UMAP$pr_graph_cell_proj_closest_vertex)
    
    # Rename the column to "pp"
    colnames(vertex) <- "pp"
    
    # Add a prefix "Y_" to the 'pp' column values
    vertex$pp <- paste("Y", vertex$pp, sep = "_")
    
    # Add a 'barcode' column using the rownames of vertex
    vertex$barcode <- rownames(vertex)
    
    # Merge the cell metadata and the vertex dataframes based on the 'barcode' column
    vertex <- merge(cell_meta, vertex, "barcode")
    
    # Filter rows where the cell type matches the specified 'cell' value (default is "HSC")
    filtered_df <- subset(vertex, predicted.celltype.l2 == cell)
    
    # Count the occurrences of each unique 'pp' value
    pp_counts <- table(filtered_df$pp)
    
    # Sort the 'pp' counts in descending order
    sorted_pp_counts <- sort(pp_counts, decreasing = TRUE)
    
    # Get the top 'pp' values based on the specified 'top_nodes' value (default is 5)
    top_pp <- names(sorted_pp_counts)[1:top_nodes]
    
    # Return the top 'pp' values
    return(top_pp) 
}
