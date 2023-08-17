make_bulk_counts <- function(counts, cluster_member_col =  "cluster.members", round = T,
                             bin_col = "bin", pseudo_bulk_profile){
    
    
    meta.info <- pseudo_bulk_profile[, c(cluster_member_col, bin_col)]
    
    # initialize an empty dataframes
    pseudo_bulk_counts_median <- pseudo_bulk_counts_sum <- pseudo_bulk_counts_mean <- data.frame(
        matrix(nrow = nrow(counts), ncol = 0))
    
    # One row at a time
    for (i in c(1:nrow(meta.info))){
        
        # Select Rows One by One
        selRow <- meta.info[i,, drop = F]
        
        # Save the column Names
        col_name <- selRow[,bin_col, drop = F]
        
        # cluster components
        clus_com <- c(str_split(selRow[,cluster_member_col], "\\|"))[[1]]
        
        # Select columns from the frame
        sel_cols_df <- counts[, colnames(counts) %in% clus_com, drop = F]
        
        # Add Mean
        sel_cols_df_mean <- as.data.frame(rowMeans(sel_cols_df))
        colnames(sel_cols_df_mean) <- col_name
        pseudo_bulk_counts_mean <- cbind(pseudo_bulk_counts_mean, sel_cols_df_mean)
        
        # Add Sum
        sel_cols_df_sum <- as.data.frame(rowSums(sel_cols_df))
        colnames(sel_cols_df_sum) <- col_name
        pseudo_bulk_counts_sum <- cbind(pseudo_bulk_counts_sum, sel_cols_df_sum)
        
        # Add median
        sel_cols_df_median <- as.data.frame(rowMedians(sel_cols_df))
        colnames(sel_cols_df_median) <- col_name
        pseudo_bulk_counts_median <- cbind(pseudo_bulk_counts_median, sel_cols_df_median)
    }
    
    if (round == T){
        pseudo_bulk_counts_mean <- round(pseudo_bulk_counts_mean)
        pseudo_bulk_counts_sum <- round(pseudo_bulk_counts_sum)
        pseudo_bulk_counts_median <- round(pseudo_bulk_counts_median)
    }
    
    # return
    return(list(pseudo_bulk_counts_mean = pseudo_bulk_counts_mean,
                pseudo_bulk_counts_sum = pseudo_bulk_counts_sum,
                pseudo_bulk_counts_median = pseudo_bulk_counts_median))
    
}