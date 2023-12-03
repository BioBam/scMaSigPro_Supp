##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

go_enrichment <- function(gene_list, scmp.obj, rep, age,sex, path, 
                          ont = "BP", pAdjustMethod = "fdr", nterms = 5, sig.level = 0.05){
    
    # Load
    suppressPackageStartupMessages(library(clusterProfiler))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    suppressPackageStartupMessages(library(enrichplot))
    suppressPackageStartupMessages(library(AnnotationDbi))
    
    # Separate Ensembl IDs and gene symbols
    ensembl_ids <- gene_list[grep("^ENSG", gene_list)]
    gene_symbols <- gene_list[!((gene_list) %in% ensembl_ids)]
    
    # Convert Ensembl IDs to Entrez IDs
    entrez_ids_from_ensembl <- mapIds(org.Hs.eg.db,
                                      keys = ensembl_ids,
                                      column = "ENTREZID",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")
    valid_entrez_ids_ensg <- entrez_ids_from_ensembl[!is.na(entrez_ids_from_ensembl)]
    
    # Convert Symbol IDs to Entrez IDs
    entrez_ids_from_symbol <- mapIds(org.Hs.eg.db,
                                     keys = gene_symbols,
                                     column = "ENTREZID",
                                     keytype = "SYMBOL",
                                     multiVals = "first")
    valid_entrez_ids_symb <- entrez_ids_from_symbol[!is.na(entrez_ids_from_symbol)]
    
    # Combine
    entrez_id <- c(valid_entrez_ids_symb, valid_entrez_ids_ensg)

    
    # Separate Ensembl IDs and gene symbols
    ensembl_ids <- gene_list[grep("^ENSG", gene_list)]
    gene_symbols <- gene_list[!((gene_list) %in% ensembl_ids)]
    
    # Set up universe
    universe_vector <- rownames(scmp.obj@compress.sce@assays@data@listData$bulk.counts)
    universe_vector_ensg <- universe_vector[grep("^ENSG", universe_vector)]
    universe_vector_sym <- universe_vector[!((universe_vector) %in% universe_vector_ensg)]
    
    # Convert Ensembl IDs to Entrez IDs
    universe_vector_ensg_entrez <- mapIds(org.Hs.eg.db,
                                          keys = universe_vector_ensg,
                                          column = "ENTREZID",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")
    valid_universe_vector_ensg_entrez <- universe_vector_ensg_entrez[!is.na(universe_vector_ensg_entrez)]
    
    # Convert Symbol IDs to Entrez IDs
    universe_vector_sym_entrez <- mapIds(org.Hs.eg.db,
                                         keys = universe_vector_sym,
                                         column = "ENTREZID",
                                         keytype = "SYMBOL",
                                         multiVals = "first")
    valid_universe_vector_sym_entrez <- universe_vector_sym_entrez[!is.na(universe_vector_sym_entrez)]
    
    # Create universe
    background_vector <- c(valid_universe_vector_sym_entrez, valid_universe_vector_ensg_entrez)
    
    # Run GO enrichment
    ego <- enrichGO(gene = entrez_id,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    universe = background_vector,
                    ont = ont,  # Biological Processes
                    pAdjustMethod = pAdjustMethod,  # Benjamini-Hochberg adjustment
                    qvalueCutoff = sig.level,  # Set threshold for q-value
                    readable = TRUE)  # To show readable gene names
    
    # Check
    if(dim(ego)[1] >= 1){
        
        # Calculate ema
        ego_similarity <- pairwise_termsim(ego)
        
        return(list(
            ego = ego,
            dot = dotplot(ego, showCategory=nterms) + ggtitle(subtitle = paste("Donor:", rep, "| Age:", age,
                                                                 "| Sex:", sex),
                                                           paste("Path: ", path, 
                                                                            "| Ontology: ",
                                                                            ont)),
            ema = emapplot(ego_similarity)
        ))
        
    }else{
        return(NULL)
    }
}
