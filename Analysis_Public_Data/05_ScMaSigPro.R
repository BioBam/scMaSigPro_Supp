##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

set.seed(007)

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(scMaSigPro))

# Load object
object.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_cds.RDS"))
})

# Extract path and create scMaSigpro Object
scmp.object.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
    
    if (rep_i == "rep1") {
        path1 <-"Early Eryth"
        path2 <-"Prog Mk"
        root = "EMP"
        individual <- "Donor-1"
        age <- "35"
        sex <- "Male"
    } else if (rep_i == "rep2") {
        path1 <-"GMP"
        path2 <-"CLP"
        individual <- "Donor-2"
        root = "HSC"
        age <- "28"
        sex <- "Female"
    } else if (rep_i == "rep3") {
        path1 <- "CLP"
        path2 <- "GMP"
        root = "LMPP"
        individual <- "Donor-3"
        age <- "19"
        sex <- "Female"
    }
    
  # get object
  rep_i_obj <- object.list[[rep_i]]
  
  # Create SCMP
  normCounts <- assay(rep_i_obj)
  cell.meta <- rep_i_obj@colData %>% as.data.frame()
  cell.meta$Pseudotime <- pseudotime(rep_i_obj)
  cell.meta <- cell.meta[cell.meta$cell_type %in% c(path1, path2, root),]
  vertex.df <- data.frame(
      vertex = paste("Y", rep_i_obj@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex, sep = "_"),
      row.names = rownames(rep_i_obj@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex),
      cell = rownames(rep_i_obj@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex)
  )
  
  # Attach Principal points
  pp <- data.frame(pp = paste("pp",
                              names(as.data.frame(rep_i_obj@principal_graph_aux$UMAP$dp_mst)), sep = "_"), 
                   vertex = names(as.data.frame(rep_i_obj@principal_graph_aux$UMAP$dp_mst)))
  
  # Merge
  vertex.df <- merge(vertex.df, pp, by =  "vertex")
  rownames(vertex.df) <- vertex.df$cell
  vertex.df <- vertex.df[rownames(cell.meta),, drop = F]
  cell.meta$vertex <- vertex.df$vertex
  
  # Set Cells
  cell.meta[cell.meta$cell_type == path1, "Lin"] <- paste0(root, "_", str_remove(path1, pattern = " "))
  cell.meta[cell.meta$cell_type == path2, "Lin"] <- paste0(root, "_", str_remove(path2, pattern = " "))
  
  # Divide EMPs
  path1.vertex <- cell.meta[cell.meta$cell_type == path1, "vertex"]
  cell.meta[cell.meta$cell_type == root & cell.meta$vertex %in% path1.vertex, "Lin"] <- paste0(root, "_", str_remove(path1, pattern = " "))
  
  path2.vertex <- cell.meta[cell.meta$cell_type == path2, "vertex"]
  cell.meta[cell.meta$cell_type == root & cell.meta$vertex %in% path2.vertex, "Lin"] <- paste0(root, "_", str_remove(path2, pattern = " "))
  
  # Drop
  cell.meta <- cell.meta[!is.na(cell.meta$Lin), ]
  
  # Subset counts
  normCounts <- normCounts[, rownames(cell.meta)]
  
  # Create scmp
  scmp.obj <- create.scmp(
      counts = normCounts,
      cell_data = cell.meta,
      path_colname = "Lin",
      pseudotime_colname = "Pseudotime"
      
  )
  
  # Return
  return(scmp.obj)
})

scmp.object.list.2 <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
    # Hard Assignment of the Cells to path
    if (rep_i == "rep1") {
        path1 <-"Early Eryth"
        path2 <-"Prog Mk"
        root = "EMP"
        individual <- "Donor-1"
        age <- "35"
        tail = F
        drop =  0.3
        bin.method = "Sturges"
        sex <- "Male"
    } else if (rep_i == "rep2") {
        path1 <-"GMP"
        path2 <-"CLP"
        individual <- "Donor-2"
        root = "HSC"
        age <- "28"
        drop =  0.3
        bin.method = "Sturges"
        tail = T
        sex <- "Female"
    } else if (rep_i == "rep3") {
        path1 <- "CLP"
        path2 <- "GMP"
        root = "LMPP"
        drop =  0.3
        tail = F
        bin.method = "Sturges"
        individual <- "Donor-3"
        age <- "19"
        sex <- "Female"
    }
    
    # Create new object with the names
    scmp.obj <- scmp.object.list[[rep_i]]

  # Sc.Squeeze
  scmp.obj <- sc.squeeze(scmp.obj,
    drop_trails = tail,
    bin_method = bin.method,
    drop_fac = drop,
    split_bins = T,
    bin_pseudotime_colname = "bPseudotime"
  )
  binPlot <- plotBinTile(scmp.obj) + ggtitle(paste(
    individual, "| Age:", age,
    "| sex:", sex
  ))
  binPlot

  # Make Design
  scmp.obj <- sc.set.poly(scmp.obj,
    poly_degree = 2,
  )
  polyGlm <- showPoly(scmp.obj)
  polyGlm

  # # Run p.vector
  scmp.obj <- sc.p.vector(
    parallel = T,
    scmpObj = scmp.obj, verbose = T,
    max_it = 10000,
    logOffset = F,
    family = gaussian(), #MASS::negative.binomial(30),
    useInverseWeights = F,
    logWeights = F,
    useWeights = F,
    offset = F
  )

  if(length(scmp.obj@scPVector@p.vector) == 0){
      return(NULL)
  }else{
  # Run Tstep
  scmp.obj <- sc.T.fit(
    scmpObj = scmp.obj, verbose = T,parallel = T,
    step.method = "backward"
  )

  # # Saving
  saveRDS(
    scmp.obj,
    paste0(outPath, rep_i, "/", "scMaSigPro_Processed_", rep_i, ".RDS")
  )

  print(paste("done", rep_i))

  return(list(
    #scmpObj = scmp.obj,
    polyGlm = polyGlm,
    binPlot = binPlot
  ))
}
})


# Plot bins
compress <- ggarrange(scmp.object.list.2$rep1$binPlot,
                      scmp.object.list.2$rep2$binPlot,
                      scmp.object.list.2$rep3$binPlot,
          nrow = 1)
compress
