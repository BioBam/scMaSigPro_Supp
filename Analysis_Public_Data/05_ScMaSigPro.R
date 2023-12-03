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
    
  # get object
  rep_i_obj <- object.list[[rep_i]]

  # Create Cell level metadata
  # Step-1: Add barcodes
  cell.data <- data.frame(
    barcodes = colnames(assay(rep_i_obj)),
    row.names = colnames(assay(rep_i_obj))
  )

  # Add Cell type
  cell.data$cell_type <- rep_i_obj@colData$cell_type

  # Add Pseudotime to cds and cell data
  cell.data$pseudotime <- pseudotime(rep_i_obj)

  # Hard Assignment of the Cells to path
  if (rep_i == "rep1") {
    path1 <-"Early Eryth"
    path2 <-"Prog Mk"
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
  } else if (rep_i == "rep2") {
    path1 <-"pre B"
    path2 <-"Early Eryth"
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
  } else if (rep_i == "rep3") {
    path1 <- "EMP"
    path2 <- "CLP"
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
  }

  # Assign paths
  cell.data[cell.data$cell_type %in% path1, "path"] <- paste(str_remove(path1, pattern = " "), collapse = "-")
  cell.data[cell.data$cell_type %in% path2, "path"] <- paste(str_remove(path2, pattern = " "), collapse = "-")
  cell.data <- cell.data[!is.na(cell.data$path), ]
  
  # Drop any infinite time
  cell.data <- cell.data[!is.infinite(cell.data$pseudotime),]

  # Extract Raw Counts and subset
  raw_counts <- assay(rep_i_obj)
  raw_counts <- raw_counts[, rownames(cell.data)]

  # Drop gene
  raw_counts <- raw_counts[rowSums(raw_counts) >= 100, ]
  
  # Create SCMP Object
  scmp.obj <- create.scmp(
    counts = raw_counts,
    cell_data = cell.data,
    pseudotime_colname = "pseudotime",
    path_colname = "path"
  )
  # scmp.obj <- selectPath.m3(rep_i_obj,
  #                           annotation_col = "cell_type")

  # Sc.Squeeze
  scmp.obj <- sc.squeeze(scmp.obj,
    drop_trails = T,
    bin_method = "Doane",
    drop_fac = 1,
    bin_pseudotime_colname = "bPseudotime"
  )
  binPlot <- plotBinTile(scmp.obj) + ggtitle(paste(
    individual, "| Age:", age,
    "| sex:", sex
  ))
  binPlot

  # Make Design
  scmp.obj <- sc.set.poly(scmp.obj,
    poly_degree = 3,
  )
  polyGlm <- showPoly(scmp.obj)

  # Run p.vector
  scmp.obj <- sc.p.vector(
    parallel = T,
    scmpObj = scmp.obj, verbose = T,
    max_it = 10000,
    logOffset = F,
    family = gaussian, #MASS::negative.binomial(30),
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
    scmpObj = scmp.obj, verbose = T,
    step.method = "backward"
  )

  # Saving
  saveRDS(
    scmp.obj,
    paste0(outPath, rep_i, "/", "scMaSigPro_Processed_", rep_i, ".RDS")
  )
  
  print(paste("done", rep_i))

  return(list(
    scmpObj = scmp.obj,
    polyGlm = polyGlm,
    binPlot = binPlot
  ))}
})
