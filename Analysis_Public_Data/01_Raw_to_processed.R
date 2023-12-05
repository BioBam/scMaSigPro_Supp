#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Pre-Process and CC Remove ##
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Read BioMart info
biomart.anno <- readRDS(paste(dirPath, "cell_cycle_data.mart", sep = "/"))
reg.out <- c(unique(biomart.anno$SYMBOL))

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Run lapply
violin.list <- mclapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
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

  # Step-2: Load the filtered Matrix
  filt_mat <- Read10X_h5(
    paste(inPath, rep_i, "filtered_feature_bc_matrix.h5", sep = "/")
  )
  # Renaming the features
  filt_mat@Dimnames[[1]] <- str_replace(filt_mat@Dimnames[[1]], "_", "-")

  # Create Raw Seurat Object
  sob.raw <- CreateSeuratObject(
    counts = filt_mat,
    min.cells = 10,
    min.features = 1000
  )

  # Calculate Percentage of Mitochondrial Reads
  sob.raw[["percent.mt"]] <- PercentageFeatureSet(sob.raw, pattern = "^MT-")

  # Remove Cells
  sob.sub <- subset(sob.raw, subset = nFeature_RNA > 100 & nCount_RNA < 30000 & percent.mt < 10)
  
  # Create Subset Plot
  violins.raw.nFeature_RNA <- VlnPlot(sob.raw, features = "nFeature_RNA") + theme(legend.position = "none")  + xlab("") +
      ggtitle(paste(
          individual, "| Age:", age,
          "| sex:", sex
      ), subtitle = "nFeature_RNA")
  violins.raw.nCount_RNA <- VlnPlot(sob.raw, features = "nCount_RNA") + theme(legend.position = "none")  + xlab("") +
      ggtitle(paste(
          individual, "| Age:", age,
          "| sex:", sex
      ), subtitle = "nCount_RNA")
  violins.raw.percent.mt <- VlnPlot(sob.raw, features = "percent.mt") + theme(legend.position = "none")  + xlab("") +
      ggtitle(paste(
          individual, "| Age:", age,
          "| sex:", sex
      ), subtitle = "percent.mt")
  violins.raw <- ggarrange(violins.raw.nFeature_RNA,
                           violins.raw.nCount_RNA,
                           violins.raw.percent.mt,
                           labels = c("A.", "B.", "C."),
                           ncol = 3, nrow = 1)
  # Create Subset Plot
  violins.sub.nFeature_RNA <- VlnPlot(sob.sub, features = "nFeature_RNA") + theme(legend.position = "none") + xlab("") +
      ggtitle(paste(
          individual, "| Age:", age,
          "| sex:", sex
      ), subtitle = "nFeature_RNA")
  violins.sub.nCount_RNA <- VlnPlot(sob.sub, features = "nCount_RNA") + theme(legend.position = "none")  + xlab("") +
      ggtitle(paste(
          individual, "| Age:", age,
          "| sex:", sex
      ), subtitle = "nCount_RNA")
  violins.sub.percent.mt <- VlnPlot(sob.sub, features = "percent.mt") + theme(legend.position = "none")  + xlab("") +
      ggtitle(paste(
          individual, "| Age:", age,
          "| sex:", sex
      ), subtitle = "percent.mt")
  violins.sub <- ggarrange(violins.sub.nFeature_RNA,
                           violins.sub.nCount_RNA,
                           violins.sub.percent.mt,
                           labels = c("D.", "E.", "F."),
                           ncol = 3, nrow = 1)
  # Violin
  violins <- ggarrange(violins.raw,
                       violins.sub,
                       ncol = 1, nrow = 2)
  # Save
  file_name <- paste(dirPath, rep_i, paste(rep_i, "violin.RDS", sep = "_"), sep = "/")
  saveRDS(
      object = violins, file = file_name)
  
  # Normalize
  sob.prs <- NormalizeData(sob.sub, verbose = F)

  # Find Variable features
  sob.prs <- FindVariableFeatures(sob.prs, verbose = F)

  # Check how many genes are present in the dataset
  indata <- rownames(sob.prs)[rownames(sob.prs) %in% reg.out]

  # Keep only those which are present in data
  biomart.anno <- biomart.anno[biomart.anno$SYMBOL %in% indata, ]

  # For seurat we need to divide genes into vectors
  s.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0006260"), "SYMBOL"])
  g2m.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0000087", "GO:0000279", "GO:0007059", "GO:0048285"), "SYMBOL"])

  # Cell Cycle Scores
  sob.prs <- CellCycleScoring(sob.prs,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )

  # Regress Cell Cycle Score
  sob.prs <- ScaleData(sob.prs,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(sob.prs), verbose = T,
  )

  # Save
  file_name <- paste0(dirPath, rep_i, "/", paste0(rep_i, "_prs.RDS"))
  saveRDS(
    object = sob.prs, file = file_name)

  # Return UMAP
  return(NULL)
}, mc.cores = availableCores())
