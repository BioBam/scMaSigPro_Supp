##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(xlsx))

set.seed(007)

# Prefix
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
dirPath <- paste0(base_string, "analysis_public_data/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")

# Load Enrichment Helper
source("R_Scripts/helper_function/go_enrichment.R")

# Create Directory if does not exist
dir.create(tabPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_hd, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_lr, showWarnings = FALSE, recursive = TRUE)

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated", "img"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec

# Set xlsx
excelFile <- paste0(tabPath, "Additional_Table_2_Mechanistic_Analysis_Results.xlsx")
if (file.exists(excelFile)) {
  file.remove(excelFile)
}
write.xlsx(as.data.frame(matrix(data = NA)), excelFile, sheetName = "dummy")

# Load
scmp_cluster_trends <- mclapply(rep_vec, function(rep_i) {
  ob <- readRDS(
    paste0(dirPath, rep_i, "/", "scMaSigPro_Clusters_", rep_i, ".rds")
  )
}, mc.cores = 8)

# Run Go and Extract important gene
scmp_go <- lapply(names(scmp_cluster_trends), function(rep_i) {
  # rep_i = "rep1"
  # get object
  scmp.obj <- scmp_cluster_trends[[rep_i]]

  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor: 1"
    age <- "35"
    sex <- "Male"
    num <- 10
    path_name <- "EMP_ProgMkv"
    gene_set_name <- "intersect"
    sel.clus <- c(2, 3)
  } else if (rep_i == "rep2") {
    individual <- "Donor: 2"
    age <- "28"
    sex <- "Female"
    num <- 10
    path_name <- "CLP_pre-mDC"
    gene_set_name <- "intersect"
    sel.clus <- c(1)
  } else if (rep_i == "rep3") {
    individual <- "Donor: 3"
    age <- "19"
    sex <- "Female"
    num <- 10
    path_name <- "EMP_GMP"
    gene_set_name <- "HSC_GMPvsHSC_EMP"
    sel.clus <- c(1, 3, 4)
  }

  # Create cluster df
  cluster.df <- data.frame(
    cluster = unlist(scmp.obj@Significant@clusters),
    gene = names(unlist(scmp.obj@Significant@clusters))
  )

  # All Clusters
  workbook <- loadWorkbook(excelFile)
  newSheet <- createSheet(workbook, sheetName = paste("All_Clusters", rep_i, sep = "_"))
  addDataFrame(cluster.df, newSheet, row.names = F)
  saveWorkbook(workbook, excelFile)

  # get genes
  gene.list <- cluster.df[cluster.df$cluster %in% sel.clus, ][["gene"]]

  workbook <- loadWorkbook(excelFile)
  newSheet <- createSheet(workbook, sheetName = paste("Gene_list_GO_input", rep_i, sep = "_"))
  addDataFrame(as.data.frame(gene.list), newSheet, row.names = F)
  saveWorkbook(workbook, excelFile)

  cat(length(gene.list))

  # Load backgound
  background.vector <- readRDS(paste0(dirPath, rep_i, "/", rep_i, "background.RDS"))

  # Perform enrichmnet
  target.path <- go_enrichment(
    background = background.vector,
    rep = rep_i,
    age = age,
    sex = sex,
    path = path_name,
    gene_list = gene.list,
    ont = "BP",
    pAdjustMethod = "BH",
    nterms = 10,
    sig.level = 0.05
  )
  target.path$dot
  target.path$ema

  # Add go results
  workbook <- loadWorkbook(excelFile)
  newSheet <- createSheet(workbook, sheetName = paste("GO_BP_Results", rep_i, sep = "_"))
  addDataFrame(as.data.frame(target.path$ego@result), newSheet, row.names = F)
  saveWorkbook(workbook, excelFile)
  return(target.path)
})

# Set names
names(scmp_go) <- rep_vec

# Dot
top <- ggarrange(scmp_go$rep1$dot,
  scmp_go$rep2$dot,
  scmp_go$rep3$dot,
  ncol = 1, nrow = 3,
  labels = c("A.", "C.", "E.")
)
bottom <- ggarrange(scmp_go$rep1$ema,
  scmp_go$rep2$ema + scale_x_continuous(
    name = "X Axis Label",
    limits = c(-5, 0), # Adjust limits as needed
    breaks = seq(-5, 0, 0.5) # Adjust breaks as needed
  ) +
    scale_y_continuous(
      name = "Y Axis Label",
      limits = c(-4, 1), # Adjust limits as needed
      breaks = seq(-4, 1, 0.5) # Adjust breaks as needed
    ),
  scmp_go$rep3$ema,
  ncol = 1, nrow = 3,
  labels = c("B.", "D.", "F.")
)

combined <- ggarrange(top,
  bottom,
  ncol = 2
)
combined

# Save the plot
ggsave(
  plot = combined,
  filename = paste0(figPath_hd, "05_Real_Data_GO_dot.png"),
  dpi = 600, height = 18, width = 20,
  bg = "white"
)
ggsave(
  plot = combined,
  filename = paste0(figPath_lr, "05_Real_Data_GO_dot.png"),
  dpi = 150, height = 18, width = 20,
  bg = "white"
)
