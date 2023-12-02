##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

library(scMaSigPro)

set.seed(007)

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Call the required libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(ggpubr))

# Load Enrichment Helper
source("R_Scripts/helper_function/go_enrichment.R")

# Define Marker gene list from Azimuth
az_emp <- c("MYCT1", "CRHBP", "NPR3", "AVP", "GATA2", "HPGDS", "CYTL1", "CRYGD", "IGSF10", "PBX1")
az_gmp <- c("SERPINB10", "RNASE3", "MS4A3", "PRTN3", "ELANE", "AZU1", "CTSG", "RNASE2", "RETN", "NPW")
az_earlyE <- c("CNRIP1", "GATA2", "ITGA2B", "TFR2", "GATA1", "KLF1", "CYTL1", "MAP7", "FSCN1", "APOC1")
az_hsc <- c("CRHBP", "AVP", "MYCT1", "BEX1", "NPR3", "CRYGD", "MSRB3", "CD34", "NPDC1", "MLLT3")
az_proB <- c("CYGB", "UMODL1", "EBF1", "MME", "VPREB1", "DNTT", "IGLL1", "UHRF1", "BLNK", "AGPS")
az_progMk <- c("CLEC1B", "SPX", "WFDC1", "ANXA3", "CMTM5", "SELP", "RBPMS2", "ARHGAP6", "GP9", "LTBP1")
az_clp <- c("ACY3", "PRSS2", "C1QTNF4", "SPINK2", "SMIM24", "NREP", "CD34", "DNTT", "FLT3", "SPNS3")

# Load
scMaSigPro.list <- lapply(rep_vec, function(rep_i) {
  ob <- readRDS(
    paste0("/supp_data/Analysis_Public_Data/", rep_i, "/", "scMaSigPro_Processed_", rep_i, ".RDS")
  )
})

# Run Go and Extract important gene
scmp_results <- lapply(rep_vec, function(rep_i) {
    #rep_i <- "rep3"
    
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    rsq <- 0.6
    num <- 5
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    rsq <- 0.6
    num <- 5
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    rsq <- 0.6
    num <- 5
  }

  # Extract the object
  scmp.obj <- scMaSigPro.list[[rep_i]]

  # Get paths
  path.vec <- unique(scmp.obj@scTFit@groups.vector)
  names(path.vec) <- c("ref", "target") 
  
  # Extract Significant Genes in path1
  sig.features.path1 <- sc.get.features(
    scmpObj = scmp.obj,
    query = "union",
    rsq = rsq,
    significant.intercept = "dummy",
    vars = "groups",
    includeInflu = T,
    union.ref = path.vec[["ref"]],
    union.target = path.vec[["target"]],
    unique.group = path.vec[["ref"]],
    union.ref.trend = "down",
    union.target.trend = "up"
  )

  # Extract Significant Genes in path2
  sig.features.path2 <- sc.get.features(
      scmpObj = scmp.obj,
      query = "unique",
      rsq = rsq,
      significant.intercept = "dummy",
      vars = "groups",
      includeInflu = T,
      union.ref = path.vec[["ref"]],
      union.target = path.vec[["target"]],
      unique.group = path.vec[["target"]],
      union.ref.trend = "down",
      union.target.trend = "up"
  )
  
  # Order
  sig.features.path1 <- sig.features.path1[order(sig.features.path1)]
  sig.features.path2 <- sig.features.path2[order(sig.features.path2)]

  path1.name <- str_split_1(path.vec[["target"]], "vs")[2]
  path2.name <-str_split_1(path.vec[["target"]], "vs")[1]
  # Run enrichmnet
  enrichmnet.list <- list(
    path1 = go_enrichment(
      scmp.obj = scmp.obj,
      rep = rep_i, age = age, sex = sex,
      path = path1.name,
      gene_list = sig.features.path1,
      ont = "BP", pAdjustMethod = "BH",
      nterms = num
    ),
    path2 = go_enrichment(
      scmp.obj = scmp.obj,
      rep = rep_i, age = age, sex = sex,
      path =path2.name,
      gene_list = sig.features.path2,
      ont = "BP", pAdjustMethod = "BH",
      nterms = num
    )
  )

  # Set Name
  names(enrichmnet.list) <- c(path1.name, path1.name)

  # Return
  return(list(
    obj = scmp.obj,
    go_results = enrichmnet.list
  ))
})

# Set names
names(scmp_results) <- rep_vec

# Dot list
dot.list <- list()
# get plot list
for (rep in rep_vec) {
  replist <- scmp_results[[rep]][["go_results"]]
  path1 <- replist[[1]][["dot"]]
  path2 <- replist[[2]][["dot"]]
  dot.list[[rep]] <- ggarrange(path1, path2, nrow = 2, ncol = 1)
}

ggarrange(dot.list$rep1,
  dot.list$rep2,
  dot.list$rep3,
  ncol = 3,
  labels = c("A.", "B.", "C.")
)
