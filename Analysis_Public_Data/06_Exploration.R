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
scmp_results <- lapply(rep_vec[1], function(rep_i) {
    rep_i ="rep1"
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    rsq <- 0.6
    num <- 10
    path1_name <-"MPP_to_GMP"
    path2_name <-"MPP_to_CLP"
    root = "EMP"
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    root = "HSC"
    sex <- "Female"
    path1_name <-"EMP_to_ProgEryth"
    path2_name <-"EMP_to_ProgMk"
    rsq <- 0.35
    num <- 10
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    path1_name <-"EMP_to_ProgMk"
    path2_name <-"EMP_to_ProgEryth"
    root = "LMPP"
    rsq <- 0.5
    num <- 10
  }

  # Extract the object
  scmp.obj <- scMaSigPro.list[[rep_i]]
  
  # Add Dummy
  scmp.obj <- sc.filter(scmp.obj,
                        rsq = rsq,
                        significant.intercept = "dummy",
                        vars = "groups")
  
  # get genes
  gene.list <- scmp.obj@sig.genes@sig.genes[[2]]
  cat(length(gene.list))
  
  # Load backgound
  background.vector <- readRDS(paste0("/supp_data/Analysis_Public_Data/", rep_i, "/", rep_i, "background.RDS"))

  # Perform enrichmnet
  target.path = go_enrichment(
      background = background.vector,
      rep = rep_i, 
      age = age, 
      sex = sex,
      path = paste0(str_split_1(names(scmp.obj@sig.genes@sig.genes)[[1]], "vs")[[1]]),
      gene_list = gene.list,
      ont = "BP",
      pAdjustMethod = "BH",
      nterms = 10,
      sig.level = 0.05
  )
  target.path$dot

  
  return(target.path)
})

# Set names
names(scmp_results) <- rep_vec

# Dot
combined.bar <- ggarrange(scmp_results$rep1$dot,
                          scmp_results$rep2$dot,
                          scmp_results$rep3$dot,
  ncol = 3,
  labels = c("A.", "B.", "C.")
)

combined.bar

ggsave(combined.bar,
       filename = paste0("Figures/SuppData/05_Real_Data-GO_dot.png"),
       dpi = 150, height = 8, width = 20
)


plotTrend(scMaSigPro.list$rep1,
          feature_id = "UBXN10",significant = F, logs = F)
