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
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec

# Call the required libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(ggpubr))

# Load Enrichment Helper
source("R_Scripts/helper_function/go_enrichment.R")

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
    num <- 10
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    root = "HSC"
    sex <- "Female"
    rsq <- 0.8
    num <- 10
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    root = "LMPP"
    rsq <- 0.8
    num <- 10
  }

  # Extract the object
  scmp.obj <- scMaSigPro.list[[rep_i]]
  
  # Add Dummy
  scmp.obj <- sc.filter(scmp.obj,
                        rsq = rsq,
                        significant.intercept = "dummy",
                        vars = "groups")
  
  if(rep_i == "rep3"){
      gene.list <- scmp.obj@sig.genes@sig.genes[[2]]
      path_name <- names(scmp.obj@sig.genes@sig.genes)[2]
  }else{
      gene.list <- scmp.obj@sig.genes@sig.genes[[1]]
      path_name <-paste0(str_split_1(names(scmp.obj@sig.genes@sig.genes)[[1]], "vs")[[1]])
  }
  # get genes
  cat(length(gene.list))
  
  # Load backgound
  background.vector <- readRDS(paste0("/supp_data/Analysis_Public_Data/", rep_i, "/", rep_i, "background.RDS"))

  # Perform enrichmnet
  target.path = go_enrichment(
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

  
  return(list(target.path = target.path,
              scmp.obj = scmp.obj))
})

# Set names
names(scmp_results) <- rep_vec

# Dot
combined.bar <- ggarrange(scmp_results$rep1$target.path$dot,
                          scmp_results$rep2$target.path$dot,
                          scmp_results$rep3$target.path$dot,
  ncol = 3,
  labels = c("A.", "B.", "C.")
)

combined.bar

ggsave(combined.bar,
       filename = paste0("Figures/SuppData/05_Real_Data-GO_dot.png"),
       dpi = 150, height = 8, width = 20
)


# Check and Plot markers
az_emp <- c("MYCT1", "CRHBP", "NPR3", "AVP", "GATA2", "HPGDS", "CYTL1", "CRYGD", "IGSF10", "PBX1")
az_gmp <- c("SERPINB10", "RNASE3", "MS4A3", "PRTN3", "ELANE", "AZU1", "CTSG", "RNASE2", "RETN", "NPW")
az_earlyE <- c("CNRIP1", "GATA2", "ITGA2B", "TFR2", "GATA1", "KLF1", "CYTL1", "MAP7", "FSCN1", "APOC1")
az_hsc <- c("CRHBP", "AVP", "MYCT1", "BEX1", "NPR3", "CRYGD", "MSRB3", "CD34", "NPDC1", "MLLT3")
az_proB <- c("CYGB", "UMODL1", "EBF1", "MME", "VPREB1", "DNTT", "IGLL1", "UHRF1", "BLNK", "AGPS")
az_progMk <- c("CLEC1B", "SPX", "WFDC1", "ANXA3", "CMTM5", "SELP", "RBPMS2", "ARHGAP6", "GP9", "LTBP1")
az_clp <- c("ACY3", "PRSS2", "C1QTNF4", "SPINK2", "SMIM24", "NREP", "CD34", "DNTT", "FLT3", "SPNS3")
az_pPDC <- c("SCT","SHD","LILRA4","LILRB4","PTPRS","TNNI2","PLD4","SPIB","IRF8","TNFRSF21")

# For rep1
scmp.ob.rep1 <- scmp_results$rep1$scmp.obj
scmp.ob.rep2 <- scmp_results$rep2$scmp.obj
scmp.ob.rep3 <- scmp_results$rep3$scmp.obj

# For ep1
plotTrend(scmp.ob.rep1,
          feature_id = "CMTM5",significant = F, logs = F,
          pseudoCount = F)

plotTrend(scmp.ob.rep2,
          feature_id = "SCT",significant = F, logs = F,
          pseudoCount = F)

plotTrend(scmp.ob.rep3,
          feature_id = "GATA1",significant = F, logs = F,
          pseudoCount = F)



# Perform clustering
plotTrendCluster(scmp.ob.rep1,geneSet = "EMP_ProgMkvsEMP_EarlyErythrocyte",
                 cluster_by = "counts",
                 logs = F, smoothness = 1,
                 includeInflu = T)

plotTrendCluster(scmp.ob.rep2,geneSet = "CLP_pre_pDCvsCLP_pre_mDC",
                 cluster_by = "counts",
                 logs = F, smoothness = 1,
                 includeInflu = T)

plotTrendCluster(scmp.ob.rep3,geneSet = "EMP_EarlyErythrocyte",
                 cluster_by = "counts",
                 logs = F, smoothness = 1,
                 includeInflu = T)

