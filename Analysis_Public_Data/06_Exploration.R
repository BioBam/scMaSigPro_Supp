##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

library(scMaSigPro)
library(stringr)
library(RColorConesa)

set.seed(007)

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec

# Call the required libraries
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(xlsx))

# Set xlsx
excelFile <- "Tables/Additional_Table_2_Mechanistic_Analysis_Results.xlsx"
file.remove(excelFile)

# Load Enrichment Helper
source("R_Scripts/helper_function/go_enrichment.R")

# Load
scMaSigPro.list <- lapply(rep_vec, function(rep_i) {
  ob <- readRDS(
    paste0("/supp_data/Analysis_Public_Data/", rep_i, "/", "scMaSigPro_Processed_", rep_i, ".RDS")
  )
})

# Perform hclust
scmp_cluster_trends <- mclapply(rep_vec, function(rep_i) {
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    rsq <- 0.7
    gene_set_name <- "union"
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    rsq <- 0.7
    gene_set_name <- "union"
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    rsq <- 0.7
    gene_set_name <- "union"
  }

  # Extract the object
  scmp.obj <- scMaSigPro.list[[rep_i]]

  # Add Dummy
  scmp.obj <- sc.filter(scmp.obj,
    rsq = rsq,
    intercept = "dummy",
    vars = "groups"
  )
  plotIntersect(scmp.obj)

  # Create trends
  scmp.obj <- sc.cluster.trend(scmp.obj,
    geneSet = gene_set_name,
    cluster_by = "counts",
    k = 6
  )

  # Plot trend
  trends <- plotTrendCluster(scmp.obj,
    summary_mode = "median",
    significant = F,
    parallel = T
  ) + ggtitle(paste(
    individual, "| Age:", age,
    "| sex:", sex
  ))

  # return
  return(list(
    trends = trends,
    scmp.obj = scmp.obj
  ))
}, mc.cores = 1)

# plot
clusters <- ggarrange(scmp_cluster_trends$rep1$trends,
  scmp_cluster_trends$rep2$trends,
  scmp_cluster_trends$rep3$trends,
  nrow = 3
)
clusters

ggsave(clusters,
  filename = paste0("Figures/SuppData/05_Real_Data_clusterTrends.png"),
  dpi = 600, height = 8, width = 6
)

# Create Dummy xlsx
write.xlsx(as.data.frame(matrix(data = NA)), excelFile, sheetName = "dummy")

# Run Go and Extract important gene
scmp_results <- lapply(names(scmp_cluster_trends), function(rep_i) {
  # get object
  scmp.obj <- scmp_cluster_trends[[rep_i]][["scmp.obj"]]

  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
    num <- 10
    path_name <- "EMP_ProgMkv"
    gene_set_name <- "intersect"
    sel.clus <- c(2, 3)
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
    num <- 10
    path_name <- "CLP_pre-mDC"
    gene_set_name <- "intersect"
    sel.clus <- c(1)
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
    num <- 10
    path_name <- "EMP_GMP"
    gene_set_name <- "HSC_GMPvsHSC_EMP"
    sel.clus <- c(1, 3)
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
  background.vector <- readRDS(paste0("/supp_data/Analysis_Public_Data/", rep_i, "/", rep_i, "background.RDS"))

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
names(scmp_results) <- rep_vec

# Dot
top <- ggarrange(scmp_results$rep1$dot,
  scmp_results$rep2$dot,
  scmp_results$rep3$dot,
  ncol = 1, nrow = 3,
  labels = c("A.", "C.", "E.")
)
bottom <- ggarrange(scmp_results$rep1$ema,
  scmp_results$rep2$ema,
  scmp_results$rep3$ema,
  ncol = 1, nrow = 3,
  labels = c("B.", "D.", "F.")
)

combined <- ggarrange(top,
  bottom,
  ncol = 2
)
combined

ggsave(combined,
  filename = paste0("Figures/SuppData/05_Real_Data_GO_dot.png"),
  dpi = 300, height = 16, width = 16
)


# Check and Plot markers
az_emp <- c("MYCT1", "CRHBP", "NPR3", "AVP", "GATA2", "HPGDS", "CYTL1", "CRYGD", "IGSF10", "PBX1")
az_gmp <- c("SERPINB10", "RNASE3", "MS4A3", "PRTN3", "ELANE", "AZU1", "CTSG", "RNASE2", "RETN", "NPW")
az_earlyE <- c("CNRIP1", "GATA2", "ITGA2B", "TFR2", "GATA1", "KLF1", "CYTL1", "MAP7", "FSCN1", "APOC1")
az_hsc <- c("CRHBP", "AVP", "MYCT1", "BEX1", "NPR3", "CRYGD", "MSRB3", "CD34", "NPDC1", "MLLT3")
az_proB <- c("CYGB", "UMODL1", "EBF1", "MME", "VPREB1", "DNTT", "IGLL1", "UHRF1", "BLNK", "AGPS")
az_progMk <- c("CLEC1B", "SPX", "WFDC1", "ANXA3", "CMTM5", "SELP", "RBPMS2", "ARHGAP6", "GP9", "LTBP1")
az_clp <- c("ACY3", "PRSS2", "C1QTNF4", "SPINK2", "SMIM24", "NREP", "CD34", "DNTT", "FLT3", "SPNS3")
az_pPDC <- c("SCT", "SHD", "LILRA4", "LILRB4", "PTPRS", "TNNI2", "PLD4", "SPIB", "IRF8", "TNFRSF21")
az_mPDC <- c("ENHO", "CLEC10A", "RNASE2", "PLBD1", "FCER1A", "IGSF6", "MNDA", "SAMHD1", "ALDH2", "PAK1")

# Plot Markers

GP9 <- plotTrend(scmp_cluster_trends$rep1$scmp.obj,
  feature_id = "GP9", significant = F, logs = F,
  pseudoCount = F
)

plotTrend(scmp_cluster_trends$rep2$scmp.obj,
  feature_id = "MNDA", significant = F, logs = F,
  pseudoCount = F
)

plotTrend(scmp_cluster_trends$rep3$scmp.obj,
  feature_id = "GATA1", significant = F, logs = F,
  pseudoCount = F
)

# Prepare figure for artcile
GP9 <- GP9 + ggtitle(
  "Glycoprotein IX (GP9), R-Square: 1",
  "Donor:1 | Age: 35 | Biological Sex: Male"
) +
  xlab("Binned Pseudotime") + ylab("Scaled Pesudobulk Expressied") +
  theme_classic(base_size = 14) +
  theme(
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
    # legend.position = "bottom"
    legend.position = c(0.1, 0.8), legend.justification = c("left", "top")
  ) +
  guides(color = guide_legend(key_width = unit(5, "cm"), key_height = unit(4, "cm"))) +
  scale_color_manual(
    labels = c(
      "Erythroid Lineage",
      "Megakaryocyte Lineage"
    ),
    values = c(
      colorConesa(6)[2],
      colorConesa(6)[1]
    ),
  )

# save
saveRDS(GP9,
  file = "Figures/MainArticle/MainArticle_FigureD.RDS"
)
