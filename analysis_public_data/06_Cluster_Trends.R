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
imgPath <- paste0(dirPath, "img/")
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
dir.create(imgPath, showWarnings = FALSE, recursive = TRUE)

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated", "img"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec

# Load
scMaSigPro.list <- mclapply(rep_vec, function(rep_i) {
  ob <- readRDS(
    paste0(dirPath, rep_i, "/", "scMaSigPro_Processed_", rep_i, ".RDS")
  )
}, mc.cores = 8)

# Perform hclust
scmp_cluster_trends <- mclapply(rep_vec, function(rep_i) {
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor: 1"
    age <- "35"
    sex <- "Male"
    rsq <- 0.7
    gene_set_name <- "union"
  } else if (rep_i == "rep2") {
    individual <- "Donor: 2"
    age <- "28"
    sex <- "Female"
    rsq <- 0.7
    gene_set_name <- "union"
  } else if (rep_i == "rep3") {
    individual <- "Donor: 3"
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

  # Create trends
  scmp.obj <- sc.cluster.trend(scmp.obj,
    geneSet = gene_set_name,
    cluster_by = "counts",
    k = 6
  )

  # Save Object
  saveRDS(
    scmp.obj,
    paste0(dirPath, rep_i, "/", "scMaSigPro_Clusters_", rep_i, ".rds")
  )
  # Plot trend
  trends <- plotTrendCluster(scmp.obj,
    summary_mode = "median",
    significant = F,
    parallel = T,
    ylab = "Scaled Pseudobulk Expression",
    xlab = "Binned Pseudotime"
  ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.75),
      text = element_text(family = "times")
    ) +
    ggtitle(paste(
      individual, "| Age:", age,
      "| Biological Sex:", sex
    ))

  # return
  return(trends)
}, mc.cores = 1)

# plot
clusters <- ggarrange(scmp_cluster_trends$rep1,
  scmp_cluster_trends$rep2,
  scmp_cluster_trends$rep3,
  nrow = 3
)
clusters

# Save Figures
ggsave(clusters,
  filename = paste0(figPath_hd, "05_Real_Data_clusterTrends.png"),
  dpi = 600, height = 9, width = 7
)
ggsave(clusters,
  filename = paste0(figPath_lr, "05_Real_Data_clusterTrends.png"),
  dpi = 150, height = 9, width = 7
)

scmp_cluster_trends$rep1

# Get the trend
trend <- scmp_cluster_trends$rep1 + theme_minimal(base_size = 18, base_family = "times") +
  theme(
    legend.position = "bottom", panel.grid.major = element_line(
      colour = "lightgrey",
      linewidth = 0.1
    ),
    panel.grid.minor = element_blank(), legend.title = element_blank(),
    strip.background = element_rect(fill = "#434343"), # Background color of the header,
    strip.background.x = element_rect(fill = "#434343"), # Background color of the header
    strip.background.y = element_rect(colour = "#434343"), # Text properties
    strip.text = element_text(colour = "white"), # Text properties
    axis.ticks = element_line(color = "black")
  ) +
  guides(color = guide_legend(
    title = "Branching Paths", # Set the title of the legend
    title.theme = element_text(size = 14, face = "bold", family = "times")
  )) +
  scale_color_manual(
    labels = c("EMP_EarlyErythrocyte" = "Erythroid Lineage", "EMP_ProgMk" = "Megakaryocyte Lineage"),
    values = c("EMP_EarlyErythrocyte" = "#15918A", "EMP_ProgMk" = "#F58A53")
  )
