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

# Load
scmp_cluster_trends <- mclapply(rep_vec, function(rep_i) {
  ob <- readRDS(
    paste0(dirPath, rep_i, "/", "scMaSigPro_Clusters_", rep_i, ".rds")
  )
}, mc.cores = 8)

# Set objects
scmp_d1 <- scmp_cluster_trends$rep1
scmp_d2 <- scmp_cluster_trends$rep2
scmp_d3 <- scmp_cluster_trends$rep3

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

# Create Markers of the genes each of the donor
donor_1_list <- list(EMP_ProgMk = az_progMk, EMP_EarlyErythrocyte = az_earlyE)
donor_2_list <- list(CLP_pre_pDC = az_pPDC, CLP_pre_mDC = az_mPDC)
donor_3_list <- list(HSC_GMP = az_gmp, HSC_EMP = az_emp)

donor_1_markers <- list()
donor_2_markers <- list()
donor_3_markers <- list()

# Get markers
for (name in names(donor_1_list)) {
  marker_vector <- donor_1_list[[name]]
  detect <- unique(unlist(scmp_d1@Significant@genes))
  detect_markers <- intersect(marker_vector, detect)
  donor_1_markers[[name]] <- detect_markers
}
for (name in names(donor_2_list)) {
  marker_vector <- donor_2_list[[name]]
  detect <- unique(unlist(scmp_d2@Significant@genes))
  detect_markers <- intersect(marker_vector, detect)
  donor_2_markers[[name]] <- detect_markers
}
for (name in names(donor_3_list)) {
  marker_vector <- donor_3_list[[name]]
  detect <- unique(unlist(scmp_d3@Significant@genes))
  detect_markers <- intersect(marker_vector, detect)
  donor_3_markers[[name]] <- detect_markers
}

# Create vectors and  sample 4 genes randomly
set.seed(321)
donor_1_markers <- sample(unlist(donor_1_markers), 3)
set.seed(123)
donor_2_markers <- sample(unlist(donor_2_markers), 3)
donor_3_markers <- sample(unlist(donor_3_markers), 3)

# Plot the genes
donor1.plots <- ggarrange(
  plotTrend(scmp_d1,
    feature_id = donor_1_markers[1], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  plotTrend(scmp_d1,
    feature_id = donor_1_markers[2], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  plotTrend(scmp_d1,
    feature_id = donor_1_markers[3], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  labels = c("A.", "B.", "C."), common.legend = TRUE,
  legend = "bottom",
  nrow = 1
)


donor2.plots <- ggarrange(
  plotTrend(scmp_d2,
    feature_id = donor_2_markers[1], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  plotTrend(scmp_d2,
    feature_id = donor_2_markers[2], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  plotTrend(scmp_d2,
    feature_id = donor_2_markers[3], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  labels = c("D.", "E.", "F."), common.legend = TRUE, legend = "bottom",
  nrow = 1
)


donor3.plots <- ggarrange(
  plotTrend(scmp_d2,
    feature_id = donor_3_markers[1], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  plotTrend(scmp_d2,
    feature_id = donor_3_markers[2], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  plotTrend(scmp_d2,
    feature_id = donor_3_markers[3], significant = F, logs = F,
    pseudoCount = F, xlab = "Binned Pseudotime"
  ),
  labels = c("G.", "H.", "I."), common.legend = TRUE, legend = "bottom",
  nrow = 1
)


markers <- ggarrange(donor1.plots,
  donor2.plots,
  donor3.plots,
  nrow = 3
)

markers

# Save the plot
ggsave(
  plot = markers,
  filename = paste0(figPath_hd, "05_Real_Data_Markers.png"),
  dpi = 600, height = 8, width = 10,
  bg = "white"
)
ggsave(
  plot = markers,
  filename = paste0(figPath_lr, "05_Real_Data_Markers.png"),
  dpi = 150, height = 8, width = 10,
  bg = "white"
)


GP9 <- plotTrend(scmp_d1,
  feature_id = "GP9", significant = F, logs = F,
  pseudoCount = F
)

# Prepare figure for Main Article
GP9 <- GP9 + ggtitle(
  "Glycoprotein IX (GP9), R-Square: 1",
  "Donor:1 | Age: 35 | Biological Sex: Male"
) +
  xlab("Binned Pseudotime") + ylab("Scaled Pesudobulk Expression") +
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
    labels = c("EMP_EarlyErythrocyte" = "Erythroid Lineage", "EMP_ProgMk" = "Megakaryocyte Lineage"),
    values = c("EMP_EarlyErythrocyte" = "#15918A", "EMP_ProgMk" = "#F58A53")
  )

GP9 <- GP9 + theme_minimal(base_size = 18, base_family = "times") +
  theme(
    legend.position = "bottom", panel.grid.major = element_line(
      colour = "lightgrey",
      linewidth = 0.1
    ),
    panel.grid.minor = element_blank(), legend.title = element_blank(),
    strip.background = element_rect(fill = "#434343"), # Background color of the header
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

# Load trend before this
combined <- ggarrange(
  trend,
  GP9,
  widths = c(1.7, 1),
  nrow = 1, ncol = 2,
  labels = c("A.", "B."),
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  plot = combined, filename = "Main_Article_Figure_3.png",
  path = figPath_hd,
  width = 16, dpi = 600, height = 6,
  bg = "white"
)

ggsave(
  plot = combined, filename = "Main_Article_Figure_3.png",
  path = figPath_lr,
  width = 16, dpi = 300, height = 6,
  bg = "white"
)
