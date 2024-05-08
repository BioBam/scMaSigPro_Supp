##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(scMaSigPro))

# Prefix
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
dirPath <- paste0(base_string, "analysis_public_data/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec

# Load object
object.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_cds.RDS"))
})

# Create ScMaSigPro Input
for (rep_i in c(rep_vec)) {
  # rep_i = "rep3"

  if (rep_i == "rep1") {
    path1_name <- "EMP_EarlyErythrocyte"
    path2_name <- "EMP_ProgMk"
    root_pp <- c("Y_39")
    path1_pp <- c("Y_39", "Y_62", "Y_85", "Y_144", "Y_160", "Y_167", "Y_168")
    path2_pp <- c("Y_2", "Y_3", "Y_39", "Y_43", "Y_97")
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
  } else if (rep_i == "rep2") {
    path1_name <- "CLP_pre_mDC"
    path2_name <- "CLP_pre_pDC"
    root_pp <- c("Y_17")
    path1_pp <- c("Y_2", "Y_17", "Y_21", "Y_24", "Y_25", "Y_26", "Y_27", "Y_30", "Y_37", "Y_38", "Y_40")
    path2_pp <- c("Y_4", "Y_8", "Y_9", "Y_10", "Y_11", "Y_13", "Y_17", "Y_20", "Y_28", "Y_29", "Y_31", "Y_41")
    age <- "28"
    sex <- "Female"
  } else if (rep_i == "rep3") {
    path1_name <- "HSC_EMP"
    path2_name <- "HSC_GMP"
    root_pp <- c("Y_206")
    path2_pp <- c(
      "Y_1", "Y_23", "Y_25", "Y_27", "Y_50", "Y_57", "Y_61", "Y_63", "Y_70", "Y_71",
      "Y_86", "Y_89", "Y_91", "Y_105", "Y_143", "Y_178", "Y_180", "Y_183", "Y_189",
      "Y_193", "Y_217", "Y_235", "Y_249", "Y_266", "Y_275", "Y_284"
    )
    path1_pp <- c("Y_42", "Y_162", "Y_213", "Y_249", "Y_254")
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
  }

  ob <- m3_select_path(
    cds = object.list[[rep_i]],
    use_shiny = F, plot_purity = F,
    anno_col = "cell_type",
    m3_pp = list(
      root_pp = root_pp,
      path1_pp = path1_pp,
      path2_pp = path2_pp,
      path1_name = path1_name,
      path2_name = path2_name
    )
  )

  saveRDS(
    ob,
    paste0(dirPath, rep_i, "/", "scMaSigPro_Input_.RDS")
  )
}
