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
rep_vec <- rep_vec[!(rep_vec %in% c("Azimuth_Human_BoneMarrow", "integrated"))]
names(rep_vec) <- rep_vec
rep_vec <- rep_vec[3]

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(scMaSigPro))

# Load object
scmp.ob.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  sob <- readRDS(file = paste0(inPath, rep_i, "/scMaSigPro_Input_.RDS"))
})

# Run ScMaSigPro
scmp.prs.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # rep_i <- "rep3"

  # Hard Assignment of the Cells to path
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    tail <- T
    drop <- 1
    split <- F
    bin.method <- "Sturges"
    sex <- "Male"
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    tail <- F
    drop <- 1
    split <- F
    bin.method <- "Sturges"
    sex <- "Female"
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    tail <- F
    drop <- 1
    split <- F
    bin.method <- "Sturges"
    sex <- "Female"
  }

  # Create new object with the names
  scmp.obj <- scmp.ob.list[[rep_i]]

  # Sc.Squeeze
  scmp.obj <- sc.squeeze(scmp.obj,
    drop_trails = tail,
    bin_method = bin.method,
    drop_fac = drop,
    split_bins = split,
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
  polyGlm

  # # Run p.vector
  scmp.obj <- sc.p.vector(
    parallel = T,
    scmpObj = scmp.obj,
    verbose = T,
    max_it = 10000,
    logOffset = F,
    family = gaussian(), # MASS::negative.binomial(30),
    offset = F
  )

  if (length(scmp.obj@profile@non.flat) < 100) {
    return(list(
      scmpObj = scmp.obj,
      polyGlm = polyGlm,
      binPlot = binPlot
    ))
  } else {
    # Run Tstep
    scmp.obj <- sc.t.fit(
      scmpObj = scmp.obj, verbose = T, parallel = T,
      step.method = "backward"
    )

    # # Saving
    saveRDS(
      scmp.obj,
      paste0(outPath, rep_i, "/", "scMaSigPro_Processed_", rep_i, ".RDS")
    )
    return(list(
      scmpObj = scmp.obj,
      polyGlm = polyGlm,
      binPlot = binPlot
    ))
  }
})

scmp.prs.list$rep1$scmpObj
scmp.prs.list$rep2$scmpObj
scmp.prs.list$rep3$scmpObj


# Plot bins
compress <- ggarrange(scmp.prs.list$rep1$binPlot,
  scmp.prs.list$rep2$binPlot,
  scmp.prs.list$rep3$binPlot,
  nrow = 1
)
compress
