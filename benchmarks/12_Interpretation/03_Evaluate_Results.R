# Title: Evaluate the results of ScMaSigPro on simulated datasets with of lenEq
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(ggpubr))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/12_Interpretation/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Load names of files
dataSets <- list.files(paste0(inPath))

# Load object
scmp.ob <- readRDS(paste0(inPath, dataSets))

# Groups
scmp.ob.none <- sc.get.siggenes(scmp.ob,
                                rsq = 0.7, 
                                vars = "groups", significant.intercept = "none")
scmp.ob.dummy <- sc.get.siggenes(scmp.ob,
                                rsq = 0.7, 
                                vars = "groups", significant.intercept = "dummy")
scmp.ob.all <- sc.get.siggenes(scmp.ob,
                                rsq = 0.7, 
                                vars = "groups", significant.intercept = "all")

# Plot Intersection
none <- sc.path.intersection(scmp.ob.none)
dummy <- sc.path.intersection(scmp.ob.dummy)
all <- sc.path.intersection(scmp.ob.all)
ggarrange(none, dummy, all, ncol = 1)


get.features(
    scmp.ob.dummy,
    query = "unique",
    unique.group = "Path1",
    unique.trend = "up",
)





 get.features(
    scmp.ob.dummy,
    query = "union",
    union.ref.trend = "up",
    union.target.trend = "down",
    union.ref = "Path1", 
    union.target = "Path2vsPath1",
    vars = "each"
)
