# Title: Run ScMaSigPro on simulated datasets with different Skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/12_Interpretation/simulated/sce/"
outPath <- "/supp_data/benchmarks/12_Interpretation/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = "\\.", i = 2),
  ".RData"
)

# Load Data
load(file = paste0(inPath, dataSets["50"]))

# Convert
scmp.obj <- as_scmp(sim.sce,
  from = "sce",
  additional_params = list(
    labels_exist = TRUE,
    existing_pseudotime_colname = "Step",
    existing_path_colname = "Group"
  ), verbose = F
)

# Compress
scmp.obj <- squeeze(
  scmpObject = scmp.obj,
  bin_method = "Sturges",
  drop_fac = 1,
  verbose = F,
  cluster_count_by = "sum",
  split_bins = F,
  prune_bins = F,
  drop_trails = F,
  fill_gaps = F
)

# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj, poly_degree = 2)

# Run p-vector
scmp.obj <- sc.p.vector(
  scmpObj = scmp.obj, verbose = F, min.obs = 1, parallel = T,
  offset = T, logOffset = T,
  useWeights = T, logWeights = F, useInverseWeights = F,
  max_it = 1000,
  globalTheta = FALSE
)

# Run-Step-2
scmp.obj <- sc.T.fit(
  parallel = T,
  scmpObj = scmp.obj, verbose = F,
  step.method = "backward"
)

# Extract gene sets
scmp.obj <- sc.get.siggenes(scmp.obj,
                significant.intercept = "dummy",
                includeInflu = T,
                vars = "groups")

sc.path.intersection(scmp.obj, keep_empty_groups = T)

#Export per cetegory
pathB_vs_pathA <- rownames(scmp.obj@sig.genes@sig.genes$Path2vsPath1$coefficients)
pathA <- rownames(scmp.obj@sig.genes@sig.genes$Path1$coefficients)
pathB_vs_pathA_and_pathA <- intersect(pathB_vs_pathA, pathA)
length(pathB_vs_pathA_and_pathA)

pathB_vs_pathA <- rownames(scmp.obj@sig.genes@sig.genes$Path2vsPath1$coefficients)
pathA <- rownames(scmp.obj@sig.genes@sig.genes$Path1$coefficients)
pathB_vs_pathA <- setdiff(pathB_vs_pathA, pathA)
length(pathB_vs_pathA)

pathB_vs_pathA <- rownames(scmp.obj@sig.genes@sig.genes$Path2vsPath1$coefficients)
pathA <- rownames(scmp.obj@sig.genes@sig.genes$Path1$coefficients)
pathA <- setdiff(pathA, pathB_vs_pathA)
length(pathA)

pathA_list <- list()
# Check patterns
for (i in pathA){
    pathA_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj, 
                  feature_id = i,
                  logs = T, logType = "log")
}

pathB_vs_pathA_list <- list()
# Check patterns
for (i in pathB_vs_pathA){
    pathB_vs_pathA_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj, 
                                    feature_id = i,
                                    logs = T, logType = "log")
}


pathB_vs_pathA_and_pathA_list <- list()
# Check patterns
for (i in pathB_vs_pathA_and_pathA){
    pathB_vs_pathA_and_pathA_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj, 
                                             feature_id = i,
                                             logs = T, logType = "log")
}

# Pvalue beta estimates
beta0 <- showSol(scmp.obj, view = F, return = T)
beta0 <- beta0[beta0$`p-value` <= 0.5, ,drop  = FALSE]
beta0 <- beta0[beta0$`R-squared` >= 0.7, ,drop  = FALSE]
beta0 <- beta0[beta0$p.valor_beta0 <= 0.05, ,FALSE]

# Pvalue beta estimates
Path2vsPath1 <- showSol(scmp.obj, view = F, return = T)
Path2vsPath1 <- Path2vsPath1[beta0$`p-value` <= 0.5, ,drop  = FALSE]
Path2vsPath1 <- Path2vsPath1[beta0$`R-squared` >= 0.7, ,drop  = FALSE]
Path2vsPath1[is.na(Path2vsPath1$p.valor_Path2vsPath1), "p.valor_Path2vsPath1"] <- as.integer(1)
Path2vsPath1[is.na(Path2vsPath1$p.valor_scmp_binned_pseudotime), "p.valor_scmp_binned_pseudotime"] <-as.integer(1)
Path2vsPath1[is.na(Path2vsPath1$p.valor_scmp_binned_pseudotimexPath2), "p.valor_scmp_binned_pseudotimexPath2"] <-as.integer(1)
Path2vsPath1[is.na(Path2vsPath1$p.valor_scmp_binned_pseudotime2), "p.valor_scmp_binned_pseudotime2"] <-as.integer(1)
Path2vsPath1[is.na(Path2vsPath1$p.valor_scmp_binned_pseudotime2xPath2), "p.valor_scmp_binned_pseudotime2xPath2"] <-as.integer(1)
Path2vsPath1 <- Path2vsPath1[(Path2vsPath1$p.valor_Path2vsPath1 <= 0.05 &
                                  Path2vsPath1$p.valor_Path2vsPath1 > 0.05), ,FALSE]
Path2vsPath1_list <- list()
for (i in rownames(Path2vsPath1)){
    Path2vsPath1_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj, 
                                                       feature_id = i,
                                                       logs = T, logType = "log")
}
