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
                significant.intercept = "none",
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

sol <- showSol(scmp.obj, view = F, return = T)
sol[is.na(sol)] <- 1
sol <- sol[sol$`p-value` <= 0.5, ,drop  = FALSE]
sol <- sol[sol$`R-squared` >= 0.7, ,drop  = FALSE]
nrow(sol)


# Group_patterns
g1 <- sol[(sol$p.valor_Path2vsPath1 > 0.05 &
               sol$p.valor_scmp_binned_pseudotimexPath2 > 0.05 &
               sol$p.valor_scmp_binned_pseudotime2xPath2 > 0.05 &
               sol$p.valor_scmp_binned_pseudotime2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime <= 0.05), , drop =FALSE]
g1_list <- list()
for (i in rownames(g1)){
    g1_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj,
                                 feature_id = i,
                                 logs = T, logType = "log")
}
g1_list[[sample(names(g1_list), 1)]]
length(g1_list)

# Group_patterns
g2 <- sol[(sol$p.valor_Path2vsPath1 > 0.05 &
               sol$p.valor_scmp_binned_pseudotimexPath2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime2xPath2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime2 > 0.05 &
               sol$p.valor_scmp_binned_pseudotime > 0.05), , drop =FALSE]
g2_list <- list()
for (i in rownames(g2)){
    g2_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj,
                                 feature_id = i,
                                 logs = T, logType = "log")
}
g2_list[[sample(names(g2_list), 1)]]
length(g2_list)


# Group_patterns
g3 <- sol[(sol$p.valor_Path2vsPath1 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotimexPath2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime2xPath2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime2 > 0.05 &
               sol$p.valor_scmp_binned_pseudotime > 0.05), , drop =FALSE]
g3_list <- list()
for (i in rownames(g3)){
    g3_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj,
                                 feature_id = i,
                                 logs = T, logType = "log")
}
g3_list[[sample(names(g3_list), 1)]]
length(g3_list)

# Group_patterns
g4 <- sol[(sol$p.valor_Path2vsPath1 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotimexPath2 > 0.05 &
               sol$p.valor_scmp_binned_pseudotime2xPath2 > 0.05 &
               sol$p.valor_scmp_binned_pseudotime2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime <= 0.05), , drop =FALSE]
g4_list <- list()

for (i in rownames(g4)){
    g4_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj,
                                 feature_id = i,
                                 logs = T, logType = "log")
}
g4_list[[sample(names(g4_list), 1)]]
length(g4_list)

# Group_patterns
g5 <- sol[(sol$p.valor_Path2vsPath1 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotimexPath2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime2xPath2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime2 <= 0.05 &
               sol$p.valor_scmp_binned_pseudotime <= 0.05), , drop =FALSE]
g5_list <- list()

for (i in rownames(g5)){
    g5_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj,
                                 feature_id = i,
                                 logs = T, logType = "log")
}
g5_list[[sample(names(g5_list), 1)]]
length(g5_list)

# remaining
sol_rem <- sol[!(rownames(sol) %in% c(rownames(g1), rownames(g2),
               rownames(g3), rownames(g4),
               rownames(g5))),]














# Pvalue beta estimates
beta0 <- sol[sol$p.valor_beta0 <= 0.05, ,FALSE]
nrow(beta0)

# Pvalue beta estimates
same_but_change_in_time <- sol[(sol$p.valor_scmp_binned_pseudotime <= 0.05 &
                     sol$p.valor_scmp_binned_pseudotime2 <= 0.05 &
                         sol$p.valor_Path2vsPath1 > 0.05 &
                         sol$p.valor_scmp_binned_pseudotimexPath2 > 0.05 &
                         sol$p.valor_scmp_binned_pseudotime2xPath2 > 0.05), , drop =FALSE]

same_but_change_in_time_list <- list()
for (i in rownames(same_but_change_in_time)){
    same_but_change_in_time_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj, 
                                                       feature_id = i,
                                                       logs = T, logType = "log")
}
same_but_change_in_time_list[[sample(names(same_but_change_in_time_list), 1)]]
nrow(same_but_change_in_time)


# Pvalue beta estimates
path1_stable_path2_change <- sol[(sol$p.valor_scmp_binned_pseudotime > 0.05 &
                         sol$p.valor_scmp_binned_pseudotime2 > 0.05 &
                         sol$p.valor_Path2vsPath1 > 0.05 &
                         sol$p.valor_scmp_binned_pseudotimexPath2 <= 0.05 &
                         sol$p.valor_scmp_binned_pseudotime2xPath2 <= 0.05), , drop =FALSE]

path1_stable_path2_change_list <- list()
for (i in rownames(path1_stable_path2_change)){
    path1_stable_path2_change_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj, 
                                                      feature_id = i,
                                                      logs = T, logType = "log")
}
path1_stable_path2_change_list[[sample(names(path1_stable_path2_change_list), 1)]]


# Opposite change
test <- sol[(sol$p.valor_scmp_binned_pseudotime > 0.05 &
                                      sol$p.valor_scmp_binned_pseudotime2 > 0.05 &
                                      sol$p.valor_Path2vsPath1 > 0.05 &
                                      sol$p.valor_scmp_binned_pseudotimexPath2 <= 0.05 &
                                      sol$p.valor_scmp_binned_pseudotime2xPath2 <= 0.05), , drop =FALSE]

test_list <- list()
for (i in rownames(test)){
    test_list[[i]]<- sc.PlotGroups(scmpObj = scmp.obj, 
                                                        feature_id = i,
                                                        logs = T, logType = "log")
}
pathA_list[[sample(names(pathA_list), 1)]]

sol_temp <- sol[setdiff(rownames(sol), c(rownames(test), rownames(path1_stable_path2_change), rownames(same_but_change_in_time))),]

sc.PlotGroups(scmpObj = scmp.obj, 
              feature_id = "Gene293",
              logs = T, logType = "log")
