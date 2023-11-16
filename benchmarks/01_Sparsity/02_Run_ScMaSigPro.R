# Title: Run ScMaSigPro on simulated datasets with different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/01_Sparsity/simulated/sce/"
outPath <- "/supp_data/benchmarks/01_Sparsity/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = "_", i = 2),
  ".RData"
)

#dataSets <- dataSets[names(dataSets) %in% c("90")]

# Set-up a for loop
for (i in names(dataSets)) {
  poly.degree <- 2
  drop_fac <- 1
  gTheta<- FALSE
  maxit<- 100
  fam <- MASS::negative.binomial(10)
  split.bins = F
  
  cat(paste("\nRunning for sparsity:", i))

  # Load Data
  load(file = paste0(inPath, dataSets[i]))
  
  if (i %in% c("70", "80", "90")) {
      
      # Extract the counts
      raw.counts <- as.matrix(sim.sce@assays@data@listData$counts) 
      # Count the number of zeros for each gene
      zero_counts_per_gene <- apply(raw.counts, 1, function(x) sum(x == 0))
      gene_names <- rownames(raw.counts)
      zero_count_table <- data.frame(Gene = gene_names, ZeroCounts = zero_counts_per_gene)
      keep_genes <- rownames(zero_count_table[zero_count_table$ZeroCounts< 2800,])
      sim.sce <- sim.sce[keep_genes,]
  }
  if (i %in% c("90", "80")){
      
      # Additinal Compression
      drop_fac = 0.3
      poly.degree <- 1
      split.bins = T
      fam =  poisson()
  }
  # 
  # if (i %in% c("70")) {
  #     
  #     # More Compression
  #     drop_fac <- 0.5
  #     
  #     # Increase number of iteration
  #     maxit <- 10000
  # }

  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as_scmp(sim.sce, from = "sce",
                          additional_params = list(
                              labels_exist = TRUE,
                              existing_pseudotime_colname = "Step",
                              existing_path_colname = "Group"), verbose = F)
      
      # Compress
      scmp.obj <- squeeze(
        scmpObject = scmp.obj,
        bin_method = "Sturges",
        drop.fac = drop_fac,
        verbose = F,
        cluster_count_by = "sum",
        split_bins = split.bins,
        prune_bins = F,
        drop_trails = F,
        fill_gaps = F
      )
      
      # Make Design
      scmp.obj <- sc.make.design.matrix(scmp.obj,
        poly_degree = poly.degree)

      # Run p-vector
      scmp.obj <- sc.p.vector(
          scmpObj = scmp.obj, verbose = F, min.obs = 1, parallel = T,
          offset = T, 
          logOffset = F,
          useWeights = T, 
          logWeights = F, 
          useInverseWeights = F,
          max_it = maxit, 
          family = fam
      )

      # Run-Step-2
      scmp.obj <- sc.T.fit(
          parallel = T,
          scmpObj = scmp.obj, verbose = F,
          step.method = "backward"
          )

      # Save Object
      save(scmp.obj, file = paste0(outPath, "scmp.obj.sparsity.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}
