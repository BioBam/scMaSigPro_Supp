# Title: Run ScMaSigPro on simulated datasets with different levels of sparsity
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/01_Sparsity/data/simulated/sce/"
dir.create("benchmarks/01_Sparsity/data/output/", showWarnings = F)
helpScriptsDir <- "R_Scripts/helper_function/"

# Load names of files
dataSets <- list.files(paste0(dirPath))
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = "_", i = 2),
  ".RData"
)

dataSets <- dataSets[names(dataSets) %in% c("70", "80", "90")]

# Set-up a for loop
for (i in names(dataSets)) {
  poly.degree <- 2
  min.gene <- 6
  drop_fac <- 1
  gTheta<- FALSE
  off <- T
  
  cat(paste("\nRunning for sparsity:", i))

  # Load Data
  load(file = paste0(dirPath, dataSets[i]))
  
  if (i %in% c("70", "80", "90")) {
      drop_fac <- 0.3
      
      # Extract the counts
      raw.counts <- as.matrix(sim.sce@assays@data@listData$counts) 
      # Count the number of zeros for each gene
      zero_counts_per_gene <- apply(raw.counts, 1, function(x) sum(x == 0))
      gene_names <- rownames(raw.counts)
      zero_count_table <- data.frame(Gene = gene_names, ZeroCounts = zero_counts_per_gene)
      keep_genes <- rownames(zero_count_table[zero_count_table$ZeroCounts< 2800,])
      sim.sce <- sim.sce[keep_genes,]
      
  }

  tryCatch(
    expr = {
      # Convert
      scmp.obj <- as_scmp(sim.sce, from = "sce",
                          additional_params = list(
                              labels_exist = TRUE,
                              existing_pseudotime_colname = "Step",
                              existing_path_colname = "Group"), verbose = F)
      
      if (i %in% c("80", "90")) {
          drop_fac <- 0.3
          maxit <- 10000
          gTheta <- T
      }
      
      if (i %in% c("70")) {
          maxit <- 10000
      }

      # Compress
      scmp.obj <- squeeze(
        scmpObject = scmp.obj,
        bin_method = "Sturges",
        drop.fac = drop_fac,
        verbose = F,
        cluster_count_by = "sum",
        split_bins = T,
        prune_bins = F,
        drop_trails = F,
        fill_gaps = F
      )
      

      # Make Design
      scmp.obj <- sc.make.design.matrix(scmp.obj,
        poly_degree = poly.degree)

      # Run p-vector
      scmp.obj <- sc.p.vector(
          scmpObj = scmp.obj, verbose = F, min.obs = 1,
          offset = T, parallel = F, useWeights = T,
          logWeights = F, logOffset = F,
          max_it = maxit, 
          useInverseWeights = F, globalTheta = gTheta
      )

      # Run-Step-2
      scmp.obj <- sc.T.fit(
          parallel = T,
        scmpObj = scmp.obj, verbose = F,
        step.method = "backward",
        offset = T
      )

      # Save Object
      save(scmp.obj, file = paste0("benchmarks/01_Sparsity/data/output/scmp.obj.sparsity.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}


# # Explore
# 
# gene.metadata <- rowData(sim.sce) %>% as.data.frame()
# raw.counts <- as.matrix(sim.sce@assays@data@listData$counts)%>% as.data.frame()
# gene.of.interest <- as.matrix(raw.counts[raw.counts$gene_short_name %in% c("Gene3087", "Gene3088"), ])
# 
# 
# # Count the number of zeros for each gene
# zero_counts_per_gene <- apply(raw.counts, 1, function(x) sum(x == 0))
# 
# # Create a data frame with gene names and their zero counts
# gene_names <- rownames(raw.counts)
# zero_count_table <- data.frame(Gene = gene_names, ZeroCounts = zero_counts_per_gene)
# 
# # Look at the table
# print(zero_count_table)
# 
# # If you want to order the table based on the number of zeros
# zero_count_table_ordered <- zero_count_table[order(zero_count_table$ZeroCounts, decreasing = TRUE), ]
# print(zero_count_table_ordered)
# 
