# Title: Run ScMaSigPro on simulated datasets with different Skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "benchmarks/11_Split_Bin_Parameter/data/simulated/sce/"
dir.create("benchmarks/11_Split_Bin_Parameter/data/output/", showWarnings = F)
helpScriptsDir <- "R_Scripts/helper_function/"

# Load names of files
dataSets <- list.files(paste0(dirPath))
names(dataSets) <- str_remove(
  str_split_i(dataSets, pattern = "_", i = 2),
  ".RData"
)

# Set-up a for loop
for (i in names(dataSets)) {
  poly.degree <- 2
  min.gene <- 6
  theta.val <- 1
  ep <- 0.00001
  
  cat(paste("\nRunning for skew:", i))
  
  stop()

  # Load Data
  load(file = paste0(dirPath, dataSets[i]))

  tryCatch(
    expr = {
      # Convert
        scmp.obj <- as_scmp(sim.sce, from = "sce",
                            additional_params = list(
                                labels_exist = TRUE,
                                existing_pseudotime_colname = "Step",
                                existing_path_colname = "Group"), verbose = F)

        # Compress
        scmp.obj.split.bin.true <- squeeze(
            scmpObject = scmp.obj,
            bin_method = "Sturges",
            drop.fac = 0.6,
            verbose = F,
            cluster_count_by = "sum",
            split_bins = TRUE,
            prune_bins = F,
            drop_trails = F,
            fill_gaps = F
        )
        sc.plot.bins.bar(scmp.obj.split.bin.true)
        
        # Split bin false
        scmp.obj.split.bin.false <- squeeze(
            scmpObject = scmp.obj,
            bin_method = "Sturges",
            drop.fac = 0.6,
            verbose = F,
            cluster_count_by = "sum",
            split_bins = FALSE,
            prune_bins = F,
            drop_trails = F,
            fill_gaps = F
        )
        sc.plot.bins.bar(scmp.obj.split.bin.false)

        # Make Design
  scmp.obj.split.bin.true <- sc.make.design.matrix(scmp.obj.split.bin.true,
                                          poly_degree = poly.degree)
  scmp.obj.split.bin.false <- sc.make.design.matrix(scmp.obj.split.bin.false,
                                          poly_degree = poly.degree)

        # Run p-vector
  scmp.obj.split.bin.true <- sc.p.vector(
            scmpObj = scmp.obj.split.bin.true, verbose = F, min.obs = 1,
            offset = T, parallel = T
        )
        scmp.obj.split.bin.false <- sc.p.vector(
            scmpObj = scmp.obj.split.bin.false, verbose = F, min.obs = 1,
            offset = T, parallel = T
        )
        
        # Run-Step-2
        scmp.obj.split.bin.true <- sc.T.fit(
            parallel = T,
            scmpObj = scmp.obj.split.bin.true, verbose = F,
            step.method = "backward",
            offset = T
        )
        scmp.obj.split.bin.false <- sc.T.fit(
            parallel = T,
            scmpObj = scmp.obj.split.bin.false, verbose = F,
            step.method = "backward",
            offset = T
        )

      # Save Object
      save(scmp.obj.split.bin.true, file = paste0("benchmarks/11_Split_Bin_Parameter/data/output/scmp.obj.skew.true.", i, ".RData"))
      save(scmp.obj.split.bin.false, file = paste0("benchmarks/11_Split_Bin_Parameter/data/output/scmp.obj.skew.false.", i, ".RData"))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i, "because", e$message))
    }
  )
}
