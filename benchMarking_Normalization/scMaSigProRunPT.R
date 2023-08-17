# Load Library
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(maSigPro))
suppressPackageStartupMessages(library(gtools))

# Load Prefix
parentFix <- "benchMarking_Normalization/"
scriptFix <- paste0("old_Scripts/")
inPrefix <- paste0(parentFix, "data/input/PB_Input/")
outPrefix <- paste0(parentFix, "data/output/scMaSigPro/")

# Load Custim Functions
load <- function() {
  source(paste0(scriptFix, "plot_simulations().R"))
  source(paste0(scriptFix, "add_gene_anno().R"))
  source(paste0(scriptFix, "calc_bin_size.R"))
  source(paste0(scriptFix, "calcNormCounts.R"))
  source(paste0(scriptFix, "entropy_discretize.R"))
  source(paste0(scriptFix, "make_bulk_design.R"))
  source(paste0(scriptFix, "make_bulk_counts.R"))
  source(paste0(scriptFix, "create_range.R"))
}

# 42: zi_shape_-0.2_
# 52: zi_shape_0_
# 62: zi_shape_0.2_
# 87: zi_shape_1_

# Load File names
all_file_names <- list.files(inPrefix)
all_file_names2 <- all_file_names[grepl(x = all_file_names, "zi_shape_1_")]
all_file_names <- all_file_names[grepl(x = all_file_names, "zi_shape_1_")]
all_file_names <- str_remove(all_file_names, ".RDS")
all_file_names <- str_split(all_file_names, "_")
names(all_file_names) <- all_file_names2
#all_file_names <- all_file_names[names(all_file_names) == "zi_shape_-0.2_TrueCounts.RDS"]


# Some Variables
poly.degree <- 2
min.gene <- 2
theta.val <- 10
ep <- 0.00001


# Run for Loop one-by-one
for (i in names(all_file_names)) {
  print(paste0("Running: ", i))

  # Set Variables
  list_name <- i

  # Read data list
  data.list <- readRDS(paste0(inPrefix, list_name))

  # Extract the design file
  exp.design.file <- data.list[["cell_meta"]]

  # Create List of all the counts
  count_list <- list(
    sum = data.list[["count_table"]][["sum"]],
    mean = data.list[["count_table"]][["mean"]],
    median = data.list[["count_table"]][["median"]]
  )

  for (j in names(count_list)) {
      
      stop()

    # Extract the Counts
    ct <- count_list[[j]]

    if (all(colnames(ct) == rownames(exp.design.file)) == F) {
      stop("MisMatch")
    }

    # Run P-vector and T- step
    tryCatch(
      expr = {
        # Create MaSigPro Design File
        design <- make.design.matrix(
          edesign = as.data.frame(exp.design.file), # bulkMeta is edesign
          degree = poly.degree,
          time.col = 1, repl.col = 2
        )

        # Run P-Vector
        gc <- capture.output(p.vector.fit <- p.vector_test(
          data = ct, design = design, Q = 0.05,
          MT.adjust = "BH", min.obs = min.gene, counts = TRUE,
          theta = theta.val, epsilon = ep
        ))

        # Runnig T Step
        gc <- capture.output(tstep.fit <- T.fit(p.vector.fit,
          step.method = "backward",
          family = p.vector.fit$family,
          epsilon = ep
        ))
        
        print(paste("Completed", j))

        # Create the directory
        fName <- paste0(outPrefix, str_remove(i, ".RDS"), "_", j, "_", "tstep", ".RDS")
        saveRDS(tstep.fit, fName)
        x <- gc()
        x <- NULL
      },
      error = function(e) {
        print("Failed")
      }
    )
  }
}

cat(paste0("\n", "Script Finished"))
