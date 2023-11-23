# Title: Run ScMaSigPro on simulated datasets with different levels of skewness
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
inPath <- "/supp_data/benchmarks/02_Skewness/simulated/sce/"
outPath <- "/supp_data/benchmarks/02_Skewness/output/"
outPath2 <- "Tables/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load helper functions
source(paste0(helpScriptsDir, "get_performance_ROCR_.R"))
source(paste0(helpScriptsDir, "calculate_metrics_binary.R"))

# Create Missing Directory
dir.create(outPath, showWarnings = F, recursive = T)

# Load names of files
dataSets <- list.files(paste0(inPath))
names(dataSets) <- str_remove(
    str_split_i(dataSets, pattern = "_", i = 2),
    ".RData"
)

# Zero-Inflation.evaluation
eval.list <- list()

# dataSets <- dataSets[names(dataSets) %in% c("90")]

# Set-up a for loop
for (i in names(dataSets)) {
    # Set variables
    poly.degree <- 2
    drop_fac <- 1
    maxit <- 100
    fam <- MASS::negative.binomial(10)
    
    cat(paste("\nRunning for Skewness:", i))
    
    # Load Data
    load(file = paste0(inPath, dataSets[i]))
    
    # Evaluate
    tryCatch(
        expr = {
            # Convert
            scmp.obj <- as_scmp(sim.sce,
                                from = "sce",
                                align_pseudotime = F,
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
                drop_fac = drop_fac,
                verbose = F,
                cluster_count_by = "sum",
                split_bins = F,
                prune_bins = F,
                drop_trails = F,
                fill_gaps = F
            )
            
            # Make Design
            scmp.obj <- sc.make.design.matrix(scmp.obj, poly_degree = poly.degree)
            
            # Run p-vector
            scmp.obj <- sc.p.vector(
                scmpObj = scmp.obj, verbose = F, min.obs = 1, parallel = F,
                offset = T,
                logOffset = F,
                useWeights = T,
                logWeights = F,
                useInverseWeights = T,
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
            save(scmp.obj, file = paste0(outPath, "scmp.obj.skew.", i, ".RData"))
            
            # Validate
            cat(paste("\nCompleted for", i))
            # Evaluate
            row_data <- as.data.frame(
                rowData(scmp.obj@sce)
            )[, c("gene_short_name", "status")]
            
            # Set binary labels
            gene.change <- rep(1, length(rownames(row_data[row_data$status != "No_Change", ])))
            gene.no.change <- rep(0, length(rownames(row_data[row_data$status == "No_Change", ])))
            
            # Add names
            names(gene.change) <- rownames(row_data[row_data$status != "No_Change", ])
            names(gene.no.change) <- rownames(row_data[row_data$status == "No_Change", ])
            
            # Ground truth
            groundTruth <- c(gene.change, gene.no.change)
            
            # Get Performance
            performance.measure <- as.data.frame(get_performance_ROCR(
                scmpObj = scmp.obj,
                groundTruth = groundTruth,
                r2_sequence = seq(0.00, 0.95, 0.05),
                include_influ = TRUE
            ))
            
            # Add to list
            performance.measure[["parameter"]] <- "Skew"
            performance.measure[["parameter.value"]] <- i
            eval.list[[i]] <- performance.measure
        },
        error = function(e) {
            cat(paste("\nFailed for", i, "because", e$message))
        }
    )
}

# Combine
evaluation.frame <- do.call(rbind, eval.list)

# Write
write.table(evaluation.frame, paste0(outPath2, "02_Skew_Performance.Table.tsv"),
            sep = "\t", row.names = F, quote = F
)
