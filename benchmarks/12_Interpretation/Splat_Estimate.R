library(splatter)
library(Seurat)
library(coop)
helpScriptsDir <- "R_Scripts/helper_function/"
source(paste0(helpScriptsDir, "plot_simulations().R"))
source(paste0(helpScriptsDir, "add_gene_anno().R"))
source(paste0(helpScriptsDir, "calc_bin_size.R"))

realData <- Read10X_h5("/supp_data/tmp/filtered_feature_bc_matrix.h5")
realData <- as.matrix(realData)
params <- splatEstimate(realData)

realData <- ReadMtx(mtx = "/supp_data/tmp/E-MTAB-3321.aggregated_filtered_counts.mtx",
                    cells = "/supp_data/tmp/E-MTAB-3321.aggregated_filtered_counts.mtx_cols",
                    features = "/supp_data/tmp/E-MTAB-3321.aggregated_filtered_counts.mtx_rows")


# Create Raw Seurat Object
sob.raw <- CreateSeuratObject(
    counts = realData, min.cells = 1500,
    min.features = 100)

# Calculate Percentage of Mitochondrial Reads
sob.raw[["percent.mt"]] <- PercentageFeatureSet(sob.raw, pattern = "^MT-")

# Remove Cells 
sob.sub <- subset(sob.raw, subset = nFeature_RNA > 100 & nCount_RNA < 40000 & percent.mt < 5)
raw_counts <- sob.sub@assays$RNA@counts

#c(sample(nrow(realData), size = 5000, replace = F)),
set.seed(123)
raw_counts <- raw_counts[,c(sample(ncol(raw_counts), size = 3000, replace = F))]
round(sparsity(as.matrix(realData)) * 100)

params <- splatEstimate(as.matrix(raw_counts))

params <- setParams(params,
                    update = list(
                        dropout.type = "experiment",
                        group.prob = c(0.5, 0.5),
                        path.from = c(0, 0),
                        path.nSteps = c(1500, 1500),
                        batchCells = 3000, # Number of Cells
                        nGenes = 5000, # Number of Genes
                        seed = 2022
                    ))

sim.sce <- splatSimulate(params = params,
                         method = "paths")

# Proportion of true Sparsity
trueSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$TrueCounts)) * 100)
simulatedSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100) - trueSparsity
totSparsity <- round(sparsity(as.matrix(sim.sce@assays@data@listData$counts)) * 100)

cat(paste("\nTotal:",totSparsity))
cat(paste("\nsimulatedSparsity:", simulatedSparsity))
cat(paste("\ntrueSparsity:", trueSparsity))


plot_simulations(sim.sce, assay_type = "counts", plot3d = F, plot2d = T,
                 colorGroupContinuous = "Step", colorGroupDiscrete = "Group",
                 frame = )
