##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: scMaSigPro Application ########################
##########################################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(Matrix))
<<<<<<< HEAD
suppressPackageStartupMessages(library(Azimuth))
=======
>>>>>>> dev
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))

# Load scMaSigPro
<<<<<<< HEAD
install.packages("../scMaSigPro_1.0.1.tar.gz")
suppressPackageStartupMessages(library(scMaSigPro))


# Prefix
prefixIn <- "benchmarks/11_RealDataSmall/data/results/"
prefixOut <- "benchmarks/11_RealDataSmall/data/output/"

# Load CDS object
load(paste0(prefixIn, "monocle3_inferred_pseudotime.RData"))

# Monocl3 3 object
cds

# Convert the ScMaSigPro Object
scmp.obj <- as_scmp(cds, from = "cell_data_set")

# Plot the Paths
Before <- plot_cells(scmp.obj@sce, color_cells_by = "predicted.celltype.l2")

# Select Path
scmp.obj <- selectPath(obj = scmp.obj, sel.path = c("Path1", "Path2"),
                       balance_paths = T, pathCol = "Path",
                       pTimeCol = "Pseudotime", plot_paths = T, verbose = F)

# Plot the Paths
after<- plot_cells(scmp.obj@sce, color_cells_by = "Path", cell_size = 2)
after.cell<- plot_cells(scmp.obj@sce, color_cells_by = "Path",cell_size = 2) + theme(legend.position = "bottom")
after+after.cell

# Compress
scmp.obj <- squeeze(
    scmp.ob = scmp.obj,
    time.col = "Pseudotime",
    path.col = "Path",
    method = "Sturges",
    drop.fac = 0.7,
    verbose = T,
    cluster.count.by = "sum"
)

# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
                                  degree = 2,
                                  time.col = "binnedTime",
                                  path.col = "path"
)

# Run p-vector
scmp.obj <- sc.p.vector(
    scmpObj = scmp.obj, verbose = T, min.obs = 10,
    counts = T, theta = 10,
    offset = T
)

# Run-Step-2
scmp.obj <- sc.T.fit(
    data = scmp.obj, verbose = T,
    step.method = "backward",
    family = scmp.obj@scPVector@family,
    offset = T
)

# Save the object for further analysis
save(scmp.obj, file = paste0(prefixOut, "scmp.obj.RData"))
=======
suppressPackageStartupMessages(library(scMaSigPro))

# Monocl3 3 object
cds <- readRDS("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/Monocle3_Processed_Donor3.RDS")

# Convert the ScMaSigPro Object
scmp.obj <- as_scmp(cds, from = "cds",
                    align_pseudotime = F,
                    annotation_colname = "cell_type")

# Compress
scmp.obj <- squeeze(scmpObject = scmp.obj,
                    split_bins = F,
                    prune_bins = F,
                    drop_trails = F)

a <- sc.plot.bins.bar(scmp.obj)
b <- sc.plot.bins.tile(scmp.obj)
a +b


# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
                                  poly_degree = 3,
)

scmp.obj@distribution <- MASS::negative.binomial(10)

# Run p-vector
scmp.obj <- sc.p.vector(
    parallel = T,
    scmpObj = scmp.obj, verbose = T,
    max_it = 10000,
    globalTheta = T,
    logOffset = T,
    useInverseWeights = F,
    logWeights = F,
    offset = T,
    min.obs = 1)

# Run-Step-2
scmp.obj <- sc.T.fit(
    scmpObj = scmp.obj, verbose = T,
    step.method = "backward"
)

# Save the object for further analysis
save(scmp.obj, file = paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/scMaSigPro_Donor2_Monocle3.RData"))
>>>>>>> dev
