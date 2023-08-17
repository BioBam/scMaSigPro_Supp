# Load EdgeR
suppressPackageStartupMessages(require(edgeR))

# Load ScSimulated Data
sce.sim <-  readRDS("simulatedData/Splatter/datasets/SCE_Objects/sim01_Non_linear.sce.RDS")

# Counts extraction
counts_raw <- as.matrix(sce.sim@assays@data@listData$counts)

# Extract Metadata
group_data <- as.data.frame(colData(sce.sim))
group_data_vector <- group_data$Group

# Create EdgeR Object
d <- DGEList(counts=counts_raw,
             group=factor(group_data_vector))

# Check Number of Genes
dim(d)

# CPM Filter
keep <- rowSums(edgeR::cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

# Refixing libraries
d$samples$lib.size <- colSums(d$counts)
d$samples

# Data Normalization
d <- calcNormFactors(d)
d

# MDS plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
d1 <- estimateTrendedDisp(d1)
d1 <- estimateDisp(d1)
names(d1)
design.mat <- model.matrix(~ 0 + d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
