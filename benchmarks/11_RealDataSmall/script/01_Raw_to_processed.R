#######################################
## Author: Priyansh Srivastava ########
## Email: spriyansh29@gmail.com #######
## Script: Pre-Process and CC Remove ##
#######################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))

num_cores <- 16

# Prefix
prefixIn <- "benchmarks/11_RealDataSmall/data/input/"
prefixOut <- "benchmarks/11_RealDataSmall/data/results/"

# Read BioMart info
biomart.anno <- readRDS(paste0(prefixIn, "cell_cycle_data.mart"))
reg.out <- c(unique(biomart.anno$SYMBOL))

# Set Variables to be used later
i <- "rep1"
individual <- "1"
age <- "35"
sex <- "Male"

# ggplot2 theme definition
custom_theme_1 <- ggplot() +
  theme_minimal() +
  theme(
    title = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.position = "none"
  ) +
  xlab(paste0(
    "Individual: ", individual,
    " Age: ", age,
    " Sex: ", sex
  ))

# Load filetered Matrix from cell ranger
filt_mat <- Read10X_h5(
  filename = paste0(
    prefixIn,
    "rep1_p1_downSampled/outs/filtered_feature_bc_matrix.h5"
  )
)

# Renaming the features
filt_mat@Dimnames[[1]] <- str_replace(filt_mat@Dimnames[[1]], "_", "-")

# Create Seurat Object
sob.raw <- CreateSeuratObject(
  counts = filt_mat, min.cells = 200,
  min.features = 1000, project = i
)

# Calculate Percentage of Mitochondrial Read
sob.raw[["percent.mt"]] <- PercentageFeatureSet(sob.raw, pattern = "^MT-")

# Violin Plots
p <- VlnPlot(sob.raw,
  features = c("nFeature_RNA"),
  pt.size = 0.5
) + theme_minimal() + xlab(paste0(
  "Individual: ", individual,
  " Age: ", age,
  " Sex: ", sex
))

ggsave(p,
  filename = paste0(prefixOut, i, "_Raw_VlnPlot_nFeature_RNA.png"),
  dpi = 1200, limitsize = FALSE
)


p <- VlnPlot(sob.raw,
  features = c("nCount_RNA"),
  pt.size = 0.5
) + theme_minimal() + xlab(paste0(
  "Individual: ", individual,
  " Age: ", age,
  " Sex: ", sex
))
ggsave(p,
  filename = paste0(prefixOut, i, "_Raw_VlnPlot_nCount_RNA.png"),
  dpi = 1200, limitsize = FALSE
)

p <- VlnPlot(sob.raw,
  features = c("percent.mt"),
  pt.size = 0.5
) + theme_minimal() + xlab(paste0(
  "Individual: ", individual,
  " Age: ", age,
  " Sex: ", sex
))
ggsave(p,
  filename = paste0(prefixOut, i, "_Raw_VlnPlot_percent.mt.png"),
  dpi = 1200, limitsize = FALSE
)

p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_smooth(method = "lm", formula = "y~x") +
  theme_minimal() + xlab(paste0(
    "Individual: ", individual,
    " Age: ", age,
    " Sex: ", sex
  ))
ggsave(p,
  filename = paste0(prefixOut, i, "_Raw_FeatureScatter_percent.mt.png"),
  dpi = 1200, limitsize = FALSE
)

p <- FeatureScatter(sob.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm", formula = "y~x") +
  theme_minimal() + xlab(paste0(
    "Individual: ", individual,
    " Age: ", age,
    " Sex: ", sex
  ))
ggsave(p,
  filename = paste0(prefixOut, i, "_Raw_FeatureScatter_nFeature_RNA.png"),
  dpi = 1200, limitsize = FALSE
)

# Subset # Check Required if-else
sob.sub <- subset(sob.raw, subset = nFeature_RNA > 200 & nCount_RNA < 30000 & percent.mt < 5)

# Add Measuremnets
file_name <- paste0(prefixOut, i, "_sub_sob")
SaveH5Seurat(
  object = sob.sub, filename = file_name, overwrite = T,
  verbose = FALSE
)

# Plots
p <- VlnPlot(sob.sub, features = c("nFeature_RNA")) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = "y~x") +
  xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
ggsave(p,
  filename = paste0(prefixOut, i, "_Sub_VlnPlot_nFeature_RNA.png"),
  dpi = 1200, limitsize = FALSE
)

p <- VlnPlot(sob.sub, features = c("nCount_RNA")) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = "y~x") +
  xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
ggsave(p,
  filename = paste0(prefixOut, i, "_Sub_VlnPlot_nCount_RNA.png"),
  dpi = 1200, limitsize = FALSE
)

p <- VlnPlot(sob.sub, features = c("percent.mt")) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = "y~x") +
  xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
ggsave(p,
  filename = paste0(prefixOut, i, "_Sub_VlnPlot_percent.mt.png"),
  dpi = 1200, limitsize = FALSE
)

p <- FeatureScatter(sob.sub, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  theme_minimal() +
  geom_smooth(method = "lm", formula = "y~x") +
  xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
ggsave(p,
  filename = paste0(prefixOut, i, "_Sub_FeatureScatter_percent.mt.png"),
  dpi = 1200, limitsize = FALSE
)

p <- FeatureScatter(sob.sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme_minimal() +
  geom_smooth(method = "lm", formula = "y~x") +
  xlab(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
ggsave(p,
  filename = paste0(prefixOut, i, "_Sub_FeatureScatter_nFeature_RNA.png"),
  dpi = 1200, limitsize = FALSE
)

# Pre-processing
sob.prs <- NormalizeData(sob.sub, verbose = F)
sob.prs <- FindVariableFeatures(sob.prs, selection.method = "vst", verbose = F)
sob.prs <- ScaleData(sob.prs, features = rownames(sob.prs), verbose = F)
sob.prs <- RunPCA(sob.prs,
  features = VariableFeatures(sob.prs),
  verbose = F, ndims.print = 0, nfeatures.print = 0
)

# Basic UMAP
p <- DimPlot(sob.prs, reduction = "pca", pt.size = 0.5) + theme_minimal() +
  theme(legend.position = "none") +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex))
ggsave(p,
  filename = paste0(prefixOut, i, "_basic_PCA.png"),
  dpi = 1200, limitsize = FALSE
)

# Save
file_name <- paste0(prefixOut, i, "_prep_sob")
SaveH5Seurat(
  object = sob.prs, filename = file_name, overwrite = T,
  verbose = FALSE
)

indata <- rownames(sob.prs)[rownames(sob.prs) %in% reg.out]

# Keep only those which are present in data
biomart.anno <- biomart.anno[biomart.anno$SYMBOL %in% indata, ]

# For seurat we need to divide genes into vectors
s.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0006260"), "SYMBOL"])
g2m.genes <- unique(biomart.anno[biomart.anno$GO %in% c("GO:0000087", "GO:0000279", "GO:0007059", "GO:0048285"), "SYMBOL"])

# Cell Cycle Scores
sob.cc <- CellCycleScoring(sob.prs,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# Re-Run PCA
sob.cc <- RunPCA(sob.cc,
  features = c(s.genes, g2m.genes),
  verbose = F, ndims.print = 0, nfeatures.print = 0
)

# Save The plot
p <- DimPlot(sob.cc, pt.size = 0.5) +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_manual(
    name = "Cell-Cycle GOs",
    breaks = c("G2M", "S", "G1"),
    values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
  )
ggsave(p,
  filename = paste0(prefixOut, i, "_CCG_PCA.png"),
  dpi = 1200, limitsize = FALSE
)

# Scale to Regress out
sob.cc <- ScaleData(sob.cc,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(sob.cc), verbose = F
)

# Re-Run PCA
test <- RunPCA(sob.cc,
  features = c(s.genes, g2m.genes),
  verbose = F, ndims.print = 0, nfeatures.print = 0
)
p <- DimPlot(test, reduction = "pca", pt.size = 0.5) +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_manual(
    name = "Cell-Cycle GOs",
    breaks = c("G2M", "S", "G1"),
    values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
  )
ggsave(p,
  filename = paste0(prefixOut, i, "_CCC_PCA.png"),
  dpi = 1200, limitsize = FALSE
)

sob.cc <- RunPCA(sob.cc,
  features = VariableFeatures(sob.cc),
  nfeatures.print = 0, verbose = F, ndims.print = 0
)

p <- DimPlot(sob.cc, reduction = "pca", pt.size = 0.5) +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_manual(
    name = "Cell-Cycle GOs",
    breaks = c("G2M", "S", "G1"),
    values = c("G2M" = "#009E73", "S" = "#D55E00", "G1" = "#56B4E9")
  )

ggsave(p,
  filename = paste0(prefixOut, i, "_CCC_all_PCA.png"),
  dpi = 1200, limitsize = FALSE
)

# Save
file_name <- paste0(prefixOut, i, "_CCC_sob")
SaveH5Seurat(
  object = sob.cc, filename = file_name, overwrite = T,
  verbose = FALSE
)

# More Pre-processing
sob.p <- FindNeighbors(sob.cc, verbose = F)
sob.p <- FindClusters(sob.p, verbose = F)
sob.p <- RunUMAP(sob.p, dims = 1:50, verbose = F)

cat(paste0("\nCheck-8 (", i, "): Clustering Completed\n"))

file_name <- paste0(prefixOut, i, "_Processed_sob")
SaveH5Seurat(
  object = sob.p, filename = file_name,
  overwrite = T, verbose = FALSE
)

# Save the UMAP with the cluster
p <- DimPlot(sob.p, reduction = "umap", pt.size = 0.5, group.by = "seurat_clusters") +
  ggtitle(paste0("Individual: ", individual, " Age: ", age, " Sex: ", sex)) +
  scale_color_hue(l = 50)
ggsave(p,
  filename = paste0(prefixOut, i, "_C_UMAP.png"),
  dpi = 1200, limitsize = FALSE
)
