
## Plot elbow
ElbowPlot(sob_prs, ndims = 200)

## Get Variance Explained by top 50 components
mat <- GetAssayData(sob_prs, layer = "RNA", slot = "scale.data")
pca <- sob_prs[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
sum(varExplained[c(1:50)] * 100)

# Create Graph
sob_prs <- FindNeighbors(sob_prs, dims = 1:50)
