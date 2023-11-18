FQnorm <- function(counts) {
  rk <- apply(counts, 2, rank, ties.method = "min")
  counts.sort <- apply(counts, 2, sort)
  refdist <- apply(counts.sort, 1, median)
  norm <- apply(rk, 2, function(r) {
    refdist[r]
  })
  rownames(norm) <- rownames(counts)
  return(norm)
}



# Seurat normCounts
calcNormCounts <- function(rawCounts, cat, size_fac = 10000) {
  # Base
  suppressPackageStartupMessages(require(Seurat))
  suppressPackageStartupMessages(require(sctransform))
  suppressPackageStartupMessages(require(Matrix))
  FQnorm <- function(counts) {
    rk <- apply(counts, 2, rank, ties.method = "min")
    counts.sort <- apply(counts, 2, sort)
    refdist <- apply(counts.sort, 1, median)
    norm <- apply(rk, 2, function(r) {
      refdist[r]
    })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

  # Make Seurat Object
  seuratObject <- CreateSeuratObject(counts = rawCounts)
  
  ## Lib.size
  if (cat == "libSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject,
      normalization.method = "RC",
      scale.factor = size_fac
    )
    return(seuratNormObject)
    # Extract Counts
    seuratNormCounts <- seuratNormObject@assays$RNA@data
    # Return dgCMatrix
    return(seuratNormCounts)
  }

  ## Log + Lib.size
  else if (cat == "logLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "LogNormalize"
    )
    # Extract Counts
    seuratNormCounts <- seuratNormObject@assays$RNA@data
    # Return dgCMatrix
    return(seuratNormCounts)
  }

  ## Centered + Lib.size
  else if (cat == "cenLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "RC"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$RNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = F)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "scLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "RC"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$RNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = T)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  }
  ## Centered + log + Lib.size
  else if (cat == "cenLogLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "LogNormalize"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$RNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = F)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "scLogLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "LogNormalize"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$RNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = T)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "sctransform") {
    # Normalize
    seuratNormObject <- SCTransform(seuratObject,
      do.center = F,
      do.scale = F,
      verbose = FALSE
    )
    # Extract Counts
    seuratNormCounts <- as.matrix(seuratNormObject@assays$SCT@data)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "cenSctransform") {
    # Normalize
    seuratNormObject <- SCTransform(seuratObject,
      do.center = T,
      do.scale = F,
      verbose = FALSE
    )
    # Extract Counts
    seuratNormCounts <- as.matrix(seuratNormObject@assays$SCT@data)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "scSctransform") {
    # Normalize
    seuratNormObject <- SCTransform(seuratObject,
      do.center = T,
      do.scale = T,
      verbose = FALSE
    )
    # Extract Counts
    seuratNormCounts <- as.matrix(seuratNormObject@assays$SCT@data)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "FQNorm") {
    # Calculate FQ norm
    fqNormCounts <- FQnorm(rawCounts)
    return(fqNormCounts)
  } else if (cat == "CLR") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject,
      normalization.method = "CLR",
      scale.factor = size_fac
    )
    # Extract Counts
    seuratNormCounts <- seuratNormObject@assays$RNA@data
    # Return dgCMatrix
    return(seuratNormCounts)
  } else {
    stop("Please use one of the availble methods")
  }
}
