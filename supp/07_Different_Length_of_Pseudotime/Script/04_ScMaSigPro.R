# Load Library
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(maSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtools))

# File-Location
inPath <- "supp/07_Different_Length_of_Pseudotime/data/input/PseudobulkInput/"
outPath <- "supp/07_Different_Length_of_Pseudotime/data/output/MaSigPro/"

# Load File names
fileNames <- list.files(inPath)
names(fileNames) <- fileNames

# Some Variables
poly.degree <- 2
min.gene <- 10
theta.val <- 1
ep <- 0.00001

# Run for Loop one-by-one
for (i in names(fileNames)) {
  # Load RDS
  input <- readRDS(paste0(inPath, i))

  # pbCounts
  pbCounts <- input$pbCounts
  pbMeta <- input$pbMeta

  # Check
  if (all(rownames(pbMeta) != colnames(pbCounts))) {
    stop("Mismatch")
  }

  # Run P-vector and T- step
  tryCatch(
    expr = {
      # Create MaSigPro Design File
      design <- make.design.matrix(
        edesign = as.data.frame(pbMeta), # bulkMeta is edesign
        degree = poly.degree,
        time.col = 1, repl.col = 2
      )

      # Run P-Vector
      gc <- capture.output(p.vector.fit <- p.vector(
        data = pbCounts, design = design, Q = 0.05,
        MT.adjust = "BH", min.obs = 6, counts = TRUE,
        theta = theta.val, epsilon = ep
      ))

      # Runnig T Step
      gc <- capture.output(tstep.fit <- T.fit(p.vector.fit,
        step.method = "backward",
        family = p.vector.fit$family,
        epsilon = ep
      ))
      # Create Out Dire
      dir.create(paste0(outPath, "pvector"), showWarnings = F, recursive = T)
      dir.create(paste0(outPath, "tstep"), showWarnings = F, recursive = T)

      # Save
      saveRDS(p.vector.fit, file = paste0(outPath, "pvector/", i))
      saveRDS(tstep.fit, file = paste0(outPath, "tstep/", i))

      # Validate
      cat(paste("\nCompleted for", i))
    },
    error = function(e) {
      cat(paste("\nFailed for", i))
    }
  )
}
