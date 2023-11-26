###########################################
## Author: Priyansh Srivastava ############
## Email: spriyansh29@gmail.com ###########
## Script: Subsampling ####################
###########################################

set.seed(007)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec

# Load data
azimuth.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
  # Step-1: Add Annotation for donors
  if (rep_i == "rep1") {
    individual <- "Donor-1"
    age <- "35"
    sex <- "Male"
  } else if (rep_i == "rep2") {
    individual <- "Donor-2"
    age <- "28"
    sex <- "Female"
  } else if (rep_i == "rep3") {
    individual <- "Donor-3"
    age <- "19"
    sex <- "Female"
  }
  # Load seurat object
  sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_azimuth.RDS"))

  # Get all cells
  all_cells <- unique(sob@meta.data$cell_type)

  # Drop cells
  keep <- all_cells[!(all_cells %in% c(
    "BaEoMa", "Stromal", "transitional B"
  ))]

  # Subset
  sob.sub <- subset(sob, cell_type %in% keep)

  # Recompute
  sob.sub <- RunPCA(sob.sub, features = VariableFeatures(object = sob.sub), verbose = F)
  sob.sub <- FindNeighbors(sob.sub, verbose = F)
  sob.sub <- FindClusters(sob.sub, resolution = 1, verbose = F)

  # Compute UMAP
  sob.sub <- RunUMAP(sob.sub,
    verbose = F, features = VariableFeatures(sob.sub)
  )

  # Return
  return(sob.sub)
})


# Subsample for donor-1
all.donor.subSample.list <- lapply(
  names(azimuth.list),
  function(rep_i, inPath = dirPath, outPath = dirPath) {
      
      sob <- azimuth.list[[rep_i]]
      
      
      # Step-1: Add Annotation for donors
      if (rep_i == "rep1") {
          individual <- "Donor-1"
          age <- "35"
          sex <- "Male"
      } else if (rep_i == "rep2") {
          individual <- "Donor-2"
          age <- "28"
          sex <- "Female"
      } else if (rep_i == "rep3") {
          individual <- "Donor-3"
          age <- "19"
          sex <- "Female"
      }

    # Create lineage sets
    lineage_set <- list(
        hsc_lineage = c("HSC", "CLP", "EMP"),
      lymphoid_lineage = c("CLP", "LMPP", "NK", "pro B", "pre B", "Memory B"),
      granulocyte_monocyte_lineage = c("HSC","GMP", "CD14 Mono", "Macrophage", "cDC2", "pre-mDC", "pre-pDC"),
      megakaryocyte_erythroid_lineage = c("EMP", "Prog Mk", "Platelet", "Early Eryth", "Late Eryth")
    )

    # Get all cells
    all_cells <- unique(sob@meta.data$cell_type)

    # Create per linegae Plot
    sub.sample.list <- lapply(names(lineage_set), function(i) {
      # Keep cells
      keep <- all_cells[(all_cells %in% lineage_set[[i]])]

      # Subset
      sob.sub <- subset(sob, cell_type %in% keep)
      sob.sub <- subset(sob.sub, subset = nFeature_RNA > 100)

      # Recompute
      sob.sub <- RunPCA(sob.sub, features = VariableFeatures(object = sob.sub), verbose = F)
      sob.sub <- FindNeighbors(sob.sub, verbose = F)
      sob.sub <- FindClusters(sob.sub, resolution = 1, verbose = F)

      # Compute UMAP
      sob.sub <- RunUMAP(sob.sub,
        verbose = F, features = VariableFeatures(sob.sub),
        min.dist = 0.1
      )

      # Plot
      plt <- DimPlot(sob.sub, group.by = "cell_type") + ggtitle(paste(
        individual, "| Age:", age,
        "| sex:", sex
      ), subtitle = paste(i))


      # Return
      return(list(
        obj = sob.sub,
        plot = plt
      ))
    })

    # Set name
    names(sub.sample.list) <- names(lineage_set)
    
    # Return UMAP
    return(sub.sample.list)
  }
)

names(all.donor.subSample.list) <-  names(azimuth.list)

sub_samples <- ggarrange(
    all.donor.subSample.list$rep1$hsc_lineage$plot,
    all.donor.subSample.list$rep1$lymphoid_lineage$plot,
    all.donor.subSample.list$rep1$granulocyte_monocyte_lineage$plot,
    all.donor.subSample.list$rep1$megakaryocyte_erythroid_lineage$plot,
    all.donor.subSample.list$rep2$hsc_lineage$plot,
    all.donor.subSample.list$rep2$lymphoid_lineage$plot,
    all.donor.subSample.list$rep2$granulocyte_monocyte_lineage$plot,
    all.donor.subSample.list$rep2$megakaryocyte_erythroid_lineage$plot,
    all.donor.subSample.list$rep3$hsc_lineage$plot,
          all.donor.subSample.list$rep3$lymphoid_lineage$plot,
          all.donor.subSample.list$rep3$granulocyte_monocyte_lineage$plot,
          all.donor.subSample.list$rep3$megakaryocyte_erythroid_lineage$plot,
    nrow=3, ncol = 4,
legend = "bottom")

ggsave(sub_samples,
       filename = "01_SubSamples.png",
       path = "Figures/SuppData",
       width = 16, height = 10)

# Save
nullList <- lapply(names(all.donor.subSample.list), function(rep_i){
    ob.list <- all.donor.subSample.list[[rep_i]]
    names <- names(ob.list)
    for (i in names){
        file_name <- paste0(dirPath, rep_i, "/", i,".RDS")
        saveRDS(ob.list[[i]]$obj, file_name)
        
    }
})
