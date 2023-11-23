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
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))

# Load scMaSigPro
suppressPackageStartupMessages(library(scMaSigPro))

# Load object
object.list <- lapply(c("don1", "don2", "don3"), function(don){
    readRDS(paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/Monocle3_Processed_",don, ".RDS"))
})
names(object.list) <- c("don1", "don2", "don3")

# Donor-1
# Root-1: Y_36
# Path1: Y_4,Y_13,Y_20,Y_26,Y_35,Y_36,Y_37,Y_38,Y_43,Y_44,Y_47,Y_48,Y_52 (HMP-> pMono)
# Path2: Y_11,Y_25,Y_33,Y_36,Y_40,Y_42,Y_46 (HMP-> pMyel)

# Donor-2
# Root-1: Y_45
# Path1: Y_4,Y_6,Y_7,Y_8,Y_18,Y_26,Y_31,Y_36,Y_44,Y_45,Y_49,Y_51,Y_52,Y_53,Y_57,Y_58,Y_59,Y_69,Y_71,Y_72,Y_78,Y_81,Y_96(HMP-> pDend)
# Path2: Y_12,Y_15,Y_16,Y_17,Y_21,Y_24,Y_28,Y_29,Y_30,Y_32,Y_33,Y_40,Y_43,Y_45,Y_54,Y_56,Y_61,Y_66,Y_68,Y_74,Y_77,Y_79,Y_91,Y_93,Y_95,Y_98 (HMP-> pMono)

# Donor-3
# Root: Y_44
# Path-1:Y_26,Y_28,Y_31,Y_40,Y_41,Y_43,Y_44,Y_46,Y_63,Y_81,Y_82,Y_103,Y_114,Y_117,Y_119,Y_163 (Myeloid)
# Path-2:Y_8,Y_12,Y_33,Y_44,Y_45,Y_64,Y_65,Y_66,Y_104,Y_107,Y_116,Y_127,Y_130,Y_135,Y_145,Y_149,Y_150 (Dendritic)

# # Load object
# scMaSigPro.list <- lapply(object.list, function(don){
#     # Convert the ScMaSigPro Object
#     scmp.obj <- as_scmp(don, from = "cds",
#                         align_pseudotime = F,
#                         annotation_colname = "cell_type")
#     return(scmp.obj)
# })

scMaSigPro.list <- lapply(c("don1", "don2", "don3"), function(don){
    readRDS(paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/scMaSigPro_Processed_",don, ".RDS"))
})
names(scMaSigPro.list) <- c("don1", "don2", "don3")

# Run scMaSigPro
scMaSigPro.list <- lapply(scMaSigPro.list, function(don){


    # Compress
    scmp.obj <- squeeze(scmpObject = don,
                        split_bins = T,
                        prune_bins = F,
                        drop_trails = T,
                        drop.fac = 0.7
                        )
    # Make Design
    scmp.obj <- sc.make.design.matrix(scmp.obj,
                                      poly_degree = 3,
    )
    
    # Run p-vector
    scmp.obj <- sc.p.vector(
        parallel = T,
        scmpObj = scmp.obj, verbose = T,
        max_it = 10000,
        globalTheta = T,
        logOffset = F,
        useInverseWeights = F,
        logWeights = F,
        useWeights = T,
        offset = T,
        min.obs = 1)
    scmp.obj@distribution <- MASS::negative.binomial(10)
    
    # Run-Step-2
    scmp.obj <- sc.T.fit(
        scmpObj = scmp.obj, verbose = T,
        step.method = "backward"
    )
    
    
    return(scmp.obj)
})

sc.plot.bins.bar(scMaSigPro.list$don2) + sc.plot.bins.tile(scMaSigPro.list$don2)

# save
# Save the object for further analysis
scMaSigPro.list <- lapply(names(scMaSigPro.list), function(don){
    saveRDS(scMaSigPro.list[[don]],
        paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/scMaSigPro_Processed_",don, ".RDS"))
})

