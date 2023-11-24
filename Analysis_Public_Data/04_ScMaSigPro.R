##################################
## Author: Priyansh Srivastava ###
## Email: spriyansh29@gmail.com ##
## Script: ScMaSigPro ############
##################################

set.seed(007)

# Prefix
dirPath <- "/supp_data/Analysis_Public_Data/"

# Get folder names
rep_vec <- list.dirs(dirPath, full.names = F, recursive = F)
rep_vec <- rep_vec[rep_vec != "Azimuth_Human_BoneMarrow"]
names(rep_vec) <- rep_vec
    
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
object.list <- lapply(rep_vec, function(rep_i, inPath = dirPath, outPath = dirPath) {
    sob <- readRDS(file = paste0(inPath, rep_i, "/", rep_i, "_cds.RDS"))
})

# Donor-1
# Root-1: Y_63
# Path1: Y_4,Y_18,Y_19,Y_24,Y_27,Y_29,Y_32,Y_37,Y_38,Y_47,Y_49,Y_50,Y_55,Y_56,Y_63,Y_67,Y_69,Y_74,Y_89,Y_100 (HSC-> Eryth)
# Path2: Y_3,Y_5,Y_8,Y_11,Y_16,Y_20,Y_25,Y_36,Y_43,Y_44,Y_51,Y_52,Y_63,Y_85,Y_90,Y_91,Y_94,Y_97,Y_98 (HSC-> Myel/Mono)

# Donor-2
# Root-1: Y_45
# Path1: Y_4,Y_6,Y_7,Y_8,Y_18,Y_26,Y_31,Y_36,Y_44,Y_45,Y_49,Y_51,Y_52,Y_53,Y_57,Y_58,Y_59,Y_69,Y_71,Y_72,Y_78,Y_81,Y_96(HMP-> pDend)
# Path2: Y_12,Y_15,Y_16,Y_17,Y_21,Y_24,Y_28,Y_29,Y_30,Y_32,Y_33,Y_40,Y_43,Y_45,Y_54,Y_56,Y_61,Y_66,Y_68,Y_74,Y_77,Y_79,Y_91,Y_93,Y_95,Y_98 (HMP-> pMono)

# Donor-3
# Root: Y_94
# Path-1:Y_5,Y_9,Y_23,Y_27,Y_28,Y_33,Y_38,Y_44,Y_50,Y_54,Y_55,Y_58,Y_66,Y_73,Y_76,Y_82,Y_86,Y_88,Y_89,Y_91,Y_94,Y_97,Y_99 (HSC-> Mono myelo)
# Path-2:Y_4,Y_10,Y_11,Y_13,Y_21,Y_22,Y_36,Y_40,Y_41,Y_49,Y_57,Y_67,Y_77,Y_81,Y_84,Y_92,Y_94,Y_95,Y_98 (HSC -> mega)

# Create ScMaSigPro
scMaSigPro.list <- lapply(object.list[1], function(don){
    # Convert the ScMaSigPro Object
    scmp.obj <- as_scmp(don, from = "cds",
                        align_pseudotime = F,
                        annotation_colname = "cell_type")
    return(scmp.obj)
})

# Run scMaSigPro
scMaSigPro.list <- lapply(scMaSigPro.list, function(don){
    
    
    # Compress
    scmp.obj <- squeeze(scmpObject = don,
                        split_bins = F,
                        prune_bins = F,
                        drop_trails = F,
                        drop_fac = 1
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
})

scmp.obj <- sc.get.siggenes(scmp.obj,
                rsq = 0.7,
                vars = "groups")

# save
# Save the object for further analysis
scMaSigPro.list <- lapply(names(scMaSigPro.list), function(don){
    saveRDS(scMaSigPro.list[[don]],
            paste0("Analysis_Public_Data/data/SingleCellExperimentAtlas/Monocle3_Input/scMaSigPro_Processed_",don, ".RDS"))
})
