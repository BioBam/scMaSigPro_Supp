# Load Library
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(maSigPro))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtools))
source("old_Scripts/classify_sig_genes.R")

# File-Location
gtPath <- "supp/07_Different_Length_of_Pseudotime/data/input/PseudobulkInput/"
inPath <- "supp/07_Different_Length_of_Pseudotime/data/output/MaSigPro/tstep/"
outPath <- "supp/07_Different_Length_of_Pseudotime/data/output/MaSigPro/ComparisonEvaluate/"

# Load File names
fileNames <- list.files(inPath)
names(fileNames) <- sapply(str_split(fileNames, "\\."), function(x) {
  return(c(x[1]))
})

# Some Variables
rsq <- 0.6

# Performance List
performance.list <- list()

## Run for Loop one-by-one
for (i in names(fileNames)) {
  # Load tstep
  tstepObj <- readRDS(paste0(inPath, fileNames[[i]]))
  # tstepObj <- readRDS(paste0(inPath, "pathA_1100_pathB_1900.RDS"))

  print(paste("Loaded for ", i))

  # Load ground truth
  gt <- readRDS(paste0(gtPath, fileNames[[i]]))
  groundTruth <- gt$gt

  print(paste("Extracted GT for ", i))

  # Checking For R-Square
  gene.change <- rownames(groundTruth[groundTruth$cobraTruth == 1, ])
  gene.no.change <- rownames(groundTruth[groundTruth$cobraTruth == 0, ])

  print(paste("Extracted Changes for for ", i))

  # Create Dataframe to plot Performance
  df.performance <- data.frame(
    VARIABLE = 0, TP = 0, FP = 0, TN = 0, FN = 0,
    TPR = 0, FPR = 0, FNR = 0, TNR = 0,
    INFLU = 0, INFLU_FP = 0, ACCURACY = 0
  )

  # Run For all Values of R square
  for (j in seq(0.05, 0.95, 0.05)) {
    print(paste("Checking for", j))
    # stop()

    # Classify Genes
    tryCatch(
      expr = {
        classifications <- classify_sig_genes(
          masigpro.tstep = tstepObj,
          rsq.value = j, draw.venn = F,
          all_genes = rownames(groundTruth),
          num_group = 2
        )
        # True Positive
        tp <- intersect(gene.change, classifications$detected.genes)

        # True Negative
        tn <- intersect(gene.no.change, classifications$undetected.gene)

        # False Positive
        fp <- intersect(gene.no.change, classifications$detected.genes)

        # False Negative
        fn <- intersect(gene.change, classifications$undetected.gene)

        # Influential Genes
        influ <- colnames(tstepObj$influ.info)

        # Influential + False Positive
        fp_influential <- intersect(fp, colnames(tstepObj$influ.info))

        # Sensitivity, True Positive Rate / Recall
        senstivity <- length(tp) / (length(tp) + length(fn))

        # False Negative Rate
        fnr <- length(fn) / (length(tp) + length(fn))

        # Specificity / True Negative Rate
        specificity <- length(tn) / (length(tn) + length(fp))

        # False Positive Rate
        fpr <- 1 - specificity

        # Accuracy
        accuracy <- (length(tp) + length(tn)) / (length(tp) + length(tn) + length(fp) + length(fn))

        # Make the vector of results
        res.perfom <- data.frame(
          VARIABLE = j,
          TP = length(tp),
          FP = length(fp),
          TN = length(tn),
          FN = length(fn),
          TPR = round(senstivity, 3),
          FPR = round(fpr, 3),
          FNR = round(fnr, 3),
          TNR = round(specificity, 3),
          INFLU = length(influ),
          INFLU_FP = length(fp_influential),
          ACCURACY = round(accuracy, 3)
        )

        print(paste("Completed for", j))

        # Add vector to the frame
        df.performance <- rbind(df.performance, res.perfom)
        gc()
      },
      error = function(e) {
        print("none found")
      }
    )
  }

  print(paste("Complete for ", i))

  # Remove the first Row
  df.performance <- df.performance[-1, ]

  # Add Info
  splits <- str_split(i, "_")[[1]]
  df.performance[[splits[1]]] <- splits[2]
  df.performance[[splits[3]]] <- splits[4]

  # Add to list
  performance.list[[i]] <- df.performance
}

# Get the Performace frame
performance.frame <- do.call(rbind, performance.list)
# Save the performace Table
dir.create(outPath, showWarnings = F, recursive = T)
saveRDS(performance.frame, paste0(outPath, "PerformanceTable.RDS"))

# Group
performance.frame$group <- paste0("pathA_", performance.frame$pathA, "_pathB_", performance.frame$pathB)
performance.frame <- performance.frame[performance.frame$group %in% c(
  # "pathA_100_pathB_2900", "pathA_200_pathB_2800", "pathA_300_pathB_2700"#,
  # "pathA_400_pathB_2600", "pathA_500_pathB_2500", "pathA_600_pathB_2400"#,
  # "pathA_700_pathB_2300", "pathA_800_pathB_2200", "pathA_900_pathB_2100"#,
  "pathA_1000_pathB_2000", "pathA_1100_pathB_1900", "pathA_1200_pathB_1800",
  "pathA_1300_pathB_1700", "pathA_1400_pathB_1600", "pathA_1500_pathB_1500"
), ]

# Comparitive ROC
acc <- ggplot(performance.frame, aes(
  x = VARIABLE, y = ACCURACY,
  color = group
)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.05)) +
  geom_path(linewidth = 1.5, alpha = 0.6) +
  scale_y_continuous(breaks = seq(0.7, 1, 0.05)) +
  # scale_color_brewer(palette = "Set1") +
  labs(
    title = "Accuracy Against Changing R Square",
    subtitle = "Red Dots: False Negatives",
    x = "Increasing Value for R-Square",
    y = "Accuracy"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
    axis.text = element_text(size = rel(1.5)),
    legend.key.size = unit(4, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.position = "bottom",
    legend.key.height = unit(1, "cm"),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
  )

# ROC
roc <- ggplot(performance.frame, aes(x = FPR, y = TPR, color = group)) +
  geom_point() +
  geom_path(linewidth = 1.5, alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_y_continuous(breaks = seq(0.5, 1, 0.05)) +
  # scale_color_brewer(palette = "Set1") +
  labs(
    title = "ROC-curve, Different Values of R-Square",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
    legend.position = "bottom",
    legend.key.size = unit(4, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = rel(1.5)),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
  ) +
  geom_vline(xintercept = 0.01, colour = "lightgrey") +
  geom_vline(xintercept = 0.05, colour = "grey") +
  geom_vline(xintercept = 0.1, colour = "darkgrey") +
  guides(color = guide_legend(
    key_width = unit(1, "cm"),
    key_height = unit(1, "cm")
  ))


ggsave(acc,
  filename = paste0(outPath, "acc.png"),
  dpi = 600, height = 10, width = 12
)
ggsave(roc,
  filename = paste0(outPath, "roc.png"),
  dpi = 600, height = 10, width = 12
)
