# Title: Create Image for Main Article (B)
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Required Packages
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(ggrepel))

# Set Paths relative to project
dirPath <- "Tables/"
outDir <- "Figures/MainArticle/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load for Zero-Inflation
zi_table <- read.table(paste0(dirPath, "01_ZI_Performance.Table.tsv"), header = T)

# Load skew table
skew_table <- read.table(paste0(dirPath, "02_SkewSplit_Performance.Table.tsv"), header = T)

# Load Length table
len_table <- read.table(paste0(dirPath, "03_Length_Performance.Table.tsv"), header = T)

# Subset tables

# Zi-at 60
zi_table <- zi_table[zi_table$parameter.value == 60, , drop = FALSE]

# Skew-at 0.1
skew_table <- skew_table[skew_table$parameter.value == 0.1, , drop = FALSE]

# Length at 1000 -2000
len_table <- len_table[len_table$parameter.value == 1000.2000, , drop = FALSE]

# Combine frame
evaluation.frame <- rbind(zi_table, skew_table, len_table)

# ROC Curve (Figure-1B)
roc_B <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = parameter)) +
  geom_path(linewidth = 1, alpha = 0.7) +
  scale_color_manual(
    labels = c(
      "Different Lengths: 1000 & 2000",
      "Skew: nCells Start >> nCells", "Zero-Inflation of 60%"
    ),
    values = c(
      colorConesa(6)[1],
      colorConesa(6)[2],
      colorConesa(6)[6]
    ),
  ) +
  scale_x_continuous(
    breaks = c(0, 0.01, 0.10, 0.05),
    limits = c(0, 0.10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    limits = c(0, 1)
  ) +

  # Zero-Inflation
  geom_point(
    data = subset(evaluation.frame, RSQ == 0.7 & parameter == "ZI"),
    color = colorConesa(6)[6],
    size = 2.5
  ) +
  geom_label_repel(
    data = subset(evaluation.frame, RSQ == 0.7 & parameter == "ZI"),
    aes(
      label = RSQ,
      color = parameter
    ),
    nudge_x = -0.01, nudge_y = 0.05, # Adjust these values as needed
    size = 4,
    segment.color = colorConesa(6)[6], # Color of the line connecting label and point
    segment.size = 0.5, # Size of the line connecting label and point
    show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  ) +
    
    
  geom_point(
    data = subset(evaluation.frame, TPR > 0.8 & FPR <= 0.05 & parameter == "ZI"),
    color = colorConesa(6)[6],
    size = 2.5
  ) +
  geom_label_repel(
    data = subset(evaluation.frame,  TPR > 0.8 & FPR <= 0.05 & parameter == "ZI"),
    aes(
      label = RSQ,
      color = parameter
    ),
    nudge_x =  -0.01, nudge_y = 0.07, # Adjust these values as needed
    size = 4,
    segment.color = colorConesa(6)[6], # Color of the line connecting label and point
    segment.size = 0.5, # Size of the line connecting label and point
    show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  ) +

  # Skew
  geom_point(
    data = subset(evaluation.frame, RSQ == 0.7 & parameter == "SkewSplit"),
    color = colorConesa(3)[2],
    size = 2.5
  ) +
  geom_label_repel(
    data = subset(evaluation.frame, RSQ == 0.7 & parameter == "SkewSplit"),
    aes(
      label = RSQ,
      color = parameter
    ),
    nudge_x = -0.005, nudge_y = 0.05, # Adjust these values as needed
    size = 4,
    segment.color = colorConesa(3)[2], # Color of the line connecting label and point
    segment.size = 0.5, # Size of the line connecting label and point
    show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  ) +
  geom_point(
    data = subset(evaluation.frame, TPR >= 0.75 & FPR <= 0.05 & parameter == "SkewSplit"),
    color = colorConesa(3)[2],
    size = 2.5
  ) +
  geom_label_repel(
    data = subset(evaluation.frame, TPR >= 0.75 & FPR <= 0.05 & parameter == "SkewSplit"),
    aes(
      label = RSQ,
      color = parameter
    ),
    nudge_x = -0.005, nudge_y = 0.05, # Adjust these values as needed
    size = 4,
    segment.color = colorConesa(3)[2], # Color of the line connecting label and point
    segment.size = 0.5, # Size of the line connecting label and point
    show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  ) +

  # Length
  geom_point(
    data = subset(evaluation.frame, RSQ == 0.7 & parameter == "len"),
    color = colorConesa(3)[1],
    size = 2.5
  ) +
  geom_label_repel(
    data = subset(evaluation.frame, RSQ == 0.7 & parameter == "len"),
    aes(
      label = RSQ,
      color = parameter
    ),
    nudge_x = 0.01, nudge_y = -0.1, # Adjust these values as needed
    size = 4,
    segment.color = colorConesa(3)[1], # Color of the line connecting label and point
    segment.size = 0.5, # Size of the line connecting label and point
    show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  ) +
  geom_point(
    data = subset(evaluation.frame, TPR >= 0.75 & FPR < 0.05 & parameter == "len"),
    color = colorConesa(3)[1],
    size = 2.5
  ) +
  geom_label_repel(
    data = subset(evaluation.frame, TPR >= 0.75 & FPR < 0.05 & parameter == "len"),
    aes(
      label = RSQ,
      color = parameter
    ),
    nudge_x = 0.01, nudge_y = -0.07, # Adjust these values as needed
    size = 4,
    segment.color = colorConesa(3)[1], # Color of the line connecting label and point
    segment.size = 0.5, # Size of the line connecting label and point
    show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  ) +
  labs(
    title = "ROC Curve: Performance on Simulated Data",
    subtitle = "Varying R-Square",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Simulated Complexity"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
    # legend.position = "bottom"
    legend.position = c(0.35, 0.5), legend.justification = c("left", "top")
  ) +
  geom_vline(xintercept = 0.01, colour = "lightgrey", linetype = "dotted") + # Highlighted the x-intercept of 0.01
  geom_vline(xintercept = 0.05, colour = "lightgrey", linetype = "dotted") +
  geom_vline(xintercept = 0.1, colour = "lightgrey", linetype = "dotted") +
  geom_hline(yintercept = 0.8, colour = "lightgrey", linetype = "dotted") +
  guides(color = guide_legend(key_width = unit(5, "cm"), key_height = unit(4, "cm")))


# Save
saveRDS(object = roc_B,file = paste0(outDir, "MainArticle_FigureB.RDS"))
