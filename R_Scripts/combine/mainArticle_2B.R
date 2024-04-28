# Title: Create Image for Main Article (B)
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Required Packages
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(ggrepel))

# Set paths
base_string <- "../scMaSigPro_supp_data/"
base_string_2 <- ""
benchmark_string_base_string <- paste0(base_string, "benchmarks/")
imgPath <- paste0(benchmark_string_base_string, "img/")
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Create Directory if does not exist
dir.create(figPath, showWarnings = FALSE, recursive = TRUE)
dir.create(imgPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_hd, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_lr, showWarnings = FALSE, recursive = TRUE)
dir.create(tabPath, showWarnings = FALSE, recursive = TRUE)
dir.create(benchmark_string_base_string, showWarnings = FALSE, recursive = TRUE)

# Load for Zero-Inflation
zi_table <- read.table(
  paste0(tabPath, "01_ZI_Performance.Table.tsv"),
  header = T, sep = "\t", stringsAsFactors = TRUE
)

# Load skew table with splot parameter
skew_split_table <- read.table(
  paste0(tabPath, "02_SkewSplit_Performance.Table.tsv"),
  header = T, sep = "\t", stringsAsFactors = TRUE
)

# Load Skew only
skew_table <- read.table(
  paste0(tabPath, "02_Skew_Performance.Table.tsv"),
  header = T, sep = "\t", stringsAsFactors = TRUE
)

# Load Length table
len_table <- read.table(
  paste0(tabPath, "03_Length_Performance.Table.tsv"),
  header = T, sep = "\t", stringsAsFactors = TRUE
)

# Zi-at 60
zi_table <- zi_table[zi_table$parameter.value == 60, , drop = FALSE]

# Skew-at 0.1
skew_table <- skew_table[skew_table$parameter.value == 0.1, , drop = FALSE]

# Skew-Split at 0.1
skew_split_table <- skew_split_table[skew_split_table$parameter.value == 0.1, , drop = FALSE]

# Length at 1000 -2000
len_table <- len_table[len_table$parameter.value == 1000.2000, , drop = FALSE]

# Combine frame
evaluation.frame <- rbind(
  len_table,
  # skew_table,
  skew_split_table,
  zi_table
)

# ROC Curve (Figure-1B)
roc_B <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = parameter)) +
  geom_path(linewidth = 1, alpha = 0.7) +
  scale_color_manual(
    labels = c(
      "Different Lengths: 1000 & 2000",
      # "Skew: nCells Start >> nCells",
      "Skew: nCells Start >> nCells (BinSplit)",
      "Zero-Inflation of 60%"
    ),
    values = c(
      colorConesa(6)[6],
      # colorConesa(6)[4],
      colorConesa(6)[2],
      colorConesa(6)[1]
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
    data = subset(evaluation.frame, RSQ == 0.6 & parameter == "ZI"),
    color = colorConesa(6)[1],
    shape = 17,
    size = 2.5
  ) +
  # geom_label_repel(
  #   data = subset(evaluation.frame, RSQ == 0.6 & parameter == "ZI"),
  #   aes(
  #     label = RSQ,
  #     color = parameter
  #   ),
  #   nudge_x = 0, nudge_y = 0.15, # Adjust these values as needed
  #   size = 4,
  #   segment.color = colorConesa(6)[1], # Color of the line connecting label and point
  #   segment.size = 0.5, # Size of the line connecting label and point
  #   show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  # ) +
  # geom_point(
  #   data = subset(evaluation.frame, TPR > 0.8 & TPR <= 0.85 & FPR <= 0.05 & parameter == "ZI"),
  #   color = colorConesa(6)[6],
  #   size = 2.5
  # ) +
  # geom_label_repel(
  #   data = subset(evaluation.frame, TPR > 0.8 & TPR <= 0.85 & FPR <= 0.05 & parameter == "ZI"),
  #   aes(
  #     label = RSQ,
  #     color = parameter
  #   ),
  #   nudge_x = -0.01, nudge_y = 0.07, # Adjust these values as needed
  #   size = 4,
  #   segment.color = colorConesa(6)[6], # Color of the line connecting label and point
  #   segment.size = 0.5, # Size of the line connecting label and point
  #   show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  # ) +

  # Skew split
  geom_point(
    data = subset(evaluation.frame, RSQ == 0.6 & parameter == "SkewSplit"),
    color = colorConesa(6)[2],
    shape = 17,
    size = 2.5
  ) +
  # geom_label_repel(
  #   data = subset(evaluation.frame, RSQ == 0.6 & parameter == "SkewSplit"),
  #   aes(
  #     label = RSQ,
  #     color = parameter
  #   ),
  #   nudge_x = 0, nudge_y = -0.2, # Adjust these values as needed
  #   size = 4,
  #   segment.color = colorConesa(6)[2], # Color of the line connecting label and point
  #   segment.size = 0.5, # Size of the line connecting label and point
  #   show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  # ) +
  # Skew
  # geom_point(
  #   data = subset(evaluation.frame, RSQ == 0.6 & parameter == "Skew"),
  #   color = colorConesa(6)[4],
  #   shape = 17,
  #   size = 2.5
  # ) +
  #     geom_label_repel(
  #         data = subset(evaluation.frame, RSQ == 0.6 & parameter == "Skew"),
  #         aes(
  #             label = RSQ,
  #             color = parameter
  #         ),
  #         nudge_x = 0.03, nudge_y = -0.08, # Adjust these values as needed
  #         size = 4,
  #         segment.color = colorConesa(6)[4], # Color of the line connecting label and point
  #         segment.size = 0.5, # Size of the line connecting label and point
  #         show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  #     ) +
  # Length
  geom_point(
    data = subset(evaluation.frame, RSQ == 0.6 & parameter == "len"),
    color = colorConesa(6)[6],
    shape = 17,
    size = 2.5
  ) +
  # geom_label_repel(
  #   data = subset(evaluation.frame, RSQ == 0.6 & parameter == "len"),
  #   aes(
  #     label = RSQ,
  #     color = parameter
  #   ),
  #   nudge_x = 0.01, nudge_y = -0.1, # Adjust these values as needed
  #   size = 4,
  #   segment.color = colorConesa(6)[6], # Color of the line connecting label and point
  #   segment.size = 0.5, # Size of the line connecting label and point
  #   show.legend = FALSE # Prevent the creation of a legend for this aesthetic
  # ) +
  labs(
    title = "ROC Curve: Performance on Simulated Data",
    subtitle = "Varying R-Square (0.1 to 0.95)",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Simulated Complexity"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
    # legend.position = "bottom"
    legend.position = c(0.2, 0.5), legend.justification = c("left", "top")
  ) +
  geom_vline(xintercept = 0.01, colour = "lightgrey", linetype = "dotted") + # Highlighted the x-intercept of 0.01
  geom_vline(xintercept = 0.05, colour = "lightgrey", linetype = "dotted") +
  geom_vline(xintercept = 0.1, colour = "lightgrey", linetype = "dotted") +
  geom_hline(yintercept = 0.8, colour = "lightgrey", linetype = "dotted") +
  guides(color = guide_legend(key_width = unit(5, "cm"), key_height = unit(4, "cm")))


roc_B

# Save
saveRDS(object = roc_B, file = paste0(imgPath, "MainFigure_2B.rds"))
