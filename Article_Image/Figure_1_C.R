# Title: Plot Evaluation Metrics
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggrepel))

# Load Evaluation
evaluation.frame.zi <- as.data.frame(read.table("benchmarks/01_Sparsity/data/output/Performance.Table.tsv", header = T))
evaluation.frame.skew <- as.data.frame(read.table("benchmarks/02_Skewness/data/output/Performance.Table.tsv", header = T))
evaluation.frame.len <- as.data.frame(read.table("benchmarks/03_Different_Length/data/output/Performance.Table.tsv", header = T))

# Subset
evaluation.frame.zi.sub <- evaluation.frame.zi[(
    evaluation.frame.zi$parameter == "ZI" & evaluation.frame.zi$parameter.value == 60
    ),]
evaluation.frame.skew.sub <- evaluation.frame.skew[(
    evaluation.frame.skew$parameter == "Skew" & evaluation.frame.skew$parameter.value == 0.9
),]
evaluation.frame.len.sub <- evaluation.frame.len[(
    evaluation.frame.len$parameter == "Length" & evaluation.frame.len$parameter.value == "400_and_2600"
),]

# Combine
evaluation.frame <- rbind(evaluation.frame.zi.sub,
                          evaluation.frame.skew.sub,
                          evaluation.frame.len.sub)

# Assuming evaluation.frame is your data frame
# Extract starting points
start_points <- subset(evaluation.frame, RSQ == 0.7 & parameter == "ZI")

# Define end coordinates (example)
end_x <- 0.8
end_y <- 0.02

# Create a new data frame for segments
segments_data <- data.frame(xstart = start_points$TPR, 
                            ystart = start_points$FPR, 
                            xend = 0.8, 
                            yend = 0.05)

# ROC Curve
roc <- ggplot(evaluation.frame, aes(x = FPR, y = TPR, color = parameter)) +
    geom_path(linewidth = 1, alpha = 0.7) +
<<<<<<< HEAD
    scale_color_manual( labels = c("Different Lengths: 400 & 2600",
                                   "Skewness: nCells Start >> nCells End", 
=======
    scale_color_manual(labels = c("Different Lengths: 400 & 2600",
                                   "Skew: nCells Start >> nCells", 
>>>>>>> dev
                                   "Zero-Inflation of 60%"),
                        
                        values = c(colorConesa(6)[1],
                                  colorConesa(6)[2],
                                  colorConesa(6)[6]),
                       ) +
    scale_x_continuous(breaks = seq(0, 0.10, 0.05), 
                       limits = c(0, 0.10)
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0, 1)
    ) +
    
    # Zero-Inflation
    geom_point(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "ZI"),
        color = colorConesa(6)[6],
<<<<<<< HEAD
        size = 2) +
=======
        size = 2.5) +
>>>>>>> dev
    
    geom_label_repel(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "ZI"),
        aes(label = RSQ,
            color = parameter),
        nudge_x = 0.007, nudge_y = -0.12, # Adjust these values as needed
<<<<<<< HEAD
            size = 3,
=======
            size = 4,
>>>>>>> dev
            segment.color = colorConesa(6)[6], # Color of the line connecting label and point
            segment.size = 0.5, # Size of the line connecting label and point
        show.legend = FALSE # Prevent the creation of a legend for this aesthetic
    )+
   
    geom_point(
        data = subset(evaluation.frame, TPR >= 0.8 & FPR < 0.01 & parameter == "ZI"),
        color = colorConesa(6)[6],
<<<<<<< HEAD
        size = 2) +
=======
        size = 2.5) +
>>>>>>> dev

    
    # Skew
    geom_point(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Skew"),
        color = colorConesa(3)[2],
<<<<<<< HEAD
        size = 2) +
=======
        size = 2.5) +
>>>>>>> dev
    
    geom_label_repel(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Skew"),
        aes(label = RSQ,
            color = parameter),
<<<<<<< HEAD
        nudge_x = 0.04, nudge_y = -0.1, # Adjust these values as needed
        size = 3,
=======
        nudge_x = 0.01, nudge_y = -0.1, # Adjust these values as needed
        size = 4,
>>>>>>> dev
        segment.color = colorConesa(3)[2], # Color of the line connecting label and point
        segment.size = 0.5, # Size of the line connecting label and point
        show.legend = FALSE # Prevent the creation of a legend for this aesthetic
    )+
    
    geom_point(
        data = subset(evaluation.frame, TPR >= 0.8 & FPR <= 0.01 & parameter == "Skew"),
        color = colorConesa(3)[2],
<<<<<<< HEAD
        size = 2) +
=======
        size = 2.5) +
>>>>>>> dev
    geom_label_repel(
        data = subset(evaluation.frame, TPR >= 0.8 & FPR <= 0.01 & parameter == "Skew"),
        aes(label = RSQ,
            color = parameter),
        nudge_x = 0.04, nudge_y = -0.2, # Adjust these values as needed
<<<<<<< HEAD
        size = 3,
=======
        size = 4,
>>>>>>> dev
        segment.color = colorConesa(3)[2], # Color of the line connecting label and point
        segment.size = 0.5, # Size of the line connecting label and point
        show.legend = FALSE # Prevent the creation of a legend for this aesthetic
    )+
   
    # Length
    geom_point(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Length"),
        color = colorConesa(3)[1],
<<<<<<< HEAD
        size = 2) +
=======
        size = 2.5) +
>>>>>>> dev
    geom_label_repel(
        data = subset(evaluation.frame, RSQ == 0.7 & parameter == "Length"),
        aes(label = RSQ,
            color = parameter),
        nudge_x = 0.0, nudge_y = 0.1, # Adjust these values as needed
<<<<<<< HEAD
        size = 3,
=======
        size = 4,
>>>>>>> dev
        segment.color = colorConesa(3)[1], # Color of the line connecting label and point
        segment.size = 0.5, # Size of the line connecting label and point
        show.legend = FALSE # Prevent the creation of a legend for this aesthetic
    )+
   
    geom_point(
        data = subset(evaluation.frame, TPR >= 0.9 & FPR < 0.01 & parameter == "Length"),
        color = colorConesa(3)[1],
<<<<<<< HEAD
        size = 2) +
=======
        size = 2.5) +
>>>>>>> dev
    geom_label_repel(
        data = subset(evaluation.frame, TPR >= 0.9 & FPR < 0.01 & parameter == "Length"),
        aes(label = RSQ,
            color = parameter),
        nudge_x = 0.01, nudge_y = 0.07, # Adjust these values as needed
<<<<<<< HEAD
        size = 3,
=======
        size = 4,
>>>>>>> dev
        segment.color = colorConesa(3)[1], # Color of the line connecting label and point
        segment.size = 0.5, # Size of the line connecting label and point
        show.legend = FALSE # Prevent the creation of a legend for this aesthetic
    )+
  
    labs(
<<<<<<< HEAD
        title = "ROC: Performance on Simulated Data",
=======
        title = "ROC Curve: Performance on Simulated Data",
>>>>>>> dev
        subtitle = "Varying R-Square",
        x = "False Positive Rate (1-Specificity)",
        y = "True Positive Rate (Sensitivity)",
        color = "Simulation") +
<<<<<<< HEAD
    theme_classic(base_size = 10) +
    theme(legend.box = "vertical",
          legend.direction = "vertical",
        panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
        #legend.position = "bottom"
        legend.position = c(0.5, 0.5), legend.justification = c("left", "top")
=======
    theme_classic(base_size = 15) +
    theme(legend.box = "vertical",
          legend.direction = "vertical",
          legend.title=element_text(size=16), 
          legend.text=element_text(size=14),
        panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"),
        #legend.position = "bottom"
        legend.position = c(0.35, 0.5), legend.justification = c("left", "top")
>>>>>>> dev
    ) +
    geom_vline(xintercept = 0.01, colour = "lightgrey", linetype = "dotted") +  # Highlighted the x-intercept of 0.01
    geom_vline(xintercept = 0.05, colour = "lightgrey", linetype = "dotted") +
    geom_vline(xintercept = 0.1, colour = "lightgrey", linetype = "dotted") +
    geom_hline(yintercept = 0.8, colour = "lightgrey", linetype = "dotted") +
    guides(color = guide_legend(key_width = unit(5, "cm"), key_height = unit(4, "cm"))) +
    geom_label_repel(
        data =subset(evaluation.frame, TPR >= 0.8 & FPR < 0.01 & parameter == "ZI"),
        aes(label = RSQ,
            color = parameter),
        nudge_x = 0.02, nudge_y = -0.12, # Adjust these values as needed
        size = 3,
        segment.color = colorConesa(6)[6], # Color of the line connecting label and point
        segment.size = 0.5, # Size of the line connecting label and point
        show.legend = FALSE # Prevent the creation of a legend for this aesthetic
    )

print(roc)

ggsave(plot = roc,
       path = "Article_Image",
       dpi = 1000,  filename = "Figure1_A.png",
       width = 5, height = 5)
saveRDS(roc, file = "Article_Image/Figure1_A.RDS")

stop()

# Plot all values against zero inflation
long_data <- melt(evaluation.frame, id.vars = c("RSQ"), measure.vars = c("Accuracy", "Precision", "FPR", "TPR", "Recall", "F1_Score"))

# Plot performance
performance <- ggplot(long_data, aes(x = RSQ, y = value, group = interaction(variable), color = variable)) +
    geom_line(linewidth = 0.6) + 
    geom_point(size = 0.8) +
    scale_color_manual(values = colorConesa(6)) +
    labs(x = "Varying R-Square", y = "Performance Metric",
         title = "Performance Metric for different levels of zero-inflation",
         color = "Measure") +
    scale_x_continuous(breaks = seq(0.1, 0.95, 0.2), limits = c(0, 0.95)) +
    theme_minimal(base_size = 10) + 
    theme(legend.position = "bottom")

print(performance)

# Save All RDS
saveRDS(performance, file = paste(dirPath, "Performance_plot.RDS"))
saveRDS(acc, file = paste(dirPath, "Accuracy_Plot.RDS"))
saveRDS(roc, file = paste(dirPath, "ROC_plot.RDS"))

# Save Images
ggsave(acc,
       filename = paste0(dirPath, "Accuracy.png"),
       dpi = 600, height = 6, width = 6
)
ggsave(roc,
       filename = paste0(dirPath, "ROC.png"),
       dpi = 600, height = 6, width = 6
)

ggsave(performance,
       filename = paste0(dirPath, "Performance.png"),
       dpi = 600, height = 6, width = 6
)
