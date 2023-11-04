# Title: Plot Raw Counts of the Bifurcating Trajectory
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Required Packages
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorConesa))

# Set Paths relative to project
dirPath <- "Article_Image/data/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load the dataset
load(file = paste0(dirPath, "sce/ArticleFigure.RData"))

# Extract random gene that is differentially expressed between two paths
gene.metadata <- rowData(sim.sce) %>% as.data.frame()

# Get differentially expressed genes
gene.metadata.diff <- gene.metadata[gene.metadata$BaseGeneMean > mean(gene.metadata$BaseGeneMean) 
                                    & gene.metadata$status2 == "High_FC" & 
                                        gene.metadata$status == "Opposite_Change",, drop = F]

# Get any gene-i
set.seed(1)
gene_i <- gene.metadata.diff[sample(c(1:nrow(gene.metadata.diff)), size = 1, ), "gene_short_name"]

# Plot Gene_i's raw counts against simulated pseudotime/steps
gene_i_counts <- sim.sce@assays@data@listData$counts[gene_i, ] %>%
    data.frame() %>%
    `colnames<-`("raw_counts") %>%
    rownames_to_column(var = "cell") 
gene_i_counts$gene_name <- gene_i

# Get Simulated steps per groups
gene_i_metadata <- sim.sce@colData %>% 
    as.data.frame() %>%
    select(c(Group, Step)) %>%
    rownames_to_column(var = "cell") 

# Left join by cell
gene_i_plot_df <- left_join(gene_i_counts, gene_i_metadata, by = "cell")

# Plot
img_A <- ggplot(data = gene_i_plot_df, aes(x = Step, y = raw_counts, group = Group, color = Group)) +
    geom_point(aes(color = Group), alpha = 0.2, stroke = 0.2) +
    geom_smooth(aes(color = Group), method = "loess", se = FALSE, linewidth = 1.5, alpha = 1) +
    theme_minimal(base_size = 15) +
    ggtitle("A. Raw Expression along Continuum", subtitle = "Relative Expression across two Paths") +
    scale_color_manual(values = c("Path1" = RColorConesa::colorConesa(2, palette = "warm")[1], 
                                  "Path2" = RColorConesa::colorConesa(2, palette = "warm")[2])) +
    labs(x = "Pseudotime like Continuum", y = "Gene Expression Counts", color = "Path") +
    scale_x_continuous(breaks = seq(from = 0, to = max(gene_i_plot_df$Step), by = 250)) +
    scale_y_continuous(breaks = seq(from = 0, to = max(gene_i_plot_df$raw_counts), by = 10)) +
    theme(legend.position = "bottom",
          panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1),
          panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.1, linetype = "dotted"),
          axis.text.x = element_text(face = "bold", color = "black", size = 12), # Highlight x-axis ticks
          axis.text.y = element_text(face = "bold", color = "black", size = 12),  # Highlight y-axis ticks
          axis.ticks = element_line(color = "grey", size = 1), # Highlight axis ticks
          axis.ticks.length = unit(0.25, "cm")
    )

# Save Images and Save Image Object
dir.create("Article_Image/data/ImageRDS/", showWarnings = F)
saveRDS(img_A, "Article_Image/data/ImageRDS/Figure_1_A.RDS")

# Write Image
ggsave(plot = img_A,
       filename = "Article_Image/data/ImageRDS/Figure_1_A.png",
       dpi = 600)
