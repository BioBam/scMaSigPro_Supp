# Title: Plot Binned Counts Counts of a Bifurcating Trajectory
# Author: Priyansh Srivastava
# Email: spriyansh29@gmail.com
# Year: 2023

# Required Packages
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorConesa))
suppressPackageStartupMessages(library(scMaSigPro))

# Set Paths relative to project
dirPath <- "Article_Image/data/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load the dataset
load(file = paste0(dirPath, "sce/ArticleFigure.RData"))

# Convert to scMaSigPro Object
scmp.sce <- as_scmp(object = sim.sce,
                    from = "sce",
                    additional_params = list(labels_exist = TRUE,
                                             existing_path_colname = "Group",
                                             existing_pseudotime_colname = "Step")
                    )

# Apply Simple binning
scmp.sce <- squeeze(
    scmpObject = scmp.sce,
    drop.fac = 1,
    verbose = T,
    split_bins = F, prune_bins = F, drop_trails = F
)

# Load the gene from figure-1A
fig1A <- readRDS("Article_Image/data/ImageRDS/Figure_1_A.RDS")

# get name of the gene plotted
fig1A_gene <- unique(fig1A$data$gene_name)

# Extract the binned counts
binned.counts.gene_i <- scmp.sce@compress.sce@assays@data@listData$bulk.counts[gene_i, ] %>%
    data.frame() %>%
    `colnames<-`("binned_counts") %>%
    rownames_to_column(var = "bin") 

# Get Simulated steps per groups
gene_i_bin_metadata <- scmp.sce@compress.sce@colData %>% 
    as.data.frame() %>%
    select(c("scmp_binned_pseudotime", "Path")) %>%
    rownames_to_column(var = "bin") 


# Left join by cell
gene_i_bin_plot_df <- left_join(binned.counts.gene_i, gene_i_bin_metadata, by = "bin")

# Plot
img_B <- ggplot(data = gene_i_bin_plot_df, aes(x = scmp_binned_pseudotime, y = binned_counts,
                                               group = Path, color = Path)) +
    geom_point(aes(color = Path), alpha = 0.7, stroke = 0.2) +
    geom_smooth(aes(color = Path), method = "loess", se = FALSE, linewidth = 1.5, alpha = 1) +
    theme_minimal(base_size = 15) +
    ggtitle("B. Pseudo-Bulked Expression along Continuum", subtitle = "Relative Expression across two Paths") +
    scale_color_manual(values = c("Path1" = RColorConesa::colorConesa(2, palette = "warm")[1], 
                                  "Path2" = RColorConesa::colorConesa(2, palette = "warm")[2])) +
    labs(x = "scMaSigPro Binned Pseudotime", y = "Pseudo-Bulked Expression Counts", color = "Path") +
    scale_x_continuous(breaks = seq(from = 0, to = max(gene_i_bin_plot_df$scmp_binned_pseudotime), by = 3)) +
    scale_y_continuous(breaks = seq(from = 0, to = max(gene_i_bin_plot_df$binned_counts), by = 500)) +
    theme(legend.position = "bottom",
          panel.grid.major = element_line(color = "lightgrey", size = 0.1),
          panel.grid.minor = element_line(color = "lightgrey", size = 0.1, linetype = "dotted"),
          axis.text.x = element_text(face = "bold", color = "black", size = 12), # Highlight x-axis ticks
          axis.text.y = element_text(face = "bold", color = "black", size = 12),  # Highlight y-axis ticks
          axis.ticks = element_line(color = "grey", size = 1), # Highlight axis ticks
          axis.ticks.length = unit(0.25, "cm")
    )
ggpubr::ggarrange(fig1A, img_B)

# Save Images and Save Image Object
dir.create("Article_Image/data/ImageRDS/", showWarnings = F)
saveRDS(img_B, "Article_Image/data/ImageRDS/Figure_1_B.RDS")

# Write Image
ggsave(plot = img_B,
       filename = "Article_Image/data/ImageRDS/Figure_1_B.png",
       dpi = 600)
