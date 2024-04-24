# Title: Analyze 60% Zero-Inflated Data with TradeSeq
# Author: Priyansh Srivastava
# Year: 2023

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(pryr))

# Set paths
dirPath <- "/supp_data/ComparisonWithTradeSeq/simulated/sce/"
resPath <- "/supp_data/ComparisonWithTradeSeq/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))

# ReadData
load(paste0(dirPath, "time_20k_cells.RData"))

# Running scMaSigPro
scmp.obj <- as_scmp(sim.sce,
                    from = "sce",
                    align_pseudotime = T,
                    additional_params = list(
                        labels_exist = TRUE,
                        exist_ptime_col = "Step",
                        exist_path_col = "Group"
                    ), verbose = F
)

# Squeeze
scmp.obj <- sc.squeeze(
    scmpObj = scmp.obj,
    bin_method = "Sturges",
    drop_fac = 200,
    verbose = F,
    aggregate = "sum",
    split_bins = F,
    prune_bins = F,
    drop_trails = F,
    fill_gaps = F
)

# Make Design
scmp.obj <- sc.set.poly(scmp.obj,
                        poly_degree = 3
)

# Benchmark time
mbm <- microbenchmark(
    "ScMaSigPro_22_CPU__sc.p.vector()" = {
        print("Computing scMaSigPro::sc.p.vector()")
        # Run p-vector
        scmp.obj <- sc.p.vector(
            scmpObj = scmp.obj, verbose = T,
            min_na = 1,
            parallel = T,
            offset = T,
            max_it = 1000
        )
        gc()
        print("Finished scMaSigPro::sc.p.vector()")
        
    },
    "ScMaSigPro_22_CPU_sc.t.fit()" = {
        print("Computing scMaSigPro::sc.t.fit()")
        # Run-Step-2
        scmp.obj <- sc.t.fit(
            scmpObj = scmp.obj, verbose = F,
            selection_method = "backward", parallel = T,
            offset = T
        )
        gc()
        print("Finished scMaSigPro::sc.t.fit()")
    },
    times = 5
)

# Process the results
data <- summary(mbm) %>% as.data.frame()
data$min_mean <- paste(round(data$mean / 60, digits = 1), "minutes")

compareBar_Time <- ggplot(data, aes(x = expr, y = mean, fill = expr)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(
        breaks = seq(0, 120, 20),
        limits = c(0, 120)
    ) +
    labs(
        title = "Execution Times for a bifurcating trajectory",
        subtitle = "Number of Cells: 1500; Number of Genes: 1000",
        x = "Method",
        y = "Time (seconds)"
    ) +
    geom_text(aes(label = min_mean),
              position = position_dodge(width = 0.9),
              size = 3,
              vjust = 0.5, hjust = -0.1
    ) +
    coord_flip() +
    scale_fill_viridis(
        discrete = TRUE, name = "Custom Legend Title",
        breaks = c("TradeSeq_1_CPU", "ScMaSigPro_1_CPU", "TradeSeq_8_CPU", "ScMaSigPro_8_CPU"),
        labels = c("Custom Label 1", "Custom Label 2", "Custom Label 3", "Custom Label 4")
    ) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "none", legend.justification = "left", legend.box.just = "left")

compareBar_Time

# Save
ggsave(
    plot = compareBar_Time,
    filename = paste0("/supp_data/Figures/SuppData/04_tradeSeq_Time.png"),
    dpi = 300, width = 10
)
