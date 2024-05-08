# Title: Evaluate Individual Functions
# Author: Priyansh Srivastava
# Year: 2024

trash <- gc()
trash <- NULL

# Create Pid
pid_pvector_1cpu <- 101
pid_tfit_1cpu <- 101
pid_pvector_8cpu <- 101
pid_tfit_8cpu <- 101

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(viridis))

# Set paths
base_string <- "../../../scMaSigPro_supp_data/"
base_string_2 <- "../../"
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
tabPath <- paste0(base_string, "tables/")
helpScriptsDir <- paste0(base_string_2, "R_Scripts/helper_function/")

# Load Base data
paramEstimates <- readRDS(paste0(base_string, "benchmarks/00_Parameter_Estimation/output/setty_et_al_d1_splatEstimates.RDS"))

# Create Directory if does not exist
dir.create(figPath, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_hd, showWarnings = FALSE, recursive = TRUE)
dir.create(figPath_lr, showWarnings = FALSE, recursive = TRUE)
dir.create(tabPath, showWarnings = FALSE, recursive = TRUE)

# Load custom function
source(paste0(helpScriptsDir, "calcNormCounts.R"))
source(paste0(helpScriptsDir, "add_gene_anno().R"))
source(paste0(helpScriptsDir, "calc_bin_size.R"))

# Defines sets
evaluate_cells <- seq(8000, 20000, 4000)
drop_factor <- c(15, 35, 45, 55)
names(drop_factor) <- evaluate_cells
names(evaluate_cells) <- evaluate_cells

cat(paste("\nInitiate Simulation:", timestamp(suffix = "", prefix = "", quiet = TRUE)))

# Create list of simulated parameters
parameter_list <- mclapply(names(evaluate_cells), function(x) {
  x <- as.numeric(x)

  # Create Base parameters/ Same for All groups
  params.groups <- newSplatParams(
    batch.rmEffect = TRUE, # No Batch affect
    batchCells = x, # Number of Cells
    nGenes = 5000, # Number of Genes
    seed = 2022, # Set seed
    mean.rate = paramEstimates@mean.rate,
    mean.shape = paramEstimates@mean.shape,
    lib.scale = paramEstimates@lib.scale,
    lib.loc = paramEstimates@lib.loc,
    bcv.common = paramEstimates@bcv.common,
    bcv.df = paramEstimates@bcv.df,
    dropout.type = "experiment",
    group.prob = c(0.6, 0.4),
    path.from = c(0, 0),
    de.prob = 0.3,
    de.facLoc = 1,
    path.nonlinearProb = 0.3,
    path.sigmaFac = 0.5,
    out.facLoc = paramEstimates@out.facLoc,
    dropout.mid = paramEstimates@dropout.mid,
    out.facScale = paramEstimates@out.facScale,
    out.prob = paramEstimates@out.prob,
    path.skew = c(0.4, 0.6),
    dropout.shape = -0.5,
    path.nSteps = c(x / 2, x / 2)
  )

  return(params.groups)
}, mc.cores = availableCores())
names(parameter_list) <- names(evaluate_cells)

cat(paste("\nSimulated_Parameters:", timestamp(suffix = "", prefix = "", quiet = TRUE)))

# Create Empty daraframe
master_data <- data.frame()

# Loop over the parameters
for (i in names(parameter_list)) {
  # Get parameter
  params.groups <- parameter_list[[i]]
  d <- drop_factor[i]

  # Simulate Object
  sim.sce <- splatSimulate(
    params = params.groups,
    method = "paths",
    verbose = F
  )

  # Add gene Info
  gene.info <- add_gene_anno(sim.sce = sim.sce)
  gene.info <- gene.info[mixedsort(gene.info$gene_short_name), ]
  rowData(sim.sce) <- DataFrame(gene.info)

  cat(paste("\nSimulation completed for ", i, "cells :", timestamp(suffix = "", prefix = "", quiet = TRUE)))

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
    drop_fac = d,
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

  # initial set-up
  scmp.obj <- sc.p.vector(scmpObj = scmp.obj, verbose = FALSE, min_na = 1, parallel = T, n_cores = 8, offset = T, max_it = 1000)

  # Cleanup
  sim.sce <- NULL
  counts <- NULL
  rm(counts)
  rm(sim.sce)
  trash <- gc()
  trash <- NULL

  cat(paste("\nClean-up for scMaSigPro completed for", i, "cells :", timestamp(suffix = "", prefix = "", quiet = TRUE)))

  mbm <- microbenchmark(
    "ScMaSigPro_1_CPU_sc.p.vector" = {
      trash <- gc()
      trash <- NULL
      scmp.obj <- sc.p.vector(scmpObj = scmp.obj, verbose = FALSE, min_na = 1, parallel = F, offset = T, max_it = 1000)
      cat(paste("\nFinished ScMaSigPro_1_CPU_sc.p.vector(", pid_pvector_1cpu, "for", i, "cells", "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
      pid_pvector_1cpu <- pid_pvector_1cpu + 1
      trash <- gc()
      trash <- NULL
    },
    "ScMaSigPro_8_CPU_sc.p.vector" = {
      trash <- gc()
      trash <- NULL
      scmp.obj <- sc.p.vector(scmpObj = scmp.obj, verbose = FALSE, min_na = 1, parallel = T, n_cores = 8, offset = T, max_it = 1000)
      cat(paste("\nFinished ScMaSigPro_8_CPU_sc.p.vector(", pid_pvector_8cpu, "for", i, "cells", "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
      pid_pvector_8cpu <- pid_pvector_8cpu + 1
      trash <- gc()
      trash <- NULL
    },
    "ScMaSigPro_1_CPU_sc.t.fit" = {
      trash <- gc()
      trash <- NULL
      scmp.obj <- sc.t.fit(scmpObj = scmp.obj, verbose = FALSE, selection_method = "backward", parallel = FALSE, offset = T)
      cat(paste("\nFinished ScMaSigPro_1_CPU_sc.t.fit(", pid_tfit_1cpu, "for", i, "cells", "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
      pid_tfit_1cpu <- pid_tfit_1cpu + 1
      trash <- gc()
      trash <- NULL
    },
    "ScMaSigPro_8_CPU_sc.t.fit" = {
      trash <- gc()
      trash <- NULL
      scmp.obj <- sc.t.fit(scmpObj = scmp.obj, verbose = FALSE, selection_method = "backward", parallel = TRUE, n_cores = 8, offset = T)
      cat(paste("\nFinished ScMaSigPro_8_CPU_sc.t.fit(", pid_tfit_8cpu, "for", i, "cells", "):", timestamp(suffix = "", prefix = "", quiet = TRUE)))
      pid_tfit_8cpu <- pid_tfit_8cpu + 1
      trash <- gc()
      trash <- NULL
    },
    times = 10
  )

  cat(paste("\nBenchmark completed for", i, "cells :", timestamp(suffix = "", prefix = "", quiet = TRUE)))
  scmp.obj <- NULL
  rm(scmp.obj)



  # Extract Data for time
  data <- mbm %>% as.data.frame()

  # Add dataset size
  data$cells <- i

  # Rbibd
  master_data <- rbind(master_data, data)

  cat(paste("\nData added for ", i, "cells :", timestamp(suffix = "", prefix = "", quiet = TRUE)))
}

# Peform Cleanup

# Convert to minutes
master_data$time_minutes <- master_data$time / 1e9 / 60

# Add function
master_data$function_name <- str_split_i(master_data$expr, "ScMaSigPro_8_CPU_sc.|ScMaSigPro_1_CPU_sc.", 2)

# Add number of CPU
master_data$cpu <- str_split_i(str_split_i(master_data$expr, "ScMaSigPro\\_", 2), "\\_CPU", 1)

# Creating a violin plot for time
master_data <- master_data[mixedorder(master_data$cells), ]
master_data$cells <- factor(master_data$cells, levels = unique(master_data$cells))
master_data$cpu <- factor(master_data$cpu, levels = unique(master_data$cpu))

# Subset data
master_data_a <- master_data %>%
  filter(cpu == "1")
master_data_b <- master_data %>%
  filter(cpu == "8")

# Create Violin plot
v1 <- ggplot(master_data_a, aes(
  x = function_name,
  color = function_name,
  y = time_minutes, fill = function_name, shape = function_name
)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_jitter(width = 0.3, height = 0.3, alpha = 0.5) +
  labs(title = "Processing Times (1 CPU)", subtitle = "Evaluated 10 times", x = "", y = "Time (minutes)") +
  theme_minimal(base_size = 12) +
  facet_grid(. ~ cells) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "top") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1))

v2 <- ggplot(master_data_b, aes(
  x = function_name,
  color = function_name,
  y = time_minutes, fill = function_name, shape = function_name
)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_jitter(width = 0.3, height = 0.3, alpha = 0.5) +
  labs(title = "Processing Times (8 CPU)", subtitle = "Evaluated 10 times", x = "", y = "Time (minutes)") +
  theme_minimal(base_size = 12) +
  facet_grid(. ~ cells) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "top") +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5))

compare_violin <- ggarrange(v1, v2,
  ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom",
  labels = c("A", "B")
)

compare_violin

# Save the plot
ggsave(
  plot = compare_violin,
  filename = paste0(figPath_hd, "04_B_scmp_fns_time.png"),
  dpi = 600, width = 12, height = 10
)
ggsave(
  plot = compare_violin,
  filename = paste0(figPath_lr, "04_B_scmp_fns_time.png"),
  dpi = 150, width = 12, height = 10
)

# Save the plot
write.table(
  x = master_data,
  file = paste0(tabPath, "Individual_Function_Time_Profiling.txt"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

cat(paste("\nScript completed:", timestamp(suffix = "", prefix = "", quiet = TRUE)))
