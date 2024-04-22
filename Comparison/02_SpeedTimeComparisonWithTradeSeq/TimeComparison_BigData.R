# Title: Evaluation of BigData for SpeedTime Comparison
# Author: Priyansh Srivastava
# Year: 2024

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tradeSeq))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))


# Set paths
dirPath <- "/supp_data/ComparisonWithTradeSeq/simulated/sce/"
resPath <- "/supp_data/ComparisonWithTradeSeq/output/"
helpScriptsDir <- "R_Scripts/helper_function/"

# We will Simulate a BigData here
load(paste0(dirPath, "time_100k_cells.RData"))

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

# Release unused memory
sim.sce <- NULL
gc()

# Define core vector
use_cores <- seq(16, 24, 2)

# Corelist
core_list <- list()

# Compute bins
scmp.obj <- sc.squeeze(
  scmpObj = scmp.obj,
  drop_fac = 200,
  verbose = F,
  aggregate = "sum",
  split_bins = F,
  prune_bins = F,
  drop_trails = F,
  fill_gaps = F
)

# Set polynomial
scmp.obj <- sc.set.poly(scmp.obj, poly_degree = 3)

# Run evaluation
for (ncpu in use_cores) {
  cat(paste("\nRunning with", ncpu, "cores..."))


  mbm_scMaSigPro <- microbenchmark(
    "scMaSigPro::sc.p.vector()" = {
      # sc.p.vector
      scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj, verbose = T,
        min_na = 1,
        parallel = T,
        offset = T,
        max_it = 1000,
        n_cores = ncpu
      )
      gc()
    },
    "scMaSigPro::sc.t.fit()" = {
      # sc.t.fit
      scmp.obj <- sc.t.fit(
        scmpObj = scmp.obj, verbose = F,
        selection_method = "backward", parallel = T,
        offset = T,
        n_cores = ncpu
      )
      gc()
    },
    times = 5
  )

  # Extract Data
  core_list[[paste0("nCPUs_", ncpu)]] <- mbm_scMaSigPro
}

benchmark_data <- bind_rows(lapply(names(core_list), function(cpu_name) {
  results <- summary(core_list[[cpu_name]])
  transform(results,
    expr = as.character(expr), # Ensuring expression names are characters for plotting
    cpu = cpu_name
  ) # Adding CPU name for each set of results
}), .id = "cpu_id") # This adds an ID, but we already have cpu, so it's redundant


# Convert time from microseconds to minutes
benchmark_data$time <- paste(round(benchmark_data$mean / 60, digits = 6), "minutes")

# Convert 'cpu' column to a factor for ordered plotting
benchmark_data$cpu <- factor(benchmark_data$cpu, levels = unique(benchmark_data$cpu))

# Creating a line plot
ggplot(benchmark_data, aes(x = cpu, y = log10(median), color = expr, group = expr)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Median Execution Times Across CPU Configurations",
    x = "CPU Configuration",
    y = "Median Execution Time (ms)",
    color = "Function"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting the violin plot
ggplot(benchmark_data, aes(x = cpu, y = median, fill = expr)) +
  geom_violin(trim = FALSE) +
  labs(
    title = "Violin Plot of Execution Times by CPU Configuration",
    x = "CPU Configuration",
    y = "Median Execution Time (ms)",
    fill = "Function"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export Table
write.table(benchmark_data,
  file = paste(resPath, "100k_time_evaluation.txt",
    sep = "/"
  ),
  quote = FALSE, row.names = FALSE
)
saveRDS(core_list, paste(resPath, "100k_time_evaluation.RDS",
  sep = "/"
))

#
# Plotting the data
ggplot(benchmark_data, aes(x = expr, y = time, fill = cpu)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
  labs(x = "Function", y = "Evaluation Time (minutes)", title = "Benchmarking scMaSigPro Functions Across Different CPUs") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # +
scale_fill_brewer(palette = "Set2") # Adds color coding for

# compareViolins <- ggplot(data, aes(x = expr, y = mean, fill = expr)) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(
#     breaks = seq(0, 120, 20),
#     limits = c(0, 120)
#   ) +
#   labs(
#     title = "Execution Times for a bifurcating trajectory",
#     subtitle = "Number of Cells: 1500; Number of Genes: 1000",
#     x = "Method",
#     y = "Time (seconds)"
#   ) +
#   geom_text(aes(label = min_mean),
#     position = position_dodge(width = 0.9),
#     size = 3,
#     vjust = 0.5, hjust = -0.1
#   ) +
#   coord_flip() +
#   scale_fill_viridis(
#     discrete = TRUE, name = "Custom Legend Title",
#     breaks = c("TradeSeq_1_CPU", "ScMaSigPro_1_CPU", "TradeSeq_8_CPU", "ScMaSigPro_8_CPU"),
#     labels = c("Custom Label 1", "Custom Label 2", "Custom Label 3", "Custom Label 4")
#   ) +
#   theme_minimal(base_size = 20) +
#   theme(legend.position = "none", legend.justification = "left", legend.box.just = "left")
#
# compareBar_Time
#
# # Save
# ggsave(
#   plot = compareBar_Time,
#   filename = paste0("/supp_data/Figures/SuppData/04_tradeSeq_Time.png"),
#   dpi = 300, width = 10
# )
