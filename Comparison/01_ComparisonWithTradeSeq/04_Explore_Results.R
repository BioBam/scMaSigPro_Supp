library(tidyverse)
library(UpSetR)
library(SingleCellExperiment)
library(ggpubr)

# Load the gene.info and the predictions
prediction.df <- readRDS("/supp_data/ComparisonWithTradeSeq/output/Prediction.RDS")
prediction.df$gene <- rownames(prediction.df)
colnames(prediction.df) <- c("ground_Truth", "TS_pattern", "TS_diffEnd", "scmp_0.6", "gene")
load("/supp_data/ComparisonWithTradeSeq/simulated/sce/testTradeSeq.RData")
gene.info <- rowData(sim.sce) %>% as.data.frame()
gene.info <- gene.info[, c("gene_short_name", "status", "status2")]
colnames(gene.info) <- c("gene", "DE", "foldChange")

# Combine
data <- merge(prediction.df, gene.info, "gene")

# Convert the columns to binary based on the p-value threshold (once)
p_value_threshold <- 0.05
data$TS_pattern <- ifelse(data$TS_pattern <= p_value_threshold, 1, 0)
data$TS_diffEnd <- ifelse(data$TS_diffEnd <= p_value_threshold, 1, 0)
data$scmp_0.6 <- ifelse(data$scmp_0.6 <= p_value_threshold, 1, 0)

# Create the UpSet plot
colnames(data) <- c(
  "GeneName",
  "Ground_Truth",
  "TradeSeq_pattern()",
  "TradeSeq_diffEnd()",
  "scMaSigPro_R2_0.06",
  "DE",
  "Fold_Change"
)

upset(
  data,
  sets = c("TradeSeq_pattern()", "TradeSeq_diffEnd()", "scMaSigPro_R2_0.06", "Ground_Truth"),
  main.bar.color = "#56B4E9", matrix.color = "#56B449",
  sets.bar.color = "#D55E00", order.by = "freq",
  text.scale = 2,
  number.angles = T,
  scale.intersections = "log2",
  keep.order = FALSE,
  show.numbers = TRUE, point.size = 3
)

# False Negative by all
fn_by_all <- data[
  (data$`TradeSeq_pattern()` == 0 &
    data$`TradeSeq_diffEnd()` == 0 &
    data$scMaSigPro_R2_0.06 == 0 &
    data$Ground_Truth == 1), ,
  drop = FALSE
]

# False Positive  by all
fp_by_all <- data[
  (data$`TradeSeq_pattern()` == 1 &
    data$`TradeSeq_diffEnd()` == 1 &
    data$scMaSigPro_R2_0.06 == 1 &
    data$Ground_Truth == 0), ,
  drop = FALSE
]

# Incorrect by scmp
incorrect_by_scmp <- data[
  (data$`TradeSeq_pattern()` == 1 &
    data$`TradeSeq_diffEnd()` == 1 &
    data$scMaSigPro_R2_0.06 == 0 &
    data$Ground_Truth == 1), ,
  drop = FALSE
]

# Incorrect by TS Pattern
incorrect_by_ts_pattern <- data[
  (data$`TradeSeq_pattern()` == 0 &
    # data$`TradeSeq_diffEnd()`== 1 &
    data$scMaSigPro_R2_0.06 == 1 &
    data$Ground_Truth == 1), ,
  drop = FALSE
]

# Incorrect by diff end
incorrect_by_ts_diffEnd <- data[
  ( # data$`TradeSeq_pattern()` == 0 &
    data$`TradeSeq_diffEnd()` == 0 &
      data$scMaSigPro_R2_0.06 == 1 &
      data$Ground_Truth == 1), ,
  drop = FALSE
]

df.list <- list(
  "fn_by_all" = fn_by_all,
  "incorrect_by_scmp" = incorrect_by_scmp,
  "incorrect_by_ts_pattern" = incorrect_by_ts_pattern,
  "incorrect_by_ts_diffEnd" = incorrect_by_ts_diffEnd
)

# Custom colors
custom_colors <- c("High_FC" = "#15918A", "Low_FC" = "#F58A53", "No_Change" = "#9F7BB8")

# Plot
bar.list <- lapply(names(df.list), function(df.name) {
  df.i <- df.list[[df.name]]

  # Calculate
  bar.df <- as.data.frame(table(df.i[, c("DE", "Fold_Change")]))
  bar.df <- bar.df[bar.df$Freq != 0, ]

  # Plot
  barplot <- ggplot(bar.df, aes(x = DE, y = Freq, fill = Fold_Change)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 5) +
    theme_minimal() +
    labs(x = "DE", y = "Frequency", fill = "Fold Change") +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = custom_colors)


  if (df.name == "fn_by_all") {
    barplot <- barplot + ggtitle("False Negative by all methods")
  } else if (df.name == "incorrect_by_scmp") {
    barplot <- barplot + ggtitle("False Negative by scMaSigPro",
      subtitle = "Correctly Identified by both 'TS_pattern' & 'TS_diffEnd'"
    )
  } else if (df.name == "incorrect_by_ts_pattern") {
    barplot <- barplot + ggtitle("False Negative by TS_pattern",
      subtitle = "Correctly Identified by scMaSigPro"
    )
  } else if (df.name == "incorrect_by_ts_diffEnd") {
    barplot <- barplot + ggtitle("False Negative by TS_diffEnd",
      subtitle = "Correctly Identified by scMaSigPro"
    )
  }
  return(barplot)
})

names(bar.list) <- names(df.list)

predictions.bar <- ggarrange(bar.list[["fn_by_all"]],
  bar.list[["incorrect_by_scmp"]],
  bar.list[["incorrect_by_ts_pattern"]],
  bar.list[["incorrect_by_ts_diffEnd"]],
  labels = c("A.", "B.", "C.", "D."),
  common.legend = F
)

predictions.bar

ggsave(
  plot = predictions.bar,
  filename = "/supp_data/Figures/SuppData/04_tradeSeq_bars.png",
  dpi = 600, width = 10, height = 8
)
