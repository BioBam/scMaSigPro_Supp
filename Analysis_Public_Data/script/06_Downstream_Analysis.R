<<<<<<< HEAD

# Extract the genes
scmp.obj <- sc.get.siggenes(scmpObj = scmp.obj, rsq = 0.7, vars = "groups")

# Plot Gene Expression
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "MPO", dis = scmp.obj@scTFit@dis,
              edesign =  scmp.obj@scTFit@edesign,
              groups.vector = scmp.obj@scTFit@groups.vector)

sc.PlotGroups(scmpObj = scmp.obj,feature_id = "EPOR", dis = scmp.obj@scTFit@dis,
              edesign =  scmp.obj@scTFit@edesign,
              groups.vector = scmp.obj@scTFit@groups.vector)


sc.PlotGroups(scmpObj = scmp.obj,feature_id = "IRF8", dis = scmp.obj@scTFit@dis,
              edesign =  scmp.obj@scTFit@edesign,
              groups.vector = scmp.obj@scTFit@groups.vector)

sc.PlotGroups(scmpObj = scmp.obj, feature_id = "EBF1", dis = scmp.obj@scTFit@dis,
              edesign =  scmp.obj@scTFit@edesign,
              groups.vector = scmp.obj@scTFit@groups.vector)

sc.plot.bins(scmpObj = scmp.obj)
sc.path.intersection(scmpObj = scmp.obj)
=======
# Load data
load("Analysis_Public_Data/data/SingleCellExperimentAtlas/scMaSigPro_Donor2_Monocle3.RData")

# Select Models with significant expression patterns per term
scmp.obj <- sc.get.siggenes(scmpObj = scmp.obj, rsq = 0.7, vars = "groups")

# Plot intersection
sc.path.intersection(scmp.obj)

sols <- showSol(scmp.obj, view = F, return = T, includeInflu = F)
sols[is.na(sols$`p-value`), "p-value"] <- 0
sols <- sols[sols$`p-value` <= 0.05, ]

sols <- sols[sols$`R-squared` >= 0.9,]

all_plots <- list()
# Plot Gene Expression
for (i in row.names(sols)){
    all_plots[[i]] <- sc.PlotGroups(scmpObj = scmp.obj,feature_id = i,
              logs = T,
              logType = "log")
}
ggsave(plot = ltb,
       filename = "Figure1_C.png",
       path = "Article_Image/",
       dpi = 1200, width = 5, height = 5)
saveRDS(
    ltb, "Article_Image/Figure1_C.RDS"
)

scmp.obj <- sc.cluster.features(scmpObj = scmp.obj,k = 9, cluster.method = "kmeans",
                                includeInflu = F)

cluster.trend.data <- sc.PlotProfiles(scmp.obj,
                groupBy = "coeff", logs = F)

cluster.trend.data <- cluster.trend.data[cluster.trend.data$feature_id %in% rownames(sols),]

ggplot(data = cluster.trend.data, aes(x = x, y = log(y), group = interaction(feature_id, path), color = path)) +
    geom_line(aes(linetype = path), size = 0.4) +
    geom_point(size = 1, shape = 21, alpha = 0.5, stroke = 1) +
    facet_wrap(~ cluster_id, scales = "free_y") +
    scale_color_manual(values = colorConesa(length(unique(cluster.trend.data$path)))) +
    theme_classic(base_size = 12) +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, angle = 0),
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = "Gene Expression over Pseudotime", color = "Path") +
    guides(color = guide_legend(title = "Path"))


>>>>>>> dev


# Get the sol
showTab <- showSol(scmpObj = scmp.obj, return = T, view = F)

View(showTab[rownames(showTab) %in% c("GATA1", "EPOR", "IRF8", "MPO"), ])
