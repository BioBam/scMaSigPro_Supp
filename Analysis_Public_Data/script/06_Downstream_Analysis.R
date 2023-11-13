
# Extract the genes
scmp.obj <- sc.get.siggenes(scmpObj = scmp.obj, rsq = 0.7, vars = "groups")

sc.path.intersection(scmp.obj)

sols <- showSol(scmp.obj, view = F, return = T, includeInflu = F)
sols[is.na(sols$`p-value`), "p-value"] <- 0
sols <- sols[sols$`p-value` <= 0.05, ]

sols <- sols[sols$`R-squared` >= 0.7,]

write.table(rownames(sols), "Analysis_Public_Data/data/rep3/sig.gene.tsv", sep = "\t",
            row.names = F, quote = F)

# Plot Gene Expression
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "SLC22A15",
              logs = T,
              logType = "log")
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "MYC",
              logs = T,
              logType = "log")
sc.PlotGroups(scmpObj = scmp.obj,feature_id ="COX7B",
              logs = T,
              logType = "log")
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "RAB44",
              logs = T,
              logType = "log")
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "RPP30",
              logs = T,
              logType = "log")
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "BLVRB",
              logs = T,
              logType = "log")

ltb <- sc.PlotGroups(scmpObj = scmp.obj,feature_id = "LTB",
              logs = T,
              logType = "log") + theme_classic(base_size = 10) +theme(
                  legend.box = "horizontal",
                  legend.direction = "horizontal",
                  legend.position = c(1, 0.5), legend.justification = c("right", "bottom"),
                  panel.grid.major = element_line(linewidth = 0.3, color = "lightgrey", linetype = "dotted"),
                  panel.grid.minor = element_line(linewidth = 0.1, color = "lightgrey", linetype = "dotted"))

ggsave(plot = ltb,
       filename = "Figure1_C.png",
       path = "Article_Image/",
       dpi = 1200, width = 5, height = 5)
saveRDS(
    ltb, "Article_Image/Figure1_C.RDS"
)

scmp.obj <- sc.cluster.features(scmpObj = scmp.obj,k = 4)

sc.PlotProfiles(scmp.obj)

# Get the sol
showTab <- showSol(scmpObj = scmp.obj, return = T, view = F)

View(showTab[rownames(showTab) %in% c("GATA1", "EPOR", "IRF8", "MPO"), ])
