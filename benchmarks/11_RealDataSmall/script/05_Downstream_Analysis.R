
# Extract the genes
scmp.obj <- sc.get.siggenes(scmpObj = scmp.obj, rsq = 0.7, vars = "groups")

# Plot Gene Expression
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "MPO", dis = scmp.obj@scTFit@dis,
              edesign =  scmp.obj@scTFit@edesign,
              groups.vector = scmp.obj@scTFit@groups.vector)

sc.plot.bins(scmpObj = scmp.obj)
sc.path.intersection(scmpObj = scmp.obj)


# Get the sol
showTab <- showSol(scmpObj = scmp.obj, return = T, view = F)

View(showTab[rownames(showTab) %in% c("GATA1", "EPOR", "IRF8", "MPO"), ])
