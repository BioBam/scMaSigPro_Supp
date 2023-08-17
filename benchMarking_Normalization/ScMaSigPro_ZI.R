# Load the library
library(scMaSigPro)

# Load the input 
sce <- readRDS("benchMarking_Normalization/data/output/SimObjectsSCE/zi_shape_-0.2.RDS")

# Select Random Genes
random_genes <- sample(1:5000, 2000, replace = FALSE)

# Subset
sce.sub <- sce[random_genes, ]
table(rowData(sce.sub)$status)

# Create Design
scmp.obj <- sc_makeDesign(sce.sub)

# P-vector ZI
scmp.obj <- sc_pVector(scmp.obj, 
                       model.type = "zip")

# Tfit ZI
scmp.obj <- sc_tFit(scmp.obj, covar.sig = 0.05)

# Extract Solution
sol <- showSol(scmp.obj, view = F)
saveRDS(sol, "benchMarking_Normalization/data/output/ZI/sol.RDS")

row.data <- as.data.frame(rowData(sce.sub))
all <- rownames(row.data)
gene.change <- rownames(row.data[row.data$status != "No_Change", ])
gene.no.change <- rownames(row.data[row.data$status == "No_Change", ])

# Running the test 
# Create Dataframe to plot Performance
df.performance <- data.frame(VARIABLE = 0, TP = 0, FP = 0, TN = 0, FN = 0,
                             TPR = 0, FPR = 0, FNR = 0, TNR = 0, 
                             ACCURACY = 0)

sol.pval <- sol[sol[,1] <= 0.05, ]

# Run For all Values of R square
for (j in seq(0.05, 0.95, 0.05)){
    
    sol.test <- sol.pval[sol.pval[,2] >= j,]
    
    detected.genes <- rownames(sol.test)
    undetected.gene <- all[!(all %in% detected.genes)]
    
    # True Positive
    tp <- intersect(gene.change, detected.genes) 
    
    # True Negative
    tn <- intersect(gene.no.change, undetected.gene)
    
    # False Positive
    fp <- intersect(gene.no.change, detected.genes) 
    
    # False Negative
    fn <- intersect(gene.change, undetected.gene)
    
    # Sensitivity, True Positive Rate / Recall
    senstivity <- length(tp)/(length(tp) + length(fn))
    
    # False Negative Rate
    fnr <- length(fn)/(length(tp) + length(fn))
    
    # Specificity / True Negative Rate
    specificity <- length(tn)/(length(tn) + length(fp))
    
    # False Positive Rate
    fpr <- 1- specificity
    
    # Accuracy
    accuracy <- (length(tp)+length(tn))/(length(tp)+length(tn)+length(fp)+length(fn))
    
    # Make the vector of results
    res.perfom <- data.frame(VARIABLE = j,
                             TP = length(tp),
                             FP = length(fp),
                             TN = length(tn),
                             FN = length(fn),
                             TPR = round(senstivity,3),
                             FPR = round(fpr,3),
                             FNR = round(fnr,3),
                             TNR = round(specificity,3),
                             ACCURACY = round(accuracy, 3)
    )
    # Add vector to the frame
    df.performance <- rbind(df.performance, res.perfom)
    gc()
}

# Remove the first Row
df.performance <- df.performance[-1,]

saveRDS(df.performance, "benchMarking_Normalization/data/output/ZI/df.performance.zinb.RDS")

sol <- readRDS("benchMarking_Normalization/data/output/ZI/sol.zip.RDS")




a <- ggplot(df.performance, aes(x = FPR, y = TPR))+
    geom_point(size = 2, alpha = 0.7, color = "red")+
    geom_path(size = 1, alpha = 0.7, color = "blue")+
    scale_x_continuous(breaks = seq(0,0.5, 0.05))+
    #scale_y_continuous(breaks = seq(0.5,1, 0.05))+
    labs(title= "ROC-curve, Different Values of R-Square", 
         x = "False Positive Rate (1-Specificity)", 
         y = "True Positive Rate (Sensitivity)") + theme_classic()+
    #geom_text(data = df.performance, aes(x = FPR, y = TPR, label = VARIABLE), 
     #         nudge_y = 0.02, size = 4, nudge_x =0) +
    theme(panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
          panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
          axis.text = element_text(size = rel(1.5)),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
    ) + geom_vline(xintercept = 0.01, colour= "green") + geom_vline(xintercept = 0.05, colour= "orange")  + geom_vline(xintercept = 0.1, colour= "red")

# GGPlot2
b<-ggplot(df.performance, aes(x = VARIABLE, y = ACCURACY))+
    geom_point(size = 2, alpha = 0.7, color = "red")+
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.05))+
    geom_path(size = 1, alpha = 0.7, color = "blue")+
    #scale_y_continuous(breaks =  seq(0.7, 1, 0.05))+
    labs(title= "Accuracy Against Chnaging R Square",
         subtitle = "Red.Dots: False Negatives",
         x = "Increasing Value for R-Square", 
         y = "Accuracy") + theme_classic()+
    #geom_text(data = df.performance, aes(x = VARIABLE, y = ACCURACY, label = FN), 
     #         nudge_y= 0.01, nudge_x = 0.01, size = 3) +
    theme(panel.grid.major = element_line(size = 0.7, color = "lightgrey", linetype = "dotted"),
          panel.grid.minor = element_line(size = 0.2, color = "grey", linetype = "dotted"),
          axis.text = element_text(size = rel(1.5)),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
    )

