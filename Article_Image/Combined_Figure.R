library(ggpubr)
# Load Data
a <- readRDS("Article_Image/Figure1_A.RDS")
b <- readRDS("Article_Image/Figure1_B.RDS")
c <- readRDS("Article_Image/Figure1_C.RDS")

combined_figure <- ggarrange(a,b,c, labels = c("A.", "B.", "C."),
                             ncol = 3, nrow = 1)

ggsave(plot = combined_figure,
       path = "Article_Image",
       dpi = 1000,  filename = "MainFigure.png",
       width = 15, height = 5)
