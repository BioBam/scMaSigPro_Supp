library(ggpubr)
# Load Data
a <- readRDS("Article_Image/Figure1_A.RDS")
b <- readRDS("Article_Image/Figure1_B.RDS")

combined_figure <- ggarrange(a,b, labels = c("A.", "B."))

ggsave(plot = combined_figure,
       path = "Article_Image",
       dpi = 1000,  filename = "MainFigure.png",
       width = 10, height = 5)