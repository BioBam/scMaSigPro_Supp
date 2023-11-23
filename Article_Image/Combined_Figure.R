library(ggpubr)
# Load Data
a <- readRDS("Article_Image/Figure1_A.RDS")
b <- readRDS("Article_Image/Figure1_B.RDS")
<<<<<<< HEAD

combined_figure <- ggarrange(a,b, labels = c("A.", "B."))

ggsave(plot = combined_figure,
       path = "Article_Image",
       dpi = 1000,  filename = "MainFigure.png",
       width = 10, height = 5)
=======
c <- readRDS("Article_Image/Figure1_C.RDS")

combined_figure <- ggarrange(a,b,c, labels = c("A.", "B.", "C."),
                             ncol = 3, nrow = 1)
combined_figure
ggsave(plot = combined_figure,
       path = "Article_Image",
       dpi = 1000,  filename = "MainFigure.png",
       width = 18, height = 5)
>>>>>>> dev
