library(ggpubr)

b<- readRDS("Figures/MainArticle/MainArticle_FigureB.RDS")
c<- readRDS("Figures/MainArticle/MainArticle_FigureC.RDS")
d<- readRDS("Figures/MainArticle/MainArticle_FigureD.RDS")


combine <- ggarrange(b,c,d, nrow = 1)
combine

ggsave(plot = combine, filename = "Figure1_BCD.png",
       path = "Figures/MainArticle/",
       width = 16, dpi = 600, height = 5)
