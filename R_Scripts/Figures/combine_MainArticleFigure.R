library(ggpubr)

b <- readRDS("/supp_data/Figures/MainArticle/MainArticle_FigureB.RDS")
c <- readRDS("/supp_data/Figures/MainArticle/MainArticle_FigureC.RDS")
d <- readRDS("/supp_data/Figures/MainArticle/MainArticle_FigureD.RDS")


combine <- ggarrange(b, c, d, nrow = 1)
combine

ggsave(
  plot = combine, filename = "Figure1_BCD.png",
  path = "/supp_data/Figures/MainArticle/",
  width = 16, dpi = 400, height = 5
)
