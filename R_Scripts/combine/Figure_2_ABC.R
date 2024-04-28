library(ggpubr)

base_string <- "/supp_data/"
a <- readRDS(paste0(base_string, "additionalFigures/Figure2_A_to_D_Ground_Truth.rds"))
b <- readRDS(paste0(base_string, "additionalFigures/Figure2_E.rds"))
c <- readRDS(paste0(base_string, "additionalFigures/Figure2_F.rds"))


combine <- ggarrange(a, b, c,
  nrow = 1,
  labels = c("", "E", "F")
)
combine

ggsave(
  plot = combine, filename = "main_article_figure_2.png",
  path = "/supp_data/additionalFigures",
  width = 18, dpi = 600, height = 5
)
