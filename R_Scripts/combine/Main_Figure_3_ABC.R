library(ggpubr)

base_string <- "/supp_data/"
a <- readRDS(paste0(base_string, "additionalFigures/Figure3_A.rds"))
b <- readRDS(paste0(base_string, "additionalFigures/Figure3_B.rds"))
c <- readRDS(paste0(base_string, "additionalFigures/Figure3_C.rds"))


combine <- ggarrange(a, b,
  nrow = 1,
  labels = c("A", "B"),
  widths = c(1.5, 2)
)
combine

ggsave(
  plot = combine, filename = "main_article_figure_3.png",
  path = "/supp_data/additionalFigures",
  width = 18, dpi = 600, height = 5
)
