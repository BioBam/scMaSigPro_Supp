library(ggpubr)

base_string <- "../scMaSigPro_supp_data/"
figPath <- paste0(base_string, "figures/")
figPath_hd <- paste0(figPath, "hd/")
figPath_lr <- paste0(figPath, "lr/")
a <- readRDS(paste0(base_string, "benchmarks/img/MainFigure2A.rds"))
b <- readRDS(paste0(base_string, "benchmarks/img/MainFigure_2B.rds"))
c <- readRDS(paste0(base_string, "comparison/out/MainFigure_2C.rds"))

# Update font family
a <- a + theme(text = element_text(family = "times"))
b <- b + theme(text = element_text(family = "times"))
c <- c + theme(text = element_text(family = "times"))

combine <- ggarrange(a, b, c,
  nrow = 1,
  labels = c("", "E.", "F.")
)
combine

ggsave(
  plot = combine, filename = "Main_Article_Figure_2.png",
  path = figPath_hd,
  width = 18, dpi = 600, height = 5,
  bg = "white"
)

ggsave(
  plot = combine, filename = "Main_Article_Figure_2.png",
  path = figPath_lr,
  width = 18, dpi = 300, height = 5,
  bg = "white"
)
