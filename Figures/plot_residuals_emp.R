"
plot residuals from empirical data
"

library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(cowplot)

# poisson
residuals_chr1_chrZ <- read.delim("~/Desktop/ECB/2DSFS_scan/emp_data/data/chr1_chrZ-4000000-4500000_residuals_poisson.txt")

lim <- max(abs(residuals_chr1_chrZ$residual), na.rm = TRUE)
plG <- ggplot(
  residuals_chr1_chrZ,
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = residual)) +
  coord_fixed() +
  scale_fill_gradient2(
    #trans = "log10",
    na.value = "transparent",
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish
  ) +
  theme_bw() +
  labs(
    x = "uv",
    y = "bv",
    fill = "residual"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    #axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    #legend.position = "none"
  )+
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.2)
plG
