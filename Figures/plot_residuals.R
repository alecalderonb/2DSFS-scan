"
plot residuals
"
library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(cowplot)

residuals <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.01.noDel.2DSFS.residuals_poisson.txt")

lim <- max(abs(residuals$residual), na.rm = TRUE)

plC <- ggplot(
  subset(residuals, generation == 3000 & !(freq_p1 == 0 & freq_p2 == 0)),
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = residual)) +
  coord_fixed() +
  
  scale_fill_gradient2(
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
    x = "p1",
    y = "p2",
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
plC

plD <- ggplot(
  subset(residuals, generation == 3600 & !(freq_p1 == 0 & freq_p2 == 0)),
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = residual)) +
  coord_fixed() +
  
  scale_fill_gradient2(
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
    x = "p1",
    y = "p2",
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
plD

plot_grid(plC,plD,labels=c("A","B"),rel_widths=c(1,1), ncol = 2)
