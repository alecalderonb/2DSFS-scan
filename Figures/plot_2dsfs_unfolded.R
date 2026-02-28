"
plot 2D SFS unfolded
"

library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(cowplot)

# migration rate = 0.01
data <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.01.2DSFS.noDel.txt", header = TRUE, sep = "\t") 

plA <- ggplot(
  subset(data, region == "bg" & generation == 3000 & !(freq_p1 == 0 & freq_p2 == 0)
  ),
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = log10(density))) +
  coord_fixed() +
  scale_fill_viridis_c(na.value = "transparent",
                       labels = scales::label_math(10^.x),
                       direction = -1
                       #name = expression(log[10](density))
  ) +
  #guides(fill = guide_colorbar(reverse = TRUE))+
  #theme_ipsum()+
  theme_bw()+
  labs(#title="bg - 3000",
    x="p1",
    y="p2") +
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

plB <- ggplot(
  subset(data, region == "fg" & generation == 3000 & !(freq_p1 == 0 & freq_p2 == 0)
  ),
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = log10(density))) +
  coord_fixed() +
  scale_fill_viridis_c(na.value = "transparent",
                       labels = scales::label_math(10^.x),
                       direction = -1,
                       #name = expression(log[10](density))
                       name = NULL
  ) +
  #guides(fill = guide_colorbar(reverse = TRUE))+
  #theme_ipsum()+
  theme_bw()+
  labs(#title="bg - 3000",
    x="p1",
    y="p2") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    #axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )+
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.2)

plC <- ggplot(
  subset(data, region == "bg" & generation == 3600 & !(freq_p1 == 0 & freq_p2 == 0)
  ),
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = log10(density))) +
  coord_fixed() +
  scale_fill_viridis_c(na.value = "transparent",
                       labels = scales::label_math(10^.x),
                       direction = -1
                       #name = expression(log[10](density))
  ) +
  #guides(fill = guide_colorbar(reverse = TRUE))+
  #theme_ipsum()+
  theme_bw()+
  labs(#title="bg - 3000",
    x="p1",
    y="p2") +
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

plD <- ggplot(
  subset(data, region == "fg" & generation == 3600 & !(freq_p1 == 0 & freq_p2 == 0)
  ),
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = log10(density))) +
  coord_fixed() +
  scale_fill_viridis_c(na.value = "transparent",
                       labels = scales::label_math(10^.x),
                       direction = -1,
                       #name = expression(log[10](density))
                       name = NULL
  ) +
  #guides(fill = guide_colorbar(reverse = TRUE))+
  #theme_ipsum()+
  theme_bw()+
  labs(#title="bg - 3000",
    x="p1",
    y="p2") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    #axis.ticks = element_line(linewidth = 0.6),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )+
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.2)

plot_grid(plA,plB,labels=c("A","B"),rel_widths=c(1,1))
plot_grid(plC,plD,labels=c("A","B"),rel_widths=c(1,1))

