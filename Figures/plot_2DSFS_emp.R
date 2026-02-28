"
plot 2D SFS for empirical data
"

library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(cowplot)

data_chr1 <- read.delim("~/Desktop/ECB/2DSFS_scan/emp_data/data/chr1.sfs.txt", header = TRUE, sep = "\t")
data <- read.delim("~/Desktop/ECB/2DSFS_scan/emp_data/data/chrZ_4000000-4500000.sfs.txt", header = TRUE, sep = "\t")

# chromosome 1
lim <- max(c(data_chr1$density, data$density), na.rm = TRUE)
min_positive <- min(
  c(data_chr1$density[data_chr1$density > 0],
    data$density[data$density > 0]),
  na.rm = TRUE
)

plA <- ggplot(
  subset(data_chr1, !(freq_p1 == 0 & freq_p2 == 0)),
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = density)) +
  coord_fixed() +
  scale_fill_viridis_c(
    trans = "log10",
    #limits = c(min_positive, lim),
    #oob = scales::squish,
    #labels = scales::label_math(10^.x),
    direction = -1,
    na.value = "transparent"
  )+
  theme_bw() +
  labs(x = "uv", y = "bv") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    #legend.position = "none"
  ) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.2)
plA


# chr Z window with highest FST
plD <- ggplot(
  data,
  aes(freq_p1, freq_p2)
) +
  geom_tile(aes(fill = density)) +
  coord_fixed() +
  scale_fill_viridis_c(
    trans = "log10",
    limits = c(min_positive, lim),
    oob = scales::squish,
    #labels = scales::label_math(10^.x),
    direction = -1,
    na.value = "transparent"
  )+
  theme_bw() +
  labs(x = "uv", y = "bv") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.position = "none"
  ) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.2)
plD
