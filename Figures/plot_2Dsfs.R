"
plot 2D SFS from LHU sims
"
library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(cowplot)

# migration rate = 0.01
data <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.01.2DSFS.noDel.txt", header = TRUE, sep = "\t") 
data <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.01.2DSFS.noDel.BGloweffects.txt", header = TRUE, sep = "\t") 
data <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.01.noDel.2DSFS.txt", header = TRUE, sep = "\t") 
data <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.01.2DSFS.txt", header = TRUE, sep = "\t") 

# migration rate = 0.0
data <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.0.noDel.BGloweffects.txt", header = TRUE, sep = "\t") 

n <- 20

folded_sfs <- data %>%
  mutate(
    freq_p1_f = ifelse(freq_p1 + freq_p2 > n, n - freq_p1, freq_p1),
    freq_p2_f = ifelse(freq_p1 + freq_p2 > n, n - freq_p2, freq_p2)
  ) %>%
  group_by(region, generation, freq_p1_f, freq_p2_f) %>%
  summarise(density = sum(density), .groups = "drop") %>%
  group_by(region, generation) %>%
  mutate(density = density / sum(density)) %>%
  ungroup()

folded_sfs <- folded_sfs %>%
  filter(freq_p1_f + freq_p2_f <= n)

plA <- ggplot(
  subset(folded_sfs, region == "bg" & generation == 3000 & !(freq_p1_f == 0 & freq_p2_f == 0)
         ),
  aes(freq_p1_f, freq_p2_f)
  #aes(freq_p1, freq_p2)
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
    legend.position = "none"
  )+
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.2)

plB <- ggplot(
  subset(folded_sfs, region == "fg" & generation == 3000 & !(freq_p1_f == 0 & freq_p2_f == 0)
         ),
  aes(freq_p1_f, freq_p2_f)
  #aes(freq_p1, freq_p2)
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
ggsave("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/figures/migRate.0.0.3300.noDel.2DSFS.BGloweffects.pdf",height=2.5,width=7)

# get 2DSFS for background but only the locus with low effects (i think it was g2)
folded_sfs_bg <- folded_sfs %>%
  filter(region == "bg")
