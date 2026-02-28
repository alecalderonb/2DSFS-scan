"
manhattan plots - ECB data
"

library(ggplot2)
library(cowplot)
library(wesanderson)

read.table("~/Desktop/ECB/2DSFS_scan/emp_data/windows.500.txt")->data
zchr <- c()
for (i in seq(1,length(data$V2))) {
  z = "Autosome/W"
  if (data$V1[i] == "NC_087119.1") {
    z = "Z"
  }
  zchr<-c(zchr,z)
}

data<-cbind(data,data.frame(zchr = zchr))

mapping <- read.table("~/Desktop/ECB/2DSFS_scan/emp_data/chromosomes.txt", header = TRUE, stringsAsFactors = FALSE)
data <- data %>%
  left_join(mapping, by = c("V1" = "chr_id")) %>%
  mutate(chr_id = chr_num) %>%
  select(-chr_num)

# remove NA chr IDs -> unmapped chromosomes
data <- data %>% drop_na()

data_manhattan <- data %>% 
  
  # Compute chromosome size
  group_by(chr_id) %>% 
  summarise(chr_len=max(V2)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(data, ., by=c("chr_id"="chr_id")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr_id, V2) %>%
  mutate( BPcum=V2+tot)

axisdf = data_manhattan %>% group_by(chr_id) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

ggplot(subset(data_manhattan, V5>0), aes(x=BPcum, y=V5)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr_id)), alpha=1, size=1.3) +
  #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  # scale_x_continuous( label = axisdf$chr_id, breaks= axisdf$center ) +
  scale_x_continuous(
    breaks = axisdf$center,
    labels = ifelse(seq_along(axisdf$chr_id) %% 2 == 1,
                    axisdf$chr_id,
                    "")
  )+
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    # panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  labs(x="Chromosome", y=expression(italic(F["ST"])))+
  coord_cartesian(
    ylim = c(-0.01, 0.15)
  ) +
  scale_color_manual(
    values = rep(c(
      wes_palette("AsteroidCity2")[5],
      wes_palette("AsteroidCity3")[4]
    ), length.out = length(unique(data_manhattan$chr_id)))
  )+
  geom_hline(yintercept=sort(data_manhattan$V5)[0.995*length(data_manhattan$V5)],lty=2)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))


ggsave("~/Desktop/ECB/2DSFS_scan/emp_data/manhattan_500.pdf",height=2.5,width=8.5)



