"
plot 1D SFS
"

library(ggplot2)
library(cowplot)
library(wesanderson)

sfs <- read.delim("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/migRate.0.01.noDel.sfs.txt", header = TRUE, sep = "\t")

plA <- ggplot(subset(sfs, generation == 3000), aes(freq_p1, density, fill=region))+
  theme_classic()+
  geom_bar(stat="identity",position="dodge")+
  xlim(c(0.5,5.5))+
  xlab("allele count")+ylab("density")+
  scale_fill_manual(values=wes_palette("Moonrise3"))+
  theme(legend.position = "NA")+theme(text = element_text(size=14))
plA

plB <- ggplot(subset(sfs, generation == 3600), aes(freq_p1, density, fill=region))+
  theme_classic()+
  geom_bar(stat="identity",position="dodge")+
  xlim(c(0.5,5.5))+
  xlab("allele count")+ylab("density")+
  scale_fill_manual(values=wes_palette("Moonrise3"), name="",labels=c("background","foreground"))+
  theme(text = element_text(size=14))
plB

plot_grid(plA,plB,labels=c("A","B"),rel_widths=c(1,1.6))
ggsave("~/Desktop/ECB/2DSFS_scan/sims_data_LHU/figures/migRate.0.01.noDel.sfs.pdf",height=2.5,width=7)

