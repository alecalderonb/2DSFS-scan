library(ggplot2)
library(wesanderson)
library(cowplot)
library(viridis)

read.table("~/projects/2DSFS-scan/data/likelihood.locus1_NC_087104.1_8930103-9010574.txt")->data

pal <- wes_palette("Zissou1", 100, type = "continuous")
# heatmap is a local dataset
pl1<-ggplot(data, aes(x = V1, y = V2, fill = -V4)) +
  geom_tile() + 
  scale_fill_viridis(limits=c(-5,8.5),guide="none")+theme_classic()+xlab(expression(italic(f["uv"])))+ylab(expression(italic(f["bv"])))
  
  #coord_equal()
 

  
  
read.table("~/projects/2DSFS-scan/data/likelihood.locus2_NC_087088.1_15081852-15342710.txt")->data

pal <- wes_palette("Zissou1", 100, type = "continuous")
# heatmap is a local dataset
pl2<-ggplot(data, aes(x = V1, y = V2, fill = -V4)) +
  geom_tile() + 
  scale_fill_viridis(name=expression(-sigma * log(italic(L))),limits=c(-5,8.5))+theme_classic()+xlab(expression(italic(f["uv"])))+ylab(expression(italic(f["bv"])))+theme()

plot_grid(pl1,pl2,labels=c("A","B"),rel_widths=c(1,1.2))


