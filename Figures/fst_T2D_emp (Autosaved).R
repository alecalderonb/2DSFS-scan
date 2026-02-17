# real pruned data

library(ggplot2)
library(cowplot)
library(wesanderson)

read.table("~/projects/2DSFS-scan/data/windows.250.txt")->data
zchr <- c()
for (i in seq(1,length(data$V2))) {
  z = "Autosome/W"
  if (data$V1[i] == "NC_087119.1") {
    z = "Z"
  }
  zchr<-c(zchr,z)
}

data<-cbind(data,data.frame(zchr = zchr))


plA<-ggplot(data,aes(V4,V5,col=as.factor(zchr)))+geom_point(size=1.7)+theme_classic()+xlab(expression(italic(T["2D"])))+theme(text=element_text(size=14))+scale_color_manual(values=wes_palette("AsteroidCity3")[3:4],name="Chromosome")+
  ylab(expression(italic(F["ST"])))+geom_hline(yintercept=sort(data$V5)[0.999*length(data$V4)],lty=3)+geom_vline(xintercept=sort(data$V4)[0.999*length(data$V3)],lty=3)+geom_hline(yintercept=sort(data$V5)[0.995*length(data$V4)],lty=2)+
  geom_vline(xintercept=sort(data$V4)[0.995*length(data$V3)],lty=2)+theme(legend.position = "none")

read.table("~/projects/2DSFS-scan/data/windows.500.txt")->data
zchr <- c()
for (i in seq(1,length(data$V2))) {
  z = "Autosome/W"
  if (data$V1[i] == "NC_087119.1") {
    z = "Z"
  }
  zchr<-c(zchr,z)
}

data<-cbind(data,data.frame(zchr = zchr))

plB<-ggplot(data,aes(V4,V5,col=as.factor(zchr)))+geom_point(size=1.7)+theme_classic()+xlab(expression(italic(T["2D"])))+theme(text=element_text(size=14))+scale_color_manual(values=wes_palette("AsteroidCity3")[3:4],name="Chromosome")+
  ylab(expression(italic(F["ST"])))+geom_hline(yintercept=sort(data$V5)[0.999*length(data$V4)],lty=3)+geom_vline(xintercept=sort(data$V4)[0.999*length(data$V3)],lty=3)+geom_hline(yintercept=sort(data$V5)[0.995*length(data$V4)],lty=2)+
  geom_vline(xintercept=sort(data$V4)[0.995*length(data$V3)],lty=2)+theme(legend.position = "none")


read.table("~/projects/2DSFS-scan/data/windows.1000.txt")->data
zchr <- c()
for (i in seq(1,length(data$V2))) {
  z = "Autosome/W"
  if (data$V1[i] == "NC_087119.1") {
    z = "Z"
  }
  zchr<-c(zchr,z)
}

data<-cbind(data,data.frame(zchr = zchr))

plC<-ggplot(data,aes(V4,V5,col=as.factor(zchr)))+geom_point(size=1.7)+theme_classic()+xlab(expression(italic(T["2D"])))+theme(text=element_text(size=14))+scale_color_manual(values=wes_palette("AsteroidCity3")[3:4],name="Chromosome")+
  ylab(expression(italic(F["ST"])))+geom_hline(yintercept=sort(data$V5)[0.999*length(data$V4)],lty=3)+geom_vline(xintercept=sort(data$V4)[0.999*length(data$V3)],lty=3)+geom_hline(yintercept=sort(data$V5)[0.995*length(data$V4)],lty=2)+
  geom_vline(xintercept=sort(data$V4)[0.995*length(data$V3)],lty=2)


plot_grid(plA,plB,plC,labels=c("A","B","C"),ncol=3,rel_widths=c(1,1,1.5))
ggsave("~/projects/2DSFS-scan/Figures/FST_T2D.pdf",height=3,width=12)


##### Coding 250 windows


read.table("~/projects/2DSFS-scan/data/T2D_and_FST.coding.250.txt")->data
zchr <- c()
for (i in seq(1,length(data$V2))) {
  z = "Autosome/W"
  if (data$V1[i] == "NC_087119.1") {
    z = "Z"
  }
  zchr<-c(zchr,z)
}

data<-cbind(data,data.frame(zchr = zchr))

plCoding<-ggplot(data,aes(V4,V5,col=as.factor(zchr)))+geom_point(size=1.7)+theme_classic()+xlab(expression(italic(T["2D"])))+theme(text=element_text(size=14))+scale_color_manual(values=wes_palette("AsteroidCity3")[3:4],name="Chromosome")+
  ylab(expression(italic(F["ST"])))+geom_hline(yintercept=sort(data$V5)[0.99*length(data$V4)],lty=2)+
  geom_vline(xintercept=sort(data$V4)[0.99*length(data$V3)],lty=2)+xlim(150,450)

plCoding<-plCoding+# Add the arrow (segment with an arrowhead)
  annotate("segment", 
           x = 415, xend = 380, 
           y = 0.1, yend = 0.106,
           colour = "black", 
           arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + # 'type="closed"' fills the arrowhead
  # Add the text label
  annotate("text", 
           x = 425, y = 0.1, 
           label = expression(italic("per")), 
           color = "black", 
           size = 5,
           fontface = "bold") #+
   #coord_cartesian(clip = "off") # if you need annotations outside the plot area
  #coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
  
 plCoding<-plCoding+# Add the arrow (segment with an arrowhead)
  annotate("segment", 
           x = 320, xend = 300, 
           y = 0.075, yend = 0.066,
           colour = "black", 
           arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + # 'type="closed"' fills the arrowhead
  # Add the text label
  annotate("text", 
           x = 332, y = 0.075, 
           label = expression(italic("pdfr")), 
           color = "black", 
           size = 5,
           fontface = "bold") #+
   #coord_cartesian(clip = "off") # if you need annotations outside the plot area
  #coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
  
   plCoding<-plCoding+# Add the arrow (segment with an arrowhead)
  annotate("segment", 
           x = 320, xend = 297, 
           y = 0.022, yend = 0.014,
           colour = "black", 
           arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + # 'type="closed"' fills the arrowhead
  # Add the text label
  annotate("text", 
           x = 340, y = 0.023, 
           label = expression(italic("rho1")), 
           color = wes_palette("AsteroidCity3")[1], 
           size = 5,
           fontface = "bold") #+
   #coord_cartesian(clip = "off") # if you need annotations outside the plot area
  #coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
   plCoding<-plCoding+# Add the arrow (segment with an arrowhead)
  annotate("segment", 
           x = 320, xend = 297, 
           y = 0.003, yend = 0.007,
           colour = "black", 
           arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + # 'type="closed"' fills the arrowhead
  # Add the text label
  annotate("text", 
           x = 340, y = 0.002, 
           label = expression(italic("prl-1")), 
           color = wes_palette("AsteroidCity3")[1], 
           size = 5,
           fontface = "bold") #+
   #coord_cartesian(clip = "off") # if you need annotations outside the plot area
  #coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))




ggsave("~/projects/2DSFS-scan/Figures/FST_T2D_coding_250.pdf",height=4/1.5,width=10/1.5)

