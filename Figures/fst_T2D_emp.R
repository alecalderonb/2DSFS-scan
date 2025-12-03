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
