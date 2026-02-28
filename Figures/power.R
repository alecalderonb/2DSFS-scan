library(ggplot2)
library(cowplot)
library(wesanderson)

#power relative to selected locus, not divergent

read.table("~/projects/2DSFS-scan/data/power.ll")->pow

plA<-ggplot(pow, aes(V2-3000,V3,group=V1,col=as.factor(V1))) + geom_smooth(method = "glm", 
                                                                           method.args = list(family = "binomial"), 
                                                                           se = FALSE) + geom_point(alpha=0.3) + theme_classic() + xlab("Generations from split") + ylab(expression("power, " * alpha * " = 5e-3"))  +  ylim(c(0,1))+ scale_color_manual(values=wes_palette("Darjeeling1"),name="migration rate")+theme(legend.position = "none")+geom_hline(yintercept = 0.8,lty=3)
)

read.table("~/projects/2DSFS-scan/data/power.fst")->pow2

plB<-ggplot(pow2, aes(V2-3000,V3,group=V1,col=as.factor(V1))) + geom_smooth(method = "glm", 
                                                                           method.args = list(family = "binomial"), 
                                                                           se = FALSE) + geom_point(alpha=0.3) + theme_classic() + xlab("Generations from split") + ylab(expression("power, " * alpha * " = 5e-3"))  +  ylim(c(0,1))+ scale_color_manual(values=wes_palette("Darjeeling1"),name="migration rate")+geom_hline(yintercept = 0.8,lty=3)


plot_grid(plA,plB,labels=c("A","B"),rel_widths=c(1,1.3))
ggsave("~/projects/2DSFS-scan/Figures/Fig2.pdf",height=2.5,width=8.5)

#power relative to neutral locus

read.table("~/projects/2DSFS-scan/data/power.ll")->pow

plA<-ggplot(pow, aes(V2-3000,V4,group=V1,col=as.factor(V1))) + geom_smooth(method = "gam", 
                                                                           se = FALSE) + geom_point(alpha=0.3) + theme_classic() + xlab("Generations from split") + ylab(expression("power, " * alpha * " = 5e-3"))  +  ylim(c(0,1))+ scale_color_manual(values=wes_palette("Darjeeling1"),name="migration rate")+theme(legend.position = "none")

read.table("~/projects/2DSFS-scan/data/power.fst")->pow2

plB<-ggplot(pow2, aes(V2-3000,V4,group=V1,col=as.factor(V1))) + geom_smooth(method = "gam", 
                                                                           se = FALSE) + geom_point(alpha=0.3) + theme_classic() + xlab("Generations from split") + ylab(expression("power, " * alpha * " = 5e-3"))  +  ylim(c(0,1))+ scale_color_manual(values=wes_palette("Darjeeling1"),name="migration rate")


plot_grid(plA,plB,labels=c("A","B"),rel_widths=c(1,1.3))

ggsave("~/projects/2DSFS-scan/Figures/FigPowerSupp.pdf",height=2.5,width=8.5)

# compare outliers at two time points

ll <- read.table('/Users/telemacher/projects/2DSFS-scan/sim_likelihood_calcs/migRate.0.05.3300.ll')
fst <-read.table('/Users/telemacher/projects/2DSFS-scan/sim_likelihood_calcs/migRate.0.05.3300_fst.txt',header=T)

comp<-data.frame(l=ll$V5,f=fst$avg_wc_fst,pos=fst$window_pos_1)
posFac = c()
for (posi in comp$pos) {
  if(posi > 490000 & posi < 500000) {
    posFac = c(posFac,"Background")
    
  } else if (posi > 990000 & posi < 1000000) {
    posFac = c(posFac,"Foreground")
  } else if (posi > 1490000 & posi < 1500000)  {
    posFac = c(posFac,"Background")
  } else {
    posFac = c(posFac,"Neutral")
  }
}

comp<-cbind(comp,data.frame(posFac=posFac))
llThresh<-sort(comp[comp$posFac == "Background",]$l)[0.995*length(comp[comp$posFac == "Background",]$l)]
fstThresh<-sort(comp[comp$posFac == "Background",]$f)[0.995*length(comp[comp$posFac == "Background",]$f)]

plA<-ggplot(comp,aes(l,f,col=as.factor(posFac)))+geom_point()+theme_classic()+ylab(expression(F[ST]))+xlab(expression(italic(T["2D"])))+theme(text = element_text(size=14))+
  scale_color_manual(values=wes_palette("AsteroidCity1"),name="Selection regime")
plA<-plA+geom_hline(yintercept=fstThresh,lty=3)
plA<-plA+geom_vline(xintercept=llThresh,lty=3)+theme(legend.position = "none")

cor.test(comp$l,comp$f)

# now a later generation 


ll <- read.table('/Users/telemacher/projects/2DSFS-scan/sim_likelihood_calcs/migRate.0.05.3600.ll')
fst <-read.table('/Users/telemacher/projects/2DSFS-scan/sim_likelihood_calcs/migRate.0.05.3600_fst.txt',header=T)

comp<-data.frame(l=ll$V5,f=fst$avg_wc_fst,pos=fst$window_pos_1)
posFac = c()
for (posi in comp$pos) {
  if(posi > 490000 & posi < 500000) {
    posFac = c(posFac,"Background")
    
  } else if (posi > 990000 & posi < 1000000) {
    posFac = c(posFac,"Foreground")
  } else if (posi > 1490000 & posi < 1500000)  {
    posFac = c(posFac,"Background")
  } else {
    posFac = c(posFac,"Neutral")
  }
}

comp<-cbind(comp,data.frame(posFac=posFac))
llThresh<-sort(comp[comp$posFac == "Background",]$l)[0.995*length(comp[comp$posFac == "Background",]$l)]
fstThresh<-sort(comp[comp$posFac == "Background",]$f)[0.995*length(comp[comp$posFac == "Background",]$f)]

plB<-ggplot(comp,aes(l,f,col=as.factor(posFac)))+geom_point()+theme_classic()+ylab(expression(F[ST]))+xlab(expression(italic(T["2D"])))+theme(text = element_text(size=14))+
  scale_color_manual(values=wes_palette("AsteroidCity1"),name="Selection regime")
plB<-plB+geom_hline(yintercept=fstThresh,lty=3)
plB<-plB+geom_vline(xintercept=llThresh,lty=3)

cor.test(comp$l,comp$f)

plot_grid(plA,plB,rel_widths = c(1,1.4),labels=c("A","B"))

ggsave("~/projects/2DSFS-scan/Figures/Fig.Comp.pdf",height=3.5,width=11)



