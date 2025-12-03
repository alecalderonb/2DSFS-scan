library(ggplot2)
library(cowplot)
library(wesanderson)

read.table("~/projects/2DSFS-scan/data/sfs.sim.txt")->sfs

plC<-ggplot(sfs[sfs$V2 == 3000,],aes(V3,V4,fill=V1))+theme_classic()+geom_bar(stat="identity",position="dodge")+xlim(c(0.5,5))+xlab("allele count")+ylab("Density")+scale_fill_manual(values=wes_palette("Moonrise3"))+theme(legend.position = "NA")+theme(text = element_text(size=14))
plD<-ggplot(sfs[sfs$V2 == 3300,],aes(V3,V4,fill=V1))+theme_classic()+geom_bar(stat="identity",position="dodge")+xlim(c(0.5,5))+xlab("allele count")+ylab("Density")+scale_fill_manual(values=wes_palette("Moonrise3"),name="",labels=c("background","foreground"))+theme(text = element_text(size=14))

read_csvs<-function(iter) {
    connection <- gzfile(paste("~/projects/2DSFS-scan/slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams/migRate.0.01/iter_",iter,"/sim_log.",iter,".txt.gz",sep=""), "rt") # "rt" for read text mode
    data <- read.csv(connection, header = TRUE)
    data<-cbind(data,data.frame(i=rep(iter,length(data$cycle))))
    close(connection)
    return(data)
}

data<-data.frame()

for (i in seq(1,1000)) {
  data<-rbind(data,read_csvs(i))
}


plA<-ggplot(data,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                    fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                       fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


plTopRow<-plot_grid(plA,plB,labels=c("A","B"))


plBottomRow<-plot_grid(plC,plD,labels=c("C","D"),rel_widths=c(1,1.4 ))

plot_grid(plTopRow,plBottomRow,ncol=1)

ggsave("~/projects/2DSFS-scan/Figures/Fig1.pdf",width=10,height=7)

# supp figure of different migration rates


read_csvs_mig<-function(iter,mig) {
  connection <- gzfile(paste("~/projects/2DSFS-scan/slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams/migRate.",as.character(mig),"/iter_",iter,"/sim_log.",iter,".txt.gz",sep=""), "rt") # "rt" for read text mode
  data <- read.csv(connection, header = TRUE)
  data<-cbind(data,data.frame(i=rep(iter,length(data$cycle))))
  close(connection)
  return(data)
}

data_mig0.0<-data.frame()
data_mig0.01<-data.frame()
data_mig0.05<-data.frame()
data_mig0.1<-data.frame()
data_mig0.2<-data.frame()

for (i in seq(1,1000)) {
  data_mig0.0<-rbind(data_mig0.0,read_csvs_mig(i,"0.0"))
  data_mig0.01<-rbind(data_mig0.01,read_csvs_mig(i,0.01))
  data_mig0.05<-rbind(data_mig0.05,read_csvs_mig(i,0.05))
  data_mig0.1<-rbind(data_mig0.1,read_csvs_mig(i,0.1))
  data_mig0.2<-rbind(data_mig0.2,read_csvs_mig(i,0.2))
}

# now plots for different mig_rates
plA<-ggplot(data_mig0.0,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                    fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.0,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                       fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.0,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.0,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p0<-plot_grid(plA,plB,labels=c("A","B"))

plA<-ggplot(data_mig0.01,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                           fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.01,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                              fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.01,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.01,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p01<-plot_grid(plA,plB,labels=c("C","D"))

plA<-ggplot(data_mig0.05,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                            fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.05,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                               fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.05,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.05,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p05<-plot_grid(plA,plB,labels=c("E","F"))

plA<-ggplot(data_mig0.1,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                            fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.1,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                               fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.1,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.1,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p1<-plot_grid(plA,plB,labels=c("G","H"))

plA<-ggplot(data_mig0.2,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                            fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.2,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                               fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.2,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.2,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p2<-plot_grid(plA,plB,labels=c("I","J"))


plot_grid(pl0p0,pl0p01,pl0p05,pl0p1,pl0p2,ncol=1)

ggsave("~/projects/2DSFS-scan/Figures/FigPMeanFst_withDel.pdf",width=8.5*1.3,height=11*1.3)

# figure of change in genetic variance over time


data_mig0.0<-cbind(data_mig0.0,data.frame(mig=rep(0.0,length(data_mig0.0$cycle))))
data_mig0.01<-cbind(data_mig0.01,data.frame(mig=rep(0.01,length(data_mig0.01$cycle))))
data_mig0.05<-cbind(data_mig0.05,data.frame(mig=rep(0.05,length(data_mig0.05$cycle))))
data_mig0.1<-cbind(data_mig0.1,data.frame(mig=rep(0.1,length(data_mig0.1$cycle))))
data_mig0.2<-cbind(data_mig0.2,data.frame(mig=rep(0.2,length(data_mig0.2$cycle))))

dataAll<-rbind(data_mig0.0,data_mig0.01,data_mig0.05,data_mig0.1,data_mig0.2)

plMigSDWithDel<-ggplot(dataAll,aes(cycle-3000,SD,color=as.factor(mig)))+
  stat_summary(geom="line", fun="mean",size=1)+theme_classic()+ylab(expression(sigma))+xlab("generations after shift")+theme(text = element_text(size=14))+
  scale_color_manual(values=wes_palette("AsteroidCity1"),name="migration rate")+theme(legend.position = "none")



# supp figure of different migration rates, no del


read_csvs_mig<-function(iter,mig) {
  connection <- gzfile(paste("~/projects/2DSFS-scan/slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams_noDel/migRate.",as.character(mig),"/iter_",iter,"/sim_log.",iter,".txt.gz",sep=""), "rt") # "rt" for read text mode
  data <- read.csv(connection, header = TRUE)
  data<-cbind(data,data.frame(i=rep(iter,length(data$cycle))))
  close(connection)
  return(data)
}

data_mig0.0<-data.frame()
data_mig0.01<-data.frame()
data_mig0.05<-data.frame()
data_mig0.1<-data.frame()
data_mig0.2<-data.frame()

for (i in seq(1,1000)) {
  data_mig0.0<-rbind(data_mig0.0,read_csvs_mig(i,"0.0"))
  data_mig0.01<-rbind(data_mig0.01,read_csvs_mig(i,0.01))
  data_mig0.05<-rbind(data_mig0.05,read_csvs_mig(i,0.05))
  data_mig0.1<-rbind(data_mig0.1,read_csvs_mig(i,0.1))
  data_mig0.2<-rbind(data_mig0.2,read_csvs_mig(i,0.2))
}

# now plots for different mig_rates
plA<-ggplot(data_mig0.0,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                           fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.0,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                              fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.0,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.0,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p0<-plot_grid(plA,plB,labels=c("A","B"))

plA<-ggplot(data_mig0.01,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                            fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.01,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                               fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.01,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.01,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p01<-plot_grid(plA,plB,labels=c("C","D"))

plA<-ggplot(data_mig0.05,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                            fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.05,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                               fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.05,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.05,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p05<-plot_grid(plA,plB,labels=c("E","F"))

plA<-ggplot(data_mig0.1,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                           fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.1,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                              fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.1,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.1,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p1<-plot_grid(plA,plB,labels=c("G","H"))

plA<-ggplot(data_mig0.2,aes(cycle-3000,Mean))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                           fun.args=list(mult = 2), fill="gray")+
  stat_summary(geom="line", fun.y=mean,size=1,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab("phenotype mean")+xlab("generations after shift")+theme(text = element_text(size=14))

plA<-plA+geom_line(data=data.frame(x=seq(1,600),y=1+0.002*seq(1,600)),aes(x,y),lty=3,col="black")

plB<-ggplot(data_mig0.2,aes(cycle-3000,FSTNeut))+stat_summary(geom="ribbon", fun.data="mean_sdl", 
                                                              fun.args=list(mult = 3), fill="gray",alpha=0.5)+
  stat_summary(geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[3])+theme_classic()+ylab(expression(F[ST]))+xlab("generations after shift")+theme(text = element_text(size=14))
plB<-plB+stat_summary(data=data_mig0.2,aes(cycle-3000,FSTBg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[1])
plB<-plB+stat_summary(data=data_mig0.2,aes(cycle-3000,FSTFg),geom="line", fun.y=mean,size=0.7,col=wes_palette("Moonrise3")[2])


pl0p2<-plot_grid(plA,plB,labels=c("I","J"))


plot_grid(pl0p0,pl0p01,pl0p05,pl0p1,pl0p2,ncol=1)

ggsave("~/projects/2DSFS-scan/Figures/FigPMeanFst_noDel.pdf",width=8.5*1.3,height=11*1.3)

# figure of change in genetic variance over time


data_mig0.0<-cbind(data_mig0.0,data.frame(mig=rep(0.0,length(data_mig0.0$cycle))))
data_mig0.01<-cbind(data_mig0.01,data.frame(mig=rep(0.01,length(data_mig0.01$cycle))))
data_mig0.05<-cbind(data_mig0.05,data.frame(mig=rep(0.05,length(data_mig0.05$cycle))))
data_mig0.1<-cbind(data_mig0.1,data.frame(mig=rep(0.1,length(data_mig0.1$cycle))))
data_mig0.2<-cbind(data_mig0.2,data.frame(mig=rep(0.2,length(data_mig0.2$cycle))))

dataAll<-rbind(data_mig0.0,data_mig0.01,data_mig0.05,data_mig0.1,data_mig0.2)

plMigSDNoDel<-ggplot(dataAll,aes(cycle-3000,SD,color=as.factor(mig)))+
  stat_summary(geom="line", fun="mean",size=1)+theme_classic()+ylab(expression(sigma))+xlab("generations after shift")+theme(text = element_text(size=14))+
  scale_color_manual(values=wes_palette("AsteroidCity1"),name="migration rate")


plot_grid(plMigSDWithDel,plMigSDNoDel,labels=c("A","B"),rel_widths = c(1,1.4))

ggsave("~/projects/2DSFS-scan/Figures/FigSigmaSupp.pdf", width=8.5,height=2.5)
