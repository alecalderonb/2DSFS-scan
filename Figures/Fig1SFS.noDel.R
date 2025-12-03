library(ggplot2)
library(cowplot)
library(wesanderson)

read.table("~/projects/2DSFS-scan/data/sfs.sim.noDel.txt")->sfs

plC<-ggplot(sfs[sfs$V2 == 3000,],aes(V3,V4,fill=V1))+theme_classic()+geom_bar(stat="identity",position="dodge")+xlim(c(0.5,5.5))+xlab("allele count")+ylab("Density")+scale_fill_manual(values=wes_palette("Moonrise3"))+theme(legend.position = "NA")+theme(text = element_text(size=14))
plD<-ggplot(sfs[sfs$V2 == 3300,],aes(V3,V4,fill=V1))+theme_classic()+geom_bar(stat="identity",position="dodge")+xlim(c(0.5,5.5))+xlab("allele count")+ylab("Density")+scale_fill_manual(values=wes_palette("Moonrise3"),name="",labels=c("background","foreground"))+theme(text = element_text(size=14))

read_csvs<-function(iter) {
    connection <- gzfile(paste("~/projects/2DSFS-scan/slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams_noDel/migRate.0.01/iter_",iter,"/sim_log.",iter,".txt.gz",sep=""), "rt") # "rt" for read text mode
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

ggsave("~/projects/2DSFS-scan/Figures/Fig1.noDel.pdf",width=10,height=7)

