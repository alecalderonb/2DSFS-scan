"
plot 1D SFS on empirical data
"

setwd("~/Desktop/ECB/2DSFS_scan/emp_data/")

# combine files in long form

library(dyplr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(wesanderson)

files <- tibble(
  file = c("data/chr1.uv_sfs.txt",
           "data/chr1.bv_sfs.txt",
           "data/chrZ_4000000-4500000.uv_sfs.txt",
           "data/chrZ_4000000-4500000.bv_sfs.txt"),
  region = c("bg", "bg", "fg", "fg"),
  pop = c("uv", "bv", "uv", "bv")
)

sfs <- files %>%
  mutate(data = lapply(file, read.delim)) %>%
  select(-file) %>% 
  unnest(data)

plA <- ggplot(subset(sfs, pop == "uv"), aes(freq, density, fill=region))+
  theme_classic()+
  geom_bar(stat="identity",position="dodge")+
  scale_x_continuous(
    limits = c(0.5, 10.5),
    breaks = 1:10,
    labels = 1:10
  ) +
  scale_y_continuous(limits = c(0, 0.25)) +
  xlab("allele count")+ylab("density")+
  scale_fill_manual(values=wes_palette("Moonrise3"))+
  theme(legend.position = "NA")+theme(text = element_text(size=14))
plA

plB <- ggplot(subset(sfs, pop == "bv"), aes(freq, density, fill=region))+
  theme_classic()+
  geom_bar(stat="identity",position="dodge")+
  scale_x_continuous(
    limits = c(0.5, 10.5),
    breaks = 1:10,
    labels = 1:10
  ) +
  scale_y_continuous(limits = c(0, 0.25)) +
  xlab("allele count")+ylab("density")+
  scale_fill_manual(values=wes_palette("Moonrise3"))+
  #theme(text = element_text(size=14))+
  theme(legend.position = "NA")+theme(text = element_text(size=14))
plB

plot_grid(plA,plB,labels=c("A","B"),rel_widths=c(1,1))