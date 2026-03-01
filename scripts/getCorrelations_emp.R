"
January 2026
calculate correlation between T2D and FST on empirical data
"

library(dplyr)

read.table("Desktop/ECB/2DSFS_scan/emp_data/windows.1000.txt")->data
zchr <- c()
for (i in seq(1,length(data$V2))) {
  z = "Autosome/W"
  if (data$V1[i] == "NC_087119.1") {
    z = "Z"
  }
  zchr<-c(zchr,z)
}

data<-cbind(data,data.frame(zchr = zchr))

# Autosome/W correlations
# V4=T2D / V5=FST
data %>%
  filter(zchr == "Autosome/W") %>%
  with(cor.test(V4, V5, method = "pearson"))

data %>%
  filter(zchr == "Z") %>%
  with(cor.test(V4, V5, method = "pearson"))
