library("ggplot2")
library("dplyr")

setwd("/home/daniil/programming/cpp/SPBUSignificanceEstimation")
df <- read.table(file="temperatureChoice.txt", header=FALSE, sep=" ", col.names = c("T", "score"))
head(df)
ggplot(df, aes(x=T, y=score)) + geom_point(size=0.5) + scale_y_continuous(limits=c(0,10^-13)) +
  geom_smooth(formula= y~x)
  

