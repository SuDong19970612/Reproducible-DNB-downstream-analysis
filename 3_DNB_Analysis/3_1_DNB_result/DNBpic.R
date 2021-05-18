rm(list = ls())
DNB<-read.csv("data/DNB.csv",header = T,row.names = 1)
#View(DNB)
library(gridExtra)
library(ggplot2)
#par(mfrow=c(2,2)) 
set.seed(2020)
ggplot(data = DNB,mapping = aes(x=timepoint,y=Std,group = 1))+
  geom_line(size=2,color="red")+
  xlab("timepoint")+ylab("Std")+
  theme_classic()
ggplot(data = DNB,mapping = aes(x=timepoint,y=Inpcc,group = 1))+
  geom_line(size=2,color="red")+
  xlab("timepoint")+ylab("Inpcc")+
  theme_classic()
ggplot(data = DNB,mapping = aes(x=timepoint,y=Outpcc,group = 1))+
  geom_line(size=2,color="red")+
  xlab("timepoint")+ylab("Outpcc")+
  theme_classic()
ggplot(data = DNB,mapping = aes(x=timepoint,y=ComplexIndex,group = 1))+
  geom_line(size=2,color="red")+
  xlab("timepoint")+ylab("CI")+
  theme_classic()
#grid.arrange(p1,p2,p3,p4,nrow=2)
