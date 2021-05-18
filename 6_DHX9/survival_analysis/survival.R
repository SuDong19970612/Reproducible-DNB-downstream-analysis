rm(list = ls())
options(stringsAsFactors = F)
# DHX9 survival analysis
library(survival)
library("survminer")
#data import
CRC<-read.csv("survival.csv",header = T,row.names = 1)
View(CRC)
