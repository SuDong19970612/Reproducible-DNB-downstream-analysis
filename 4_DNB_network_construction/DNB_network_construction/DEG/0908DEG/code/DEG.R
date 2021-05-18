rm(list = ls())
setwd("G:/biology/CRC/result/CNP/exct_DNB/c1/")
Stage<-read.csv(file = "G:/biology/CRC/result/CNP/c1/count/CNP_Stage/P/Pstage.csv",header = T,row.names = 1)
#使用tidyverse去除版本号
options(stringsAsFactors = F)
my_ensembl_gene_id<-row.names(Stage)
a<-data.frame(my_ensembl_gene_id)
library(org.Hs.eg.db)
library(tidyverse)
g2s <- toTable(org.Hs.egSYMBOL)
g2e <- toTable(org.Hs.egENSEMBL)
a1<-a
a1$V1 = apply(a1[1], 1,function(x){
  str_split(x,'[.]')[[1]][1]
}) %>% unlist()
my_ensembl_gene_id<-a1$V1


#差异基因注释，用bioMart对差异表达基因进行注释
library(biomaRt)
library(curl)
#使用人的注释数据
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mms_symbols<- getBM(attributes=c('ensembl_gene_id','chromosome_name',"hgnc_symbol","hgnc_id"),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
#合并
Stage<-cbind(a1$V1,Stage)
names(Stage)[1]<-"ensembl_gene_id"
Stage<-merge(Stage,mms_symbols,by="ensembl_gene_id")
Stage<-Stage[,c(2,3,4,5,6,7,8,9,11)]
write.csv(Stage,"Stage.csv")

#提取
Stage<-read.csv("Stage.csv",header = T)
DNB<-read.table("DNB.txt",header = T)
DNB<-DNB$Symbol
Gene<-Stage$hgnc_symbol
index<-which(Gene%in%DNB)
heat_Stage<-Stage[index,]
heat_Stage<-heat_Stage[,-1]
write.csv(heat_Stage,file = "heat_Stage.csv")
#过滤
heat_Stage<-read.csv("heat_Stage.csv",header = T,row.names = 10)
heat_Stage<-heat_Stage[,-1]
index<-which(apply(heat_Stage,1,sum)!=0)
heat_Stage<-heat_Stage[index,]
write.csv(heat_Stage,file = "heat_stage.csv")
#分期
stage1<-apply(heat_Stage[,1:2],1,mean)
stage2<-heat_Stage$LC01
stage3<-apply(heat_Stage[,4:6],1,mean)
stage4<-apply(heat_Stage[,7:8],1,mean)
stage<-round(cbind(stage1,stage2,stage3,stage4))
names(stage)<-c("T1","T2","T3","T4")
card<-row.names(stage)
write.csv(stage,"DNB_stage.csv")
write(card,"card.txt")





#绘制热图
