rm(list = ls())
options(stringsAsFactors = F)

#need package GSVA
library(GSVA)
#用于读取gmt文件格式
library(GSEABase)
#用于基于通路ES（富集分数）矩阵的差异分析
#当ES值大于0时，表示某一功能基因富集在排序序列的前端，若为小于0时，则某一功能基因富集在排序序列的后端，
#ES值越高说明这些基因在通路中有富集，非散在分布
library(limma)
#读取基因集文件
geneSets <- getGmt("Data/immune.gmt")
#geneSets <- getGmt("Data/kegg.v7.1.symbols.gmt")
#读取表达量文件并去除重复(我的表达矩阵Gene没有重复，走个过程)
mydata <- read.csv(file = "Data/B_P_DNB_genes_stage_1.csv",header = T)
name=mydata[,1]
index <- duplicated(mydata[,1])
fildup=mydata[!index,]
exp=fildup[,-1]
row.names(exp)=name
#将数据框转换成矩阵
mydata= as.matrix(exp)
mydata=round(mydata)
#使用gsva方法进行分析，默认mx.diff=TRUE，min.sz=1,max.zs=Inf，这里是设置最小值和最大值
#kcdf=Poisson 使用的为count数据
Es <- gsva(mydata, geneSets, parallel.sz=1,kcdf="Poisson")
pheatmap::pheatmap(Es,angle_col = 0,cluster_cols = F,cluster_rows = T,color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
write.csv(Es,"GSVA_ES.csv")

#limma 差异基因集分析
adjPvalueCutoff <- 0.05
logFCcutoff <- log2(2)
#设置分组
class<-c("treat","treat","control","control")
design<-model.matrix(~factor(class))
colnames(design)<-c("control","treat")

fit <- lmFit(Es, design)
fit <- eBayes(fit)

allgeneSets<-topTable(fit,coef = 2,adjust="fdr",number = 200000)
DEgeneSets <- topTable(fit, coef=2, number=Inf,p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)


#提取差异基因集索引可视化
#DEgeneSetsgenes<-row.names(DEgeneSets)
#index<-which(row.names(leukemia_es@assayData$exprs)%in%DEgeneSetsgenes)
#DEgeneSetsexp<-leukemia_es@assayData$exprs[index,]
#pheatmap::pheatmap(DEgeneSetsexp)


mydata<-read.csv("Data/Bcell_GSVA.csv",row.names = 1)
#name=mydata[,1]
#index <- duplicated(mydata[,1])
#fildup=mydata[!index,]
#mydata<-fildup
#write.csv(mydata,"Data/Bcell_GSVA.csv")
pheatmap::pheatmap(mydata[20:49,],cluster_cols = F,angle_col = 0)
