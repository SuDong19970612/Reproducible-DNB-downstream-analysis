#确定颜色填充node
stage<-read.csv("B_P_stage.csv",header = T)
View(stage)
#从stage挑选出颜色属性
genes<-read.table("DHX9node.txt",header = T)
View(genes)
genes<-genes$name
all_genes<-stage$hgnc_symbol
index<-which(all_genes%in%genes)
stage_DHX9node<-stage[index,]
View(stage_DHX9node)
write.csv(stage_DHX9node,"nodedis.csv")

node<-read.csv("nodedis.csv",header = T,row.names = 1)
View(node)
node<-as.matrix(node)

#标准化
node<-scale(node)
View(node)
write.csv(node,"nodescale.csv")
