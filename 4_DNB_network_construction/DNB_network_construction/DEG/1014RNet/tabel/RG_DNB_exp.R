# 读取表达矩阵
exp<-read.csv("B_P_stage_dump.csv",row.names = 1)
View(exp)
Genes<-row.names(exp)
class(Genes)
#RG_DNB
RG.DNB_New<-read.csv("net/RG_DNB_edge.txt default node.csv")
RG_DNB<-RG.DNB_New$name
#提取RG_DNB node 表达矩阵
index<-which(Genes%in%RG_DNB)
RG_DNB_exp<-exp[index,]
dim(RG_DNB_exp)
#标准化
RG_DNB_exp_scale<-scale(RG_DNB_exp,center=T,scale=T)
write.csv(RG_DNB_exp,"RG_DNB_exp.csv")
write.csv(RG_DNB_exp_scale,"RG_DNB_exp_scale.csv")
#log2
RG_DNB_lg<-log2(RG_DNB_exp)
write.csv(RG_DNB_lg,"RG_DNB_lg.csv")
RG_DNB_lg_scale<-scale(RG_DNB_lg,center = T,scale = T)
write.csv(RG_DNB_lg_scale,"RG_DNB_lg_scale.csv")

#吐了，发现这里的标准化可以用到之前的一个结果
# 读取表达矩阵
exp<-read.csv("exp_standardise.csv",row.names = 1)
View(exp)
Genes<-row.names(exp)
class(Genes)
#RG_DNB
RG.DNB_New<-read.csv("net/RG_DNB_edge.txt default node.csv")
RG_DNB<-RG.DNB_New$name
#提取RG_DNB node 表达矩阵
index<-which(Genes%in%RG_DNB)
RG_DNB_exp<-exp[index,]
dim(RG_DNB_exp)
write.csv(RG_DNB_exp,"RG_DNB_exp_anotation.csv")
