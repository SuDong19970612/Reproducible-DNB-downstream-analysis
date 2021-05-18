#读net数据
DNB_top50_genes_edge_list <- read_csv("DNB_top50_genes_edge_list.csv")
DNB_edge<-DNB_top50_genes_edge_list[,1:2]
#去自成环
s<-DNB_edge$source
t<-DNB_edge$target
index = which(s == t) 
DNB_edge_noLoop<-DNB_edge[-index,]  
View(DNB_edge_noLoop)  
write.csv(DNB_edge_noLoop,"DNB_edge_noLoop.csv")
# 读取RG基因
RG<-read.table("reverseGene.txt")
# 129 个有翻转模式的基因
(length(RG$V1))
# 从DNB_edge_noLoop中提取关系对
index<-which(DNB_edge_noLoop$target%in%RG$V1)
# 共有310个关系对
(length(index))
#存出RG_DNB关系对
RG_DNB_edge<-DNB_edge_noLoop[index,]
write.csv(RG_DNB_edge,"RG_DNB_edge.csv")

#RG.DNB,手动合并的DNB与RG基因 179
RG.DNB<-read.table("RG.DNB.txt")
dim(RG.DNB)
#真实存在的RG_DNB基因 168 
RG.DNB_New<-read.csv("net/RG_DNB_edge.txt default node.csv")
dim(RG.DNB_New)

