对筛选出来的翻转基因（RG），即stage1 vs stage2 无差异，stage2 vs stage3 差异显著的基因；以及DNB基因，在Cytoscape（version:3.7.1）中按照表达量可视化。

首先，提取RG 与 DNB 基因的网络关系对。

```R
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
```

Cytoscape(version:3.7.1)可视化RG_DNB关系对：

<img src="https://user-images.githubusercontent.com/38640955/95998801-cf767780-0e67-11eb-8750-e204abbc1999.png" alt="RG_DNB_edge" style="zoom:70%;" />

第二步，提取网络注释信息，先提取RG_DNB基因表达矩阵。

```R
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
```
DNB基因在DNB_RG网络中的核心地位，如下图：

![DNB_net](https://user-images.githubusercontent.com/38640955/96229634-cc939800-0fc9-11eb-9b98-ab6af17ba8f5.png)

Cytoscope（version:3.7.1）中，使用基因表达注释矩阵“RG_DNB_exp.csv”，可获得时序性RG_DNB翻转网络图：

![RG_DNB](https://user-images.githubusercontent.com/38640955/96235275-29dd1880-0fcd-11eb-976b-ab52e6c95da0.png)



