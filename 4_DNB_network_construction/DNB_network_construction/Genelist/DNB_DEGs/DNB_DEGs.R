DNB<-read.csv("DNB_top50_genes_edge_list.csv")
DEG<-read.table("DNBDEGs.txt")
DEG<-DEG$V1
tar<-DNB$target
index<-which(tar%in%DEG)
DNB_DEG<-DNB[index,]
write.csv(DNB_DEG,"DNB_DEG.csv")
