options(stringsAsFactors = F)
B_P_stage<-read.csv("B_P_stage.csv")
View(B_P_stage)
dim(B_P_stage)
B_P_DNB<-read.table("B_P_DNB.txt")
index<-which(B_P_stage$hgnc_symbol%in% B_P_DNB$V1)
B_P_stage_DNBexp<-B_P_stage[index,]
row.names(B_P_stage_DNBexp)<-B_P_stage_DNBexp$hgnc_symbol
B_P_stage_DNBexp<- B_P_stage_DNBexp[,-1]
View(B_P_stage_DNBexp)
dim(B_P_stage_DNBexp)
B_P_stage_DNBexp_scale<-scale(B_P_stage_DNBexp,center = T,scale = T)
View(B_P_stage_DNBexp_scale)

pheatmap::pheatmap(B_P_stage_DNBexp_scale,cluster_rows = F,cluster_cols = T)
