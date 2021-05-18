#绘制Venn图
#install.packages("VennDiagram")
library(VennDiagram)
library(dplyr)
#载入数据
NoDEG<-read.csv("result/anno/NDEG(1vs2)_anno.csv")
DEG<-read.csv("result/anno/DEG（2vs3)_anno.csv")
intersec<-intersect(NoDEG$hgnc_symbol,DEG$hgnc_symbol)
(length(intersec))
write(intersec,"result/venn/CrulGene.txt")
#绘制韦恩图
venn.diagram(list(NODEG=NoDEG$hgnc_symbol,DEG=DEG$hgnc_symbol),fill=c("red","blue"),filename = "result/venn/CrulGene.tiff")
