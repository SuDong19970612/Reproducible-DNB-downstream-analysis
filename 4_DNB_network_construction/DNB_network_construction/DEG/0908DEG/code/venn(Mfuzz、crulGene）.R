#绘制Venn图
#install.packages("VennDiagram")
library(VennDiagram)
library(dplyr)
#载入数据
DG<-read.table("data/venn（mfuzz、crulGene）/c3c4.txt")
cru<-read.table("result/venn/CrulGene.txt")
a<-intersect(DG$V1,cru$V1)
(length(a))
write(a,"result/venn（mfuzz、cru）/venn_result.txt")
#绘制韦恩图
venn.diagram(list(NB=DG$V1,CR=cru$V1),height=2500,width=2500,resolution=500,
             fill=c("red","blue"),filename = "result/venn（mfuzz、cru）/venn.tiff")
