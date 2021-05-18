genes<-read.csv("Gene.csv",header = T)
genes$DEG<-30*genes$DEG
genes$Cluster<-30*genes$Cluster
genes$Immune<-30*genes$Immune
fix(genes)

View(genes)
library(reshape2)
genes <- melt(genes, id="Gene")
write.csv(genes,"genes_melt.csv")

gene<-read.csv("genes_melt.csv",header = T)



library(ggplot2)
set.seed(2020)
p <- ggplot(gene,aes(x=Gene,y=value,fill=variable))+
  geom_bar(position="dodge",stat="identity",width=0.8)+
  guides(fill=guide_legend(label = T,reverse=F,title="Standard"))+
  theme_classic()
p

