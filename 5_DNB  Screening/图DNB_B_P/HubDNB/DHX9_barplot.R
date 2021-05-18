genes<-read.csv("Gene.csv",header = T)

fix(genes)

View(genes)
library(reshape2)
genes <- melt(genes, id="Gene")
write.csv(genes,"genes_melt.csv")

genes<-read.csv("genes_melt.csv",header = T)

library(ggplot2)
set.seed(2020)
p <- ggplot(genes,aes(x=Gene,y=value,fill=variable))+
  geom_bar(position="dodge",stat="identity",width=0.8,colour="black")+
  guides(fill=guide_legend(reverse=F,title="Standard"))

p+theme_classic()
