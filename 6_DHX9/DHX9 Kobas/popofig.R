#DHX9kobas<-read.csv("data/DHX9 Kobas.csv",header = T)
#View(DHX9kobas)
#DHX9kobas$p<- round(-log10(DHX9kobas$Corrected.P.Value))
#popo<-data.frame(DHX9kobas$Fun,DHX9kobas$number,DHX9kobas$p)
#View(popo)
#write.csv(popo,"popo.csv")

popo<-read.csv("popo.csv",header = T)

library(ggplot2)
ggplot(data = popo)+
  geom_point(aes(x=number,y=Fun,size=p,color=number))+
  scale_color_gradient(low="green",high ="red")+
  labs(x="Genecounts",y="Pathway name",title="Pathway enrichment",size=expression(-log[10](Qvalue)))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)
      )

