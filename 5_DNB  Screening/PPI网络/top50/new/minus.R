node<-read.csv("nodedis.csv",header = T,row.names = 1)
View(node)
node<-scale(node,center = T,scale = T)
write.csv(node,"nodeminus.csv")
