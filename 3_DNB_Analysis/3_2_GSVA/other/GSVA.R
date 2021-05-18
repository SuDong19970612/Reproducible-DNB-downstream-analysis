# 学习GSVA
# https://www.cnblogs.com/raisok/p/11039239.html

#GSVA本质原理在于：研究感兴趣的基因集在不同样品间的差异，或者寻找比较重要的基因集，
#作为一种分析方法，主要是是为了从生物信息学的角度去解释导致表型差异的原因。

library(GSVA)

#模拟数据

p<-20000 #基因数
n<-30#样本数
nGS<-100 #基因集数目
#最大与最小基因集大小
min.sz<-10 
max.sz<-100

X <- matrix(rnorm(p*n), nrow=p, dimnames=list(1:p, 1:n))
View(X)
dim(X)
#样本基因集大小（size）
gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE))
#样本基因集
gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p)

es.max <- gsva(X, gs, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
es.dif <- gsva(X, gs, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
pheatmap::pheatmap(es.max)
pheatmap::pheatmap(es.dif)


