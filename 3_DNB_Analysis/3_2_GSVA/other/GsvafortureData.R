rm(list = ls())
options(stringsAsFactors = F)

#need package GSVA
library(GSVA)
#读取基因集文件
geneSets <- getGmt("test.geneset")
#读取表达量文件并去除重复
mydata <- read.table(file = "all.genes.fpkm.xls",header=T)
name=mydata[,1]
index <- duplicated(mydata[,1])
fildup=mydata[!index,]
exp=fildup[,-1]
row.names(exp)=name
#将数据框转换成矩阵
mydata= as.matrix(exp)

#使用gsva方法进行分析，默认mx.diff=TRUE，min.sz=1,max.zs=Inf，这里是设置最小值和最大值
res_es <- gsva(mydata, geneSets, min.sz=10, max.sz=500, verbose=FALSE, parallel.sz=1)
pheatmap(res_es)

#mx.diff=FALSE es值是一个双峰的分布
es.max <- gsva(mydata, geneSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
pheatmap(es.max)

#mx.diff=TURE es值是一个近似正态分布
es.dif <- gsva(mydata, geneSets, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
pheatmap(es.dif)

#可以看一下两种不同分布的效果,前者是高斯分布，后者是二项分布
par(mfrow=c(1,2), mar=c(4, 4, 4, 1))
plot(density(as.vector(es.max)), main="Maximum deviation from zero",xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
plot(density(as.vector(es.dif)), main="Difference between largest\npositive and negative deviations",xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)



#设置分组
col = names(exp)
sample=col[1:10]
group=c(rep('control',6),rep('treat',4))
phno=data.frame(sample,group)

Group<-factor(phno$group,levels=levels(phno$group))
design<-model.matrix(~0+Group)
colnames(design) <- c("control", "treat")

#获取需要进行差异分析的分组
res=es.max[,1:10]
#定义阈值
logFCcutoff <- log2(1.5)
adjPvalueCutoff <- 0.001
#进行差异分析
fit <- lmFit(res, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="treat", number=Inf)
DEgeneSets <- topTable(fit, coef="treat", number=Inf,p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
#summary(res)

#画火山图
DEgeneSets$significant="no"
DEgeneSets$significant=ifelse(DEgeneSets$logFC>0|DEgeneSets$logFC<0,"up","down")
ggplot(DEgeneSets,aes(logFC,-1*log10(adj.P.Val)))+geom_point(aes(color =significant)) + xlim(-4,4) + ylim(0,30)+labs(title="Volcanoplot",x="log[2](FC)", y="-log[10](FDR)")+scale_color_manual(values =c("#00ba38","#619cff","#f8766d"))+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
#获取差异基因集的表达量
DEgeneSetspkm = merge(DEgeneSets,es.max,by=0,all.x=TRUE)[,c(1,11:20)]
degsetsp=DEgeneSetspkm[,-1]
name=DEgeneSetspkm[,1]
row.names(degsetsp)=name
pheatmap(degsetsp)

