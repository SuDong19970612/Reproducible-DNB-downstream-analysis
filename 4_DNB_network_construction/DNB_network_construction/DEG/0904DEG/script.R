#require library
library(DESeq2)
#import data
exp <- read.csv("B_P_stage.csv",row.names = NULL)
#去重
which(duplicated(exp[,1]))
#6028
mean_POL<-apply(exp[c(6027,6028),-1],2,mean)
exp2<-exp[!duplicated(exp[,1]),]
row.names(exp2)<-exp2[,1]
exp2<-exp2[,-1]
write.csv(exp2,"B_P_stage_dump.csv")
exp_dump<-read.csv("B_P_stage_dump.csv",row.names = 1)
exp_dump<-round(as.matrix(exp_dump))
write.csv(exp_dump,"exp_dum2.csv")
condition <- factor(c(rep("control",2),rep("treat",2)), levels = c("control","treat"))
coldata<-data.frame(row.names=colnames(exp_dump), condition)
#差异分析
dds <- DESeqDataSetFromMatrix(exp_dump, coldata, design= ~ condition)
dds<-DESeq(dds)
res = results(dds, contrast=c("condition", "control", "treat"))
res = res[order(res$pvalue),]
write.csv(res,file="result.csv")
diff_gene <-subset(res, padj < 0.05 & abs(log2FoldChange) > 2)
write.csv(diff_gene,file = "diff_gene.csv")



