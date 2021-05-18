

#edgeR
library(edgeR)
Pstage<-read.csv(file = "stage2_B_PC_exp.csv",header = T,row.names = 1)
Pstage<-round(Pstage)
counts<-Pstage
group<-1:2
y <- DGEList(counts=counts, group = group)
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y_bcv <- y
bcv <- 0.4
et <- exactTest(y_bcv, dispersion = bcv ^ 2)
gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)
write.csv(et$table,"gene_exp.csv")
#FDR矫正
exp<-read.csv("gene_exp.csv",header = T,row.names = 1)
p<-exp$PValue

#install.packages("fdrtool")
library(fdrtool)
fdr=fdrtool(p,statistic="pvalue")
qval<-fdr$qval
exp<-cbind(exp,qval)
write.csv(exp,"exp_qvalue.csv")

#差异基因log2FC 大于2，padj<0.05
exp<-read.csv("exp_qvalue.csv",header = T,row.names = 1)
diffgene <-subset(exp, qval < 0.05 & abs(logFC) > 2)
dim(diffgene)
head(diffgene)
write.csv(diffgene,file = "diffgene.csv")
#注释差异基因
#exp<-read.csv("diffgene.csv",header = T,row.names = 1)
#ESM注释前需要去除ESM版本号，使用bioMart注释
options(stringsAsFactors = F)
my_ensembl_gene_id<-row.names(exp)
a<-data.frame(my_ensembl_gene_id)
library(org.Hs.eg.db)
library(tidyverse)
g2s <- toTable(org.Hs.egSYMBOL)
g2e <- toTable(org.Hs.egENSEMBL)
a2<-a
a2$V2 = apply(a2[1], 1,function(x){
  str_split(x,'[.]')[[1]][1]
}) %>% unlist()
my_ensembl_gene_id<-a2$V2
#bioMart注释
library(biomaRt)
library(curl)
#使用人的注释数据
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mms_symbols<- getBM(attributes=c('ensembl_gene_id','chromosome_name',"hgnc_symbol","hgnc_id"),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
#合并
exp<-cbind(a2$V2,exp)
names(exp)[1]<-"ensembl_gene_id"
exp<-merge(exp,mms_symbols,by="ensembl_gene_id")
write.csv(exp,"exp_annotation.csv")

#绘图
exp<-read.csv("exp_qvalue.csv",header = T,row.names = 1)
pdf(file = "vol.pdf")
plot(-log10(exp$qval),exp$logFC,xlab="adj.P.Val",ylab="logFC",main="Volcano",pch=20,cex=0.4)
diffSub=subset(exp,exp$qval<0.05 & abs(exp$logFC)>2 )
points(-log10(diffSub$qval),diffSub$logFC,pch=20,col="red",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()



