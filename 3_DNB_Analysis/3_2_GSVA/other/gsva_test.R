library(GSEABase)
library(GSVAdata)
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
#官方文档有数据的预处理过程
cacheDir <- system.file("extdata", package="GSVA")
cachePrefix <- "cache4vignette_"
file.remove(paste(cacheDir, list.files(cacheDir, pattern=cachePrefix), sep="/"))
data(leukemia)
leukemia_eset
head(pData(leukemia_eset))
table(leukemia_eset$subtype)
#过滤
data(leukemia)
leukemia_eset
#package="genefilter"
filtered_eset <- nsFilter(leukemia_eset, require.entrez=TRUE, 
                          remove.dupEntrez=TRUE,var.func=IQR, 
                          var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE,
                          feature.exclude="^AFFX")
leukemia_filtered_eset <- filtered_eset$eset

cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,min.sz=10, max.sz=500,
                          verbose=TRUE),
      dir=cacheDir, prefix=cachePrefix)
#pheatmap::pheatmap(leukemia_es)
#取具有差异的基因集来可视化
#GSVA的计算差异基因集 
#limma包 对象为Gsva下的leukemia_es对象
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)

design <- model.matrix(~ factor(leukemia_es$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, design)
fit <- eBayes(fit)

allgeneSets<-topTable(fit,coef = "MLLvsALL",number = Inf)
DEgeneSets <- topTable(fit, coef="MLLvsALL", number=Inf,
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)

#提取差异基因集索引可视化
DEgeneSetsgenes<-row.names(DEgeneSets)
index<-which(row.names(leukemia_es@assayData$exprs)%in%DEgeneSetsgenes)
DEgeneSetsexp<-leukemia_es@assayData$exprs[index,]
pheatmap::pheatmap(DEgeneSetsexp)
# limma计算差异表达基因 对象为leukemia_filtered_eset
logFCcutoff <- log2(2)
design <- model.matrix(~ factor(leukemia_eset$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_filtered_eset, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgenes <- topTable(fit, coef="MLLvsALL", number=Inf,
                    p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)
res <- decideTests(fit, p.value=adjPvalueCutoff, lfc=logFCcutoff)
summary(res)
