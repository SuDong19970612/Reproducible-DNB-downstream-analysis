library(GSVA)
#RNAseq的数据

#如果是RNA-seq的原始整数的read count 在使用gsva时需要设置:
#kcdf="Possion",
#如果是取过log的RPKM,TPM等结果可以使用默认的值，
#所以如果我们在使用fpkm进行分析的时候使用默认参数局可以了
data(commonPickrellHuang)
stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset),
                    featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))
#创建gmt文件
canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
data(genderGenesEntrez)
MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),collectionType=BroadCollection(category="c2"), setName="MSY")

XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),collectionType=BroadCollection(category="c2"), setName="XiE")

canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))
#使用GSVA方法进行计算

esmicro <- gsva(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500,mx.diff=TRUE, 
                verbose=FALSE, parallel.sz=1)

esrnaseq <- gsva(pickrellCountsArgonneCQNcommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500,
                 kcdf="Poisson", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
#pheatmap::pheatmap(esrnaseq)
