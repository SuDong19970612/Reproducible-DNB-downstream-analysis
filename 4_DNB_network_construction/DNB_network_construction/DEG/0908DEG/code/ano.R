

#去版本号
options(stringsAsFactors = F)
data1<-read.table("data/B_P_stage/DEG（2vs3）.txt")
my_ensembl_gene_id<-data1$V1

a<-data.frame(my_ensembl_gene_id)

library(org.Hs.eg.db)
library(tidyverse)

g2s <- toTable(org.Hs.egSYMBOL)
g2e <- toTable(org.Hs.egENSEMBL)
a1<-a
a1$V1 = apply(a1[1], 1,function(x){
  str_split(x,'[.]')[[1]][1]
}) %>% unlist()
my_ensembl_gene_id<-a1$V1


#差异基因注释，用bioMart对差异表达基因进行注释
library(biomaRt)
library(curl)
#使用人的注释数据
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mms_symbols<- getBM(attributes=c('ensembl_gene_id','chromosome_name',"hgnc_symbol","hgnc_id"),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
write.csv(mms_symbols,"result/anno/DEG（2vs3)_anno.csv")
