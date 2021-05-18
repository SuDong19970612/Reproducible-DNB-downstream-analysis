Mfuzz_exp<-read.csv("Mfuzz_standardize.csv",header = T,row.names = 1)
col_stage<-read.table("colo.txt")
col_stage<-col_stage$V1
index<-which(rownames(Mfuzz_exp)%in%col_stage)
DHX9_col_stage<-Mfuzz_exp[index,]
write.csv(DHX9_col_stage,"DHX9_col_stage.csv")
