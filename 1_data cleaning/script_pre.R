options(digits = 2,scipen = 100,stringsAsFactors = F)


all.ReadCount.gene.matrix<-read.table("All.ReadCount.gene.matrix.txt",header = T, sep='\t')
names(all.ReadCount.gene.matrix)
#make index 
phenoMatrix<-read.csv("phenoMatrix.csv",header = T)
B_pheno<-phenoMatrix[phenoMatrix$Cell_Type==2,]
index_id<-B_pheno$ID
index<-which( names(all.ReadCount.gene.matrix)%i
# segregated B cell from C n%index_id )
# segregated B cell
B.Readcount.gene.matrix<-all.ReadCount.gene.matrix[,c(1,index)]
dim(B.Readcount.gene.matrix)
write.csv(B.Readcount.gene.matrix,"B.Readcount.gene.matrix.csv")
write.csv(B_pheno,"B_pheno.csv")

index_tissue<-which(names(B.Readcount.gene.matrix)%in%grep("_C",index_id,value = T))  
B_C_Matrix<-B.Readcount.gene.matrix[,c(1,index_tissue)]
write.csv(B_C_Matrix,"B_C_Matrix.csv")

# segregated B cell from P
index_tissue<-which(names(B.Readcount.gene.matrix)%in%grep("_P",index_id,value = T))  
B_P_Matrix<-B.Readcount.gene.matrix[,c(1,index_tissue)]
write.csv(B_P_Matrix,"B_P_Matrix.csv")


# By patient ID
# C
names<-names(B_C_Matrix)
index_C473<-grep("C473",names)
index_C481<-grep("C481",names)
index_LC03<-grep("LC03",names)
index_RC02<-grep("RC02",names)
index_RC03<-grep("RC03",names)
index_LC01<-grep("LC01",names)
index_R03<-grep("R03",names)
index_R02<-grep("R02",names)
# aver
C473_Exp_ave<-B_C_Matrix[,index_C473]
C481_Exp_ave<-apply(B_C_Matrix[,index_C481],1,mean)
LC03_Exp_ave<-apply(B_C_Matrix[,index_LC03],1,mean)
RC02_Exp_ave<-apply(B_C_Matrix[,index_RC02],1,mean)
RC03_Exp_ave<-apply(B_C_Matrix[,index_RC03],1,mean)
LC01_Exp_ave<-apply(B_C_Matrix[,index_LC01],1,mean)
R03_Exp_ave<-apply(B_C_Matrix[,index_R03],1,mean)
R02_Exp_ave<-apply(B_C_Matrix[,index_R02],1,mean)

Patient_B_C<-data.frame(B_C_Matrix$Ensembl_ID,C473_Exp_ave,C481_Exp_ave,
                        LC03_Exp_ave,RC02_Exp_ave,RC03_Exp_ave,
                        LC01_Exp_ave,R03_Exp_ave,R02_Exp_ave)
write.csv(Patient_B_C,"Patient_B_C.csv")

# P
names<-names(B_P_Matrix)
index_C473<-grep("C473",names)
index_C481<-grep("C481",names)
index_LC03<-grep("LC03",names)
index_RC02<-grep("RC02",names)
index_RC03<-grep("RC03",names)
index_LC01<-grep("LC01",names)
index_R03<-grep("R03",names)
index_R02<-grep("R02",names)
# aver
C473_Exp_ave<-apply(B_P_Matrix[,index_C473],1,mean)
C481_Exp_ave<-apply(B_P_Matrix[,index_C481],1,mean)
LC03_Exp_ave<-apply(B_P_Matrix[,index_LC03],1,mean)
RC02_Exp_ave<-apply(B_P_Matrix[,index_RC02],1,mean)
RC03_Exp_ave<-apply(B_P_Matrix[,index_RC03],1,mean)
LC01_Exp_ave<-apply(B_P_Matrix[,index_LC01],1,mean)
R03_Exp_ave<-apply(B_P_Matrix[,index_R03],1,mean)
R02_Exp_ave<-apply(B_P_Matrix[,index_R02],1,mean)

Patient_B_P<-data.frame(B_P_Matrix$Ensembl_ID,C473_Exp_ave,C481_Exp_ave,
                        LC03_Exp_ave,RC02_Exp_ave,RC03_Exp_ave,
                        LC01_Exp_ave,R03_Exp_ave,R02_Exp_ave)
write.csv(Patient_B_P,"Patient_B_P.csv")

#修改Excel名字
#by stage

# C
Patient_B_C<-read.csv("result/Patient_B_C.csv",header = T,row.names = 1)
stage1<-apply(Patient_B_C[,c("R02","RC03")],1,mean)
stage2<-Patient_B_C[,c("LC01")]
stage3<-apply(Patient_B_C[,c("R03","RC02","C473")],1,mean)
stage4<-apply(Patient_B_C[,c("C481","LC03")],1,mean)
Bystage_Matrix<-cbind(stage1,stage2,stage3,stage4)
write.csv(Bystage_Matrix,"result/B_C_Stage_Matrix.csv")
# P
Patient_B_P<-read.csv("result/Patient_B_P.csv",header = T,row.names = 1)
stage1<-apply(Patient_B_P[,c("R02","RC03")],1,mean)
stage2<-Patient_B_P[,c("LC01")]
stage3<-apply(Patient_B_P[,c("R03","RC02","C473")],1,mean)
stage4<-apply(Patient_B_P[,c("C481","LC03")],1,mean)
Bystage_Matrix<-cbind(stage1,stage2,stage3,stage4)
write.csv(Bystage_Matrix,"result/B_P_Stage_Matrix.csv")


#minus zero
#B_C
patient_no_0 <- patient_df[!(patient_df$C481==0|patient_df$LC03==0|patient_df$LC01==0
                             |patient_df$RC02==0|patient_df$R02==0|patient_df$R03==0|patient_df$RC03==0),]

stage_no_0 <- stage_df[!(stage_I==0|stage_II==0|stage_III==0|stage_IV==0),]

write.csv(stage_no_0,file = "stage_no_0.csv")
write.csv(patient_no_0,file = "patient_no_0.csv")

#B_P
