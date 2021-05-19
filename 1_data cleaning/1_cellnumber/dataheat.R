cellnum<-read.table("1_cellnumber/numcell.txt")
cellnum<-as.matrix(cellnum)
annotation_col<-data.frame(Patient= c("RC03","R02","LC01","RC02","R03","C473","C481","LC03"),
                           Stage= c("I","I","II","III","III","III","IV","IV"),row.names = 1)
pheatmap::pheatmap(cellnum,cluster_rows = F,cluster_cols = F,legend_breaks =c(0,20,50,100,200),legend_labels = c("0","20","50","100","200") ,display_numbers = T,
                   number_color = "blue",number_format = "%.f",annotation_col = annotation_col,gaps_col = c(2,3,6))

