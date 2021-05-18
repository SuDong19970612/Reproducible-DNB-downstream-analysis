options(stringsAsFactors = F)
library(pheatmap)
library(Cairo)
rt<-read.csv("data/pheatmap .csv",row.names = 1)
rt<-t(rt)
pheatmap(rt,cluster_cols = F,cluster_rows = F,cellwidth = 8, cellheight = 8, fontsize = 10,
         legend = FALSE,
         angle_col = "90",
         color = c("white","brown"),
         gaps_row = 1:15,
         gaps_col = 1:24)
