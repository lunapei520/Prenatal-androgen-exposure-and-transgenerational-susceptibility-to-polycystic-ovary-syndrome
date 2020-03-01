test.df <- read.table("~/Desktop/oocytes_analyses/190226/counts_597DEGcvcdhv.txt",header = T,
                      row.names = NULL,sep=",")
test.df <- test.df[!duplicated(data.frame(test.df[,1])),]
rownames(test.df) <- test.df[,1]
test.df <- test.df[,2:ncol(test.df)]
metadata <- read.table("~/Desktop/oocytes_analyses/190225/sampleinfoF1F2F3_CVHVCDnew.txt",header = T,
                       row.names = NULL,sep="\t")
pca_data=prcomp(t(test.df))
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
metadata$celltype <- factor(metadata$celltype, c('CV','CD','HV'))
df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(test.df), condition=metadata$celltype, generation=metadata$tech)
head(df_pca_data)
library(ggplot2)
ggplot(df_pca_data, aes(PC1,PC2, color = condition, shape= generation))+
  geom_point(size=1.8,alpha = 1)+
  labs(x=paste0("PC1 (",pca_data_perc[1],"%",")"), y=paste0("PC2 (",pca_data_perc[2],"%",")")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c('#B4AF46','#B4464B','#56B4E9'),
                     guide=T) +scale_shape_manual(values = c(1, 2, 4))+
  guides(color = guide_legend(order = 1),
         size = guide_legend(order = 2),
         shape = guide_legend(order = 3))




