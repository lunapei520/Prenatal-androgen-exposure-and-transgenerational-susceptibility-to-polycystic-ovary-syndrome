x <- read.table("~/Desktop/oocytes_analyses/counts_CVCDnew.txt",header = T, row.names =NULL,
                sep="\t")
x <- x[!duplicated(data.frame(x[,1])),]
rownames(x) <- x[,1]
x <- x[,2:ncol(x)]
library(dplyr)
commongene <- read.table("~/Desktop/oocytes_analyses/top20CVCD.txt",header = F, row.names = NULL,
                         sep="\t")
data <- x[which(rownames(x) %in% commongene), ]
metadata <- read.table(file = "~/Desktop/oocytes_analyses/sampleinfoF1F2F3_CVCDnew.txt", header = T,
                       row.names = 1, sep="\t")
metadata$celltype <- factor(metadata$celltype, c('CV','CD'))
condition <- data.frame(class= metadata$celltype)
data.log10 <- log10(data + 1)

pdf("~/Desktop/oocytes_analyses/190307/CVCDtop20_ave.pdf",height=5, width = 16)

rownames(condition) <- colnames(data)
heatplot(data, dend = "row",
         cols.default = T, lowcol = "blue", highcol = "yellow", scale = "row",
         method="ave", dualScale=T, zlim=c(-2,2),scaleKey=TRUE, classvec2 = condition$class,
         classvecCol = c("green4","orange3"), main = "TOP20 DEGs CD+Veh vs CD+DHT", margins = c(9,12))
legend("topright",legend = c("CV","CD"), col=c("green4","orange3"),lty=1.5,cex = 0.8, lwd = 5)
dev.off()



y <- read.table("~/Desktop/oocytes_analyses/counts_CVHVnew.txt",header = T, row.names =NULL,
                sep="\t")
y <- y[!duplicated(data.frame(y[,1])),]
rownames(y) <- y[,1]
y <- y[,2:ncol(y)]
library(dplyr)
commongene1 <- read.table("~/Desktop/oocytes_analyses/top20CVHV.txt",header = F, row.names = NULL,
                         sep="\t")
data1 <- y[which(rownames(y) %in% commongene1), ]
metadata1 <- read.table(file = "~/Desktop/oocytes_analyses/sampleinfoF1F2F3_CVHVnew.txt", header = T,
                       row.names = 1, sep="\t")
metadata1$celltype <- factor(metadata1$celltype, c('CV','HV'))
condition1 <- data.frame(class= metadata1$celltype)
data1.log10 <- log10(data1 + 1)

pdf("~/Desktop/oocytes_analyses/190307/CVHVtop20_s.pdf",height=5, width = 16)

rownames(condition1) <- colnames(data1)
heatplot(data1, dend = "row",
         cols.default = T, lowcol = "blue", highcol = "yellow", scale = "row",
         method="ave", dualScale=T, zlim=c(-2,2),scaleKey=TRUE, classvec2 = condition1$class,
         classvecCol = c("green4","purple4"), main = "TOP20 DEGs HFHS+Veh vs HFHS+DHT", margins = c(9,12))
legend("topright",legend = c("CV","HV"), col=c("green4","purple4"),lty=1.5,cex = 0.8, lwd = 5)
dev.off()


