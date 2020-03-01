install.packages("Seurat")
library(Seurat)
df <- read.table("~/Desktop/oocytes_analyses/190225/counts_CVCDnew.txt",header = T,
                      row.names = NULL,sep="\t")
df <- df[!duplicated(data.frame(df[,1])),]
rownames(df) <- df[,1]
df <- df[,2:ncol(df)]
qc_star_log <- read.table(file = "~/Desktop/oocytes_analyses/190225/sampleinfoF1F2F3_CVCDnew.txt", header = T,
                          row.names = 1, sep="\t")
seurat.obj <- CreateSeuratObject(counts = df, min.cells = 3, min.features = 2, project = "anna")
seurat.obj <- AddMetaData(object = seurat.obj, metadata = qc_star_log)
seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- ScaleData(object = seurat.obj, vars.to.regress = c("lab",'nCount_RNA','tech'))
seurat.obj <- FindVariableFeatures(object = seurat.obj, mean.function = ExpMean, 
                                        dispersion.function = LogVMR, x.low.cutoff = 0.0125, 
                                        x.high.cutoff = 5, y.cutoff = 0.5)
seurat.obj <- RunPCA(object =seurat.obj, pc.genes = seurat.obj@var.genes,npcs = 5,
                          do.print = TRUE, pcs.print = 1:5,genes.print = 3, pcs.compute=5)
x <- factor(seurat.obj@meta.data$celltype,c('CV','CD'))

seurat.obj@meta.data$celltype <- x
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:5)
DimPlot(object = seurat.obj, reduction = "pca", group.by = 'celltype', pt.size = 1.5)
seurat.obj<- SetIdent(seurat.obj,value = 'celltype')
diff_table <- FindAllMarkers(seurat.obj,logfc.threshold = 0.2,only.pos = T,test.use = "DESeq2")
VlnPlot(test.seurat.obj,diff_table$gene[diff_table$cluster=='CV'][1:12])
write.table(diff_table,"~/Desktop/oocytes_analyses/DEGsCVCD.txt", row.names = T, sep = "\t")

df <- read.table("~/Desktop/oocytes_analyses/190225/counts_CVHVnew.txt",header = T,
                 row.names = NULL,sep="\t")
df <- df[!duplicated(data.frame(df[,1])),]
rownames(df) <- df[,1]
df <- df[,2:ncol(df)]
qc_star_log <- read.table(file = "~/Desktop/oocytes_analyses/190225/sampleinfoF1F2F3_CVHVnew.txt", header = T,
                          row.names = 1, sep="\t")
seurat.obj <- CreateSeuratObject(counts = df, min.cells = 3, min.features = 2, project = "anna")
seurat.obj <- AddMetaData(object = seurat.obj, metadata = qc_star_log)
seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- ScaleData(object = seurat.obj, vars.to.regress = c("lab",'nCount_RNA','tech'))
seurat.obj <- FindVariableFeatures(object = seurat.obj, mean.function = ExpMean, 
                                   dispersion.function = LogVMR, x.low.cutoff = 0.0125, 
                                   x.high.cutoff = 5, y.cutoff = 0.5)
seurat.obj <- RunPCA(object =seurat.obj, pc.genes = seurat.obj@var.genes,npcs = 5,
                     do.print = TRUE, pcs.print = 1:5,genes.print = 3, pcs.compute=5)
x <- factor(seurat.obj@meta.data$celltype,c('CV','HV'))

seurat.obj@meta.data$celltype <- x
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:5)
DimPlot(object = seurat.obj, reduction = "pca", group.by = 'celltype', pt.size = 1.5)
seurat.obj<- SetIdent(seurat.obj,value = 'celltype')
diff_table <- FindAllMarkers(seurat.obj,logfc.threshold = 0.2,only.pos = T,test.use = "DESeq2")
VlnPlot(test.seurat.obj,diff_table$gene[diff_table$cluster=='CV'][1:12])
write.table(diff_table,"~/Desktop/oocytes_analyses/DEGsCVHV.txt", row.names = T, sep = "\t")






















