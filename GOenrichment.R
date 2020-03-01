install.packages('org.Mm.eg.db')
install.packages('clusterProfiler')
library(clusterProfiler)
library(org.Mm.eg.db)
df <- read.table("~/Desktop/oocytes_analyses/DEGsCVCD.txt",sep = "\t", header =F, 
                 row.names = NULL)
geneList <- as.matrix(df)
symbols <- df[,1]
gene.df <- bitr(gene = symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
ego <- enrichGO(gene          = gene.df$SYMBOL,
                
                OrgDb         = org.Mm.eg.db,
                keyType = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "none",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
barplot(ego2,drop=T, showCategory = 30)
dotplot(ego,showCategory=20)
cnetplot(ego, categorySize="pvalue", foldChange=geneList)

































