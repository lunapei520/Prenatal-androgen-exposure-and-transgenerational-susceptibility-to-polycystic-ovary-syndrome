qc_star_log <- read.table(file = "~/Desktop/oocytes_analyses/190225/sampleinfoF1F2F3_CVHVCDnew.txt", header = T,
                          row.names = 1, sep="\t")
rpkm_matrix <- read.table(file = "~/Desktop/oocytes_analyses/190225/counts_F1F2F3CVCDHVnew.txt",header = T,
                          row.names = NULL,sep="\t")
map_ratio_reads <- (qc_star_log[,3]>60)&(qc_star_log[,1]>1000000)
rpkm_matrix_f1 <- t(t(rpkm_matrix)[map_ratio_reads,])
qc_star_log_f1 <- qc_star_log[map_ratio_reads,]
gene_expression_cell_num_1 <- NULL
for (i in 1:nrow(rpkm_matrix_f1)) {
  gene_expression_cell_num_1[i] <- sum(rpkm_matrix_f1[i,]>0)
}
rpkm_matrix_f2 <- rpkm_matrix_f1[(gene_expression_cell_num_1 > 1),]
spear_cor <- 0.4
cell_cor_spearman <- cor(rpkm_matrix_f2, method = 'spearman')
cell_cor_spearman_lower <- cell_cor_spearman
for (i in 1:ncol(cell_cor_spearman_lower)) {
  cell_cor_spearman_lower[i,i] = 0
}
qc_star_log_f3 <- qc_star_log_f1[as.matrix(apply(cell_cor_spearman_lower, 1, max))>spear_cor,]
rpkm_matrix_f3 <- rpkm_matrix_f2[,as.matrix(apply(cell_cor_spearman_lower, 1, max))>spear_cor]
gene_per_cell <- 5000
gene_exp_num <- NULL
for (i in 1:ncol(rpkm_matrix_f3)) {
  gene_exp_num[i] <- sum(rpkm_matrix_f3[,i]>0)
} 
qc_star_log_f4 <- qc_star_log_f3[as.matrix(gene_exp_num>gene_per_cell)]
rpkm_matrix_f4 <- rpkm_matrix_f3[,as.matrix(gene_exp_num>gene_per_cell)]
qc_star_log_QC <- qc_star_log_f4
rpkm <- rpkm_matrix_f4

