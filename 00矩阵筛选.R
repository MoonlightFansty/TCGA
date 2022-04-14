count_survival_filter <- function(fpkm_file='', survival_file=''){
  # 过滤缺失生存信息FPKM矩阵中的样本
  library(stringr)
  fpkm <- read.table(fpkm_file, header = T, row.names = 1)
  survival <- read.table(survival_file, header = T)
  survival$count_sample <- apply(matrix(survival[,1]), 1, function(str){str_replace_all(str, c('-' = '.'))})
  survival <- survival[survival$count_sample %in% colnames(fpkm),]
  fpkm <- fpkm[,colnames(fpkm) %in% survival$count_sample]
  Ensembl_ID <- row.names(fpkm)
  fpkm <- cbind(Ensembl_ID, fpkm)
  survival_filter <- list(fpkm = fpkm, survival = survival)
  
  return(survival_filter)
}

survival_filter <- count_survival_filter('TCGA-BRCA.htseq_fpkm.tsv', 'TCGA-BRCA.survival.tsv')
fpkm <- survival_filter$fpkm
survival <- survival_filter$survival

write.table(survival, file='TCGA-BRCA.survival.tsv',sep='\t', quote = F, row.names = F, col.names = T)
write.table(fpkm, file='TCGA-BRCA.htseq_fpkm.tsv',sep='\t', quote = F, row.names = F, col.names = T)











