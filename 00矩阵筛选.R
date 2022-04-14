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
