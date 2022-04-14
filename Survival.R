count_survival_filter <- function(fpkm_file='', survival_file=''){
  library(data.table)
  
  # 过滤缺失生存信息FPKM矩阵中的样本
  fpkm <- fread(fpkm_file, data.table = F)
  survival <- fread(survival_file, data.table = F)
  survival <- survival[survival$sample %in% colnames(fpkm),]
  fpkm <- fpkm[,c(TRUE, colnames(fpkm)[-1] %in% survival$sample)]
  survival_filter <- list(fpkm = fpkm, survival = survival)
  
  return(survival_filter)
}
