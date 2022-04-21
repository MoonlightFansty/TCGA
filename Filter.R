sample_filter <- function(sample_file=''){
  # 删除缺失信息的样本
  sample_id <- readxl::read_excel(sample_file)
  sample_list <- sample_id[[1]]
  for (i in length(sample_id[1, ])){
    if (i != 1){
      sample_list <- intersect(sample_list, sample_id[[i]])
    }
  }
  
  return(sample_list)
}


survival_filter <- function(sample_list, survival_file=''){
  library(data.table)
  
  # 删除缺失信息的样本
  survival <- fread(survival_file, data.table = F)
  survival <- survival[survival$sample %in% sample_list,]
  
  return(survival)
}


count_filter <- function(sample_list, fpkm_file=''){
  library(data.table)
  
  # 删除缺失信息的样本
  fpkm <- fread(fpkm_file, data.table = F)
  fpkm <- fpkm[,c(TRUE, colnames(fpkm)[-1] %in% sample_list)]
  
  return(fpkm)
}


# Old Version
# count_survival_filter <- function(fpkm_file='', survival_file=''){
#   library(data.table)
#   
#   # 过滤缺失生存信息FPKM矩阵中的样本
#   fpkm <- fread(fpkm_file, data.table = F)
#   survival <- fread(survival_file, data.table = F)
#   survival <- survival[survival$sample %in% colnames(fpkm),]
#   fpkm <- fpkm[,c(TRUE, colnames(fpkm)[-1] %in% survival$sample)]
#   survival_filter <- list(fpkm = fpkm, survival = survival)
#   
#   return(survival_filter)
# }
