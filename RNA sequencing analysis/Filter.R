sample_filter <- function(sample_file=''){
  # 删除缺失信息的样本
  sample_id <- readxl::read_excel(sample_file)
  sample_list <- sample_id[[1]]
  for (i in 1:length(sample_id[1, ])){
    if (i != 1){
      sample_list <- intersect(sample_list, sample_id[[i]])
    }
  }
  
  return(sample_list)
}


count_sample_filter <- function(sample_list, fpkm_file=''){
  library(data.table)
  
  # 删除缺失信息的样本
  fpkm <- fread(fpkm_file, data.table = F)
  fpkm <- fpkm[,c(TRUE, colnames(fpkm)[-1] %in% sample_list)]
  
  return(fpkm)
}


survival_sample_filter <- function(sample_list, survival_file=''){
  library(data.table)
  
  # 删除缺失信息的样本
  survival <- fread(survival_file, data.table = F)
  survival <- survival[survival$sample %in% sample_list,]
  
  return(survival)
}


phenotype_sample_filter <- function(sample_list, phenotype_file=''){
  library(data.table)
  
  # 删除缺失信息的样本
  phenotype <- fread(phenotype_file, data.table = F)
  phenotype <- phenotype[phenotype$submitter_id.samples %in% sample_list,]
  
  return(phenotype)
}


count_gene_filter <- function(fpkm, filter_method = 'half_zero'){
  if (filter_method == 'half_zero'){
    # 保留在一半样本以上表达的基因
    fpkm <- fpkm[apply(fpkm[-1], 1, function(x) sum(x > 0) > 0.5*ncol(fpkm)), ]
  } else if (filter_method == 'zero'){
    # 去除所有样本表达全为0的基因
    fpkm <- fpkm[rowSums(fpkm[-1])>0, ]
  }
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
