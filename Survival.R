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


gene_set_xls <- function(DEG_file='', gene_file='', merge_method = 'intersect'){
  # 差异基因 Gene Symbol
  degs <- data.table::fread(DEG_file, data.table = F)
  degs <- degs[[1]]
  degs <- gsub(degs, pattern = '-', replacement = '_')
  
  # 特定基因集 Gene Symbol
  if (gene_file != ''){
    gene_set <- readxl::read_excel(gene_file)
    gene_set <- gene_set[[1]]
    gene_set <- gsub(gene_set, pattern = '-', replacement = '_')
  } else{
    gene_set <- NULL
  }
  
  # gene_list
  if (merge_method == 'intersect') {
    gene_list <- intersect(degs, gene_set)
  } else if (merge_method == 'union') {
    gene_list <- union(degs, gene_set)
  } else if (merge_method == 'setdiff') {
    gene_list <- setdiff(degs, gene_set)
  } else {
    print('Please select valid merge_method: intersect, union or setdiff')
  }
  
  return(gene_list)
}


gene_set_txt <- function(DEG_file='', gene_file='', merge_method = 'intersect'){
  # 差异基因 Gene Symbol
  degs <- data.table::fread(DEG_file, data.table = F)
  degs <- degs[[1]]
  degs <- gsub(degs, pattern = '-', replacement = '_')
  
  # 特定基因集 Gene Symbol
  if (gene_file != ''){
    gene_set <- data.table::fread(gene_file, data.table = F)
    gene_set <- gene_set[[1]]
    gene_set <- gsub(gene_set, pattern = '-', replacement = '_')
  } else{
    gene_set <- NULL
  }
  
  # gene_list
  if (merge_method == 'intersect') {
    gene_list <- intersect(degs, gene_set)
  } else if (merge_method == 'union') {
    gene_list <- union(degs, gene_set)
  } else if (merge_method == 'setdiff') {
    gene_list <- setdiff(degs, gene_set)
  } else {
    print('Please select valid merge_method: intersect, union or setdiff')
  }
  
  return(gene_list)
}


survival_tpm <- function(gene_list=gene_list, survival_file='', tpm_file=''){
  library(data.table)
  library(stringr)
  
  # 读取生存信息,差异基因,TPM矩阵
  # 构建生存信息和TPM的矩阵
  survival <- fread(survival_file, data.table = F)
  tpm <- fread(tpm_file, data.table = F)
  row.names(tpm) <- tpm[, 1]
  tpm <- tpm[, -1]
  t_tpm <- as.data.frame(t(tpm))
  survival <- survival[match(row.names(t_tpm), survival$sample),]
  colnames(survival) <- c('sample', 'status', 'patient', 'time')
  survival_matrix <- cbind(survival, t_tpm)
  row.names(survival_matrix) <- survival_matrix[, 1]
  survival_matrix <- survival_matrix[, -c(1,3)]
  colnames(survival_matrix) <- gsub(colnames(survival_matrix), pattern = '-', replacement = '_')
  
  # 删除正常样本
  group_list <- ifelse(as.numeric(str_sub(row.names(survival_matrix), 14, 15))<10, TRUE, FALSE)
  survival_matrix <- survival_matrix[group_list, ]
  
  return(survival_matrix)
}
