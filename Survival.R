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


survival_uni_cox <- function(gene_list, survival_matrix){
  library(survival)
  
  # 对应gene list的矩阵
  # 构建生存信息和gene list的矩阵
  gene_list <- intersect(colnames(survival_matrix), gene_list)
  
  uni_cox <- function(single_gene) {
    formula <- as.formula(paste('Surv(time, status)~', single_gene))
    surv_uni_cox <- summary(coxph(formula, survival_matrix))
    ph_hypothesis_p <- try({cox.zph(coxph(formula, survival_matrix))$table[1, 3]}, silent = TRUE)
    if('try-error' %in% class(ph_hypothesis_p)){
      ph_hypothesis_p <- 0
    }
    if (surv_uni_cox$coefficients[, 5] < 0.05 & ph_hypothesis_p > 0.05){
      single_gene_report <- data.frame('Gene_symbol'=single_gene,
                                       'beta'=surv_uni_cox$coefficients[, 1],
                                       'Hazard_Ratio'=exp(surv_uni_cox$coefficients[, 1]),
                                       'z_pvalue'=surv_uni_cox$coefficients[, 5],
                                       'Wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
                                       'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
      
      single_gene_report
    }
  }
  
  gene_report <- lapply(gene_list, uni_cox)
  uni_cox_report <- do.call(rbind, gene_report)
  
  return(uni_cox_report)
}


# Old Version
# survival_uni_cox <- function(gene_list, survival_matrix){
#   library(survival)
# 
#   # 对应gene list的矩阵
#   # 构建生存信息和gene list的矩阵
#   gene_list <- intersect(colnames(survival_matrix), gene_list)
#   formula_list <- sapply(gene_list, function(gene_symbol) as.formula(paste('Surv(time, status)~', gene_symbol)))
#   uni_cox_list <- lapply(formula_list, function(formula){coxph(formula, survival_matrix)})
#   
#   ph_hypothesis_p <- c()
#   for (i in 1:length(uni_cox_list)){
#     ph_hypothesis <- try({cox.zph(uni_cox_list[[i]])$table[1, 3]}, silent = TRUE)
#     if('try-error' %in% class(ph_hypothesis)){
#       ph_hypothesis <- 0
#     }
#     ph_hypothesis_p <- append(ph_hypothesis_p, ph_hypothesis)
#   }
# 
#   uni_cox <- function(x) {
#     surv_uni_cox <- summary(x)
#     gene_report <- data.frame('beta'=surv_uni_cox$coefficients[, 1],
#                               'Hazard_Ratio'=exp(surv_uni_cox$coefficients[, 1]),
#                               'z_pvalue'=surv_uni_cox$coefficients[, 5],
#                               'Wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
#                               'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]),
#                               'ph_hypothesis_p'=ph_hypothesis_p)
# 
#     return(gene_report)
#   }
#   
#   gene_report <- lapply(uni_cox_list, uni_cox)
#   uni_cox_report <- do.call(rbind, gene_report)
# 
#   return(uni_cox_report)
# }
