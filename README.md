# TCGA-RNAseq
## FPKM
***pre-processing TCGA or GDC RNAseq data(FPKM) from UCSC xena \
including FPKM normalization and Wilcoxon rank sum test***

代码示例: \
*matrix get max gene expression which have multi-rows* \
count_norm <- count_max_tpm('TCGA-BRCA.htseq_fpkm.tsv', 'gencode.v22.annotation.gene.probeMap') \
out <- count_Wilcoxon_test(count_norm) 

*matrix get mean gene expression which have multi-rows* \
count_norm <- count_mean_tpm('TCGA-BRCA.htseq_fpkm.tsv', 'gencode.v22.annotation.gene.probeMap') \
out <- count_Wilcoxon_test(count_norm) 

count_norm <- as.data.frame(count_norm) \
fwrite(count_norm, file = 'TPM.txt', sep = '\t', row.names = T) \
fwrite(out, file = 'DEG.txt', sep = '\t', row.names = T)

## Survival
***pre-processing TCGA or GDC RNAseq data(FPKM) from UCSC xena \
including survival sample filtering***

代码示例: \
*filter samples which have both survival and expression data* \
survival_filter <- count_survival_filter('TCGA-BRCA.htseq_fpkm.tsv', 'TCGA-BRCA.survival.tsv') \
fpkm <- survival_filter$fpkm \
survival <- survival_filter$survival

fwrite(survival, file='TCGA-BRCA.survival1.tsv', sep='\t') \
fwrite(fpkm, file='TCGA-BRCA.htseq_fpkm1.tsv', sep='\t')

***Gene set and survival matrix***

代码示例: \
*select genes from Gene set and DEG* \
gene_list <- gene_set_xls('DEG.txt', 'Gene_set.xls', merge_method = 'union')

*matrix with survival and TPM* \
survival_matrix <- survival_tpm(survival_file = 'TCGA-BRCA.survival.tsv', tpm_file = 'TPM.txt')

***Univariate Cox regression***

代码示例: \
uni_cox_report <- survival_uni_cox(gene_list, survival_matrix)

fwrite(uni_cox_report, file='uni_cox.txt', sep='\t')
