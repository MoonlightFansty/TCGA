# TCGA-RNAseq
## Filter
***filtering TCGA or GDC data from UCSC xena \
including FPKM and survival filtering***

代码示例: \
sample_list <- sample_filter(sample_file = 'TCGA-sample.xlsx') \
fpkm <- count_filter(sample_list, fpkm_file='TCGA-BRCA.htseq_fpkm.tsv') \
survival <- survival_filter(sample_list, survival_file='TCGA-BRCA.survival.tsv')

fwrite(fpkm, file='TCGA-BRCA.htseq_fpkm.tsv', sep='\t') \
fwrite(survival, file='TCGA-BRCA.survival.tsv', sep='\t')

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
***Univariate Cox regression \
including Gene set and survival matrix***

代码示例: \
gene_list <- gene_set_xls('DEG.txt', 'Gene_set.xls', merge_method = 'union') \
survival_matrix <- survival_tpm(survival_file = 'TCGA-BRCA.survival.tsv', tpm_file = 'TPM.txt') \
uni_cox_report <- survival_uni_cox(gene_list, survival_matrix)

fwrite(uni_cox_report, file='Uni_Cox.txt', sep='\t')
