# TCGA-RNAseq
## FPKM
***pre-processing TCGA or GDC RNAseq data(FPKM) from UCSC xena \
including FPKM normalization and Wilcoxon rank sum test***

代码示例: \
*matrix get max genes expression which have multi-rows* \
count_norm <- count_max_tpm('TCGA-BRCA.htseq_fpkm.tsv', 'gencode.v22.annotation.gene.probeMap') \
out <- count_Wilcoxon_test(count_norm) 

*matrix get mean genes expression which have multi-rows* \
count_norm <- count_mean_tpm('TCGA-BRCA.htseq_fpkm.tsv', 'gencode.v22.annotation.gene.probeMap') \
out <- count_Wilcoxon_test(count_norm) 

# Survival
***pre-processing TCGA or GDC RNAseq data(FPKM) from UCSC xena \
including survival sample filtering***

代码示例: \
*filter samples which have both survival and expression data* \
survival_filter <- count_survival_filter('TCGA-BRCA.htseq_fpkm.tsv', 'TCGA-BRCA.survival.tsv') \
fpkm <- survival_filter$fpkm \
survival <- survival_filter$survival

write.table(survival, file='TCGA-BRCA.survival.tsv',sep='\t', quote = F, row.names = F, col.names = T) \
write.table(fpkm, file='TCGA-BRCA.htseq_fpkm.tsv',sep='\t', quote = F, row.names = F, col.names = T)
