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

