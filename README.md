# TCGA-RNAseq
## FPKM
pre-processing TCGA or GDC RNAseq data(FPKM) from UCSC xena \
including FPKM normalization and Wilcoxon rank sum test
-------------
代码示例: \
count_norm <- count_max_tpm('TCGA-BRCA.htseq_fpkm.tsv', 'gencode.v22.annotation.gene.probeMap') \
out <- count_Wilcoxon_test(count_norm1) \

