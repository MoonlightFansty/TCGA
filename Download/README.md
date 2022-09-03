# GDC Download
## 一、数据下载
***download TCGA or GDC data from UCSC xena***

代码示例: \
sample_list <- sample_filter(sample_file = 'TCGA-LUAD-sample.xlsx') \
fpkm <- count_sample_filter(sample_list, fpkm_file='TCGA-LUAD.htseq_fpkm.tsv') \
fpkm <- count_gene_filter(fpkm, filter_method = 'half_zero') \
survival <- survival_sample_filter(sample_list, survival_file='TCGA-LUAD.survival.tsv') \
phenotype <- phenotype_sample_filter(sample_list, phenotype_file='TCGA-LUAD.GDC_phenotype.tsv' )

fwrite(fpkm, file='TCGA-LUAD.htseq_fpkm.tsv', sep='\t') \
fwrite(survival, file='TCGA-LUAD.survival.tsv', sep='\t') \
fwrite(phenotype, file='TCGA-LUAD.GDC_phenotype.tsv', sep='\t')




