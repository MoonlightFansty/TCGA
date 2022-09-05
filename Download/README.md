# GDC Download
## 一、数据下载
***download TCGA or GDC data from UCSC xena*** 
TCGA的33种癌症数据列表
```
mkdir TCGA
cd TCGA

cat >cancer_list.txt
GDC TCGA Acute Myeloid Leukemia (LAML) (15 datasets)
GDC TCGA Adrenocortical Cancer (ACC) (14 datasets)
GDC TCGA Bile Duct Cancer (CHOL) (14 datasets)
GDC TCGA Bladder Cancer (BLCA) (14 datasets)
GDC TCGA Breast Cancer (BRCA) (20 datasets)
GDC TCGA Cervical Cancer (CESC) (14 datasets)
GDC TCGA Colon Cancer (COAD) (15 datasets)
GDC TCGA Endometrioid Cancer (UCEC) (15 datasets)
GDC TCGA Esophageal Cancer (ESCA) (14 datasets)
GDC TCGA Glioblastoma (GBM) (15 datasets)
GDC TCGA Head and Neck Cancer (HNSC) (14 datasets)
GDC TCGA Kidney Chromophobe (KICH) (14 datasets)
GDC TCGA Kidney Clear Cell Carcinoma (KIRC) (15 datasets)
GDC TCGA Kidney Papillary Cell Carcinoma (KIRP) (15 datasets)
GDC TCGA Large B-cell Lymphoma (DLBC) (14 datasets)
GDC TCGA Liver Cancer (LIHC) (14 datasets)
GDC TCGA Lower Grade Glioma (LGG) (14 datasets)
GDC TCGA Lung Adenocarcinoma (LUAD) (15 datasets)
GDC TCGA Lung Squamous Cell Carcinoma (LUSC) (15 datasets)
GDC TCGA Melanoma (SKCM) (14 datasets)
GDC TCGA Mesothelioma (MESO) (14 datasets)
GDC TCGA Ocular melanomas (UVM) (14 datasets)
GDC TCGA Ovarian Cancer (OV) (15 datasets)
GDC TCGA Pancreatic Cancer (PAAD) (14 datasets)
GDC TCGA Pheochromocytoma & Paraganglioma (PCPG) (14 datasets)
GDC TCGA Prostate Cancer (PRAD) (14 datasets)
GDC TCGA Rectal Cancer (READ) (15 datasets)
GDC TCGA Sarcoma (SARC) (14 datasets)
GDC TCGA Stomach Cancer (STAD) (15 datasets)
GDC TCGA Testicular Cancer (TGCT) (14 datasets)
GDC TCGA Thymoma (THYM) (14 datasets)
GDC TCGA Thyroid Cancer (THCA) (14 datasets)
GDC TCGA Uterine Carcinosarcoma (UCS) (14 datasets)
#ctrl+D 结束文件

perl -alne '{/\((.*?)\)/;print $1}' cancer_list.txt |while read id;do 
nohup wget https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-$id.{mutect2_snv,GDC_phenotype,survival,htseq_counts}.tsv.gz &
done
```
下载编码和非编码信息 
```
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz
zcat gencode.v22.annotation.gtf.gz|perl -alne '{print if $F[2] eq "gene"}'|cut -d"\"" -f 2,4|sed 's/"/\t/g' > ensembID2type.txt
```
批量读取下载好的表达量矩阵，并且根据基因注释信息，拆分成为蛋白编码基因和非编码基因这两个不同的表达量矩阵
```
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table)  
dir.create('Rdata')

# dataset: gene expression RNAseq - HTSeq - Counts 
# unitlog2(count+1)
fs=list.files('ucsc_xena',pattern = 'htseq_counts.tsv.gz')
fs
lapply(fs, function(x){
  # x=fs[1]
  pro=gsub('tsv.gz','',x)
  print(pro)
  a=fread(file.path('ucsc_xena',x),data.table = F)
  head(a[ ,1:4])
  tail(a[ ,1:4])
  mat=a[1:60483,] 
  rownames(mat)=mat$Ensembl_ID
  mat[1:4,1:4]
  mat=mat[,-1]

  exprSet=floor(2^mat - 1 )
  exprSet[1:4,1:4] 

  n=floor(ncol(exprSet)/20)
  n
  keep_feature <- rowSums (exprSet > 2) > n
  table(keep_feature)
  mat <- exprSet[keep_feature, ]
  mat=mat[, colSums(mat) > 1000000]
  mat[1:4,1:4]
  dim(mat) 
  colnames(mat)

  b=read.table('ucsc_xena/gencode.v22.annotation.gene.probeMap',header = T)
  head(b)
  d=read.table('ucsc_xena/ensembID2type.txt')
  head(d)
  b=merge(b,d,by.x='id',by.y='V1')
  head(b)

  length(unique(b[match(rownames(mat),b$id),2]))
  # 可以看到，多个ensembl的ID对应同一个基因的symbol，这个是需要处理掉的。
  # 下面的代码略微有一点点复杂
  dat=mat
  ids <-  b[match(rownames(dat),
                  b$id),1:2] #取需要的列 
  head(ids)
  colnames(ids)=c('probe_id','symbol')  
  ids=ids[ids$symbol != '',]
  ids=ids[ids$probe_id %in%  rownames(dat),]
  dat[1:4,1:4]   
  dat=dat[ids$probe_id,] 

  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  dat[1:4,1:4]  #保留每个基因ID第一次出现的信息

  mat=dat
  mat[1:4,1:4]
  print(fivenum(mat['GAPDH',]))
  print(fivenum(mat['ACTB',]))

  tp=b[match(rownames(mat),b$gene),7]
  print(  tail(sort(table(tp))) )

  pd_mat=mat[tp=='protein_coding',]
  non_mat=mat[tp !='protein_coding',]

  save(pd_mat,non_mat,
       file = file.path('Rdata', paste0(pro,'Rdata')))

})
```
查看表达矩阵里面的存储的蛋白编码基因和非编码基因这两个不同的表达量矩阵的基因数量和病人数量
```
fs=list.files('Rdata/',pattern = 'htseq_counts')
fs
do.call(rbind,lapply(fs, function(x){
  load(file =  file.path('Rdata/',x))  
  return(c(x,dim(pd_mat),dim(non_mat)))
}))
```
