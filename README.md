# pyclone_clonevol_analysis
## Intra-tumor heterogeneity with Pyclone and ClonEvol

## Description

R scripts to perform analysis of intra-tumor heterogeneity using Pyclone to estimate the clusters (from both CNVs estimated with facets and SNVs estimates with MuTect2 previously) and ClonEvol to reconstruct the clonal evolution. The outputs are both clusters of mutations and tree of evolution.  

## Parameters

To obtain details on input parameters, user should ask the help, _e.g._:

```
Rscript reformat4pyclone.r --help
```

## Dependencies

1. [R](https://www.r-project.org/)
2. [PyClone](https://bitbucket.org/aroth85/pyclone/wiki/Installation)
3. [ClonEvol](https://github.com/hdng/clonevol)

## Script reformat4pyclone.r

### Usage
  ```
  Rscript reformat4pyclone.r --VCF=SAMPLE1_PASS.vcf.hg38_multianno.vcf.gz --CNV_tumor1=TUMOR1.csv.gz_CNV.txt --CNV_tumor2=TUMOR2.csv.gz_CNV.txt --normal_id=XB00JALY --tumor1_id=XB00JALZ --tumor2_id=XB00JAM0 --kept_genes=LATS2,BAP1,NF2,TP53,SETD2,SETDB1,DDX3X,ULK2,RYR2,DDX51 --min_cnv_length=100000 --min_dp=70
  ```

### Output
  | Type      | Description     |
  |-----------|---------------|
  | file    | tab file (e.g. reformated4pyclone.csv) with columns: mutation_id	ref_counts	var_counts	normal_cn	minor_cn	major_cn	gene |

### After running the script reformat4pyclone.r, user should run Pyclone to obtain a folder input for clonevol (e.g. a folder pyclone_output_SAMPLE1)

## Script clonevol.r

### Usage
  ```
  Rscript clonevol.r --pyclone_res_folder=pyclone_output_SAMPLE1 --tumor1_reformated=XB00JALZ_reformated4pyclone.csv --tumor1_id=XB00JALZ tumor2_id=XB00JAM0 --kept_genes=LATS2,BAP1,NF2,TP53,SETD2,SETDB1,DDX3X,ULK2,RYR2,DDX51
  ```

### Output
  | Type      | Description     |
  |-----------|---------------|
  | folder    | folder (e.g. output_pyclone_XB00JAJB) containing the cluster plot and the tree model plots |




## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Tiffany Delhomme*    |            delhommet@students.iarc.fr | Developer to contact for support |
