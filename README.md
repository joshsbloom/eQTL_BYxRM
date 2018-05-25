# Genetics of trans-regulatory variation in gene expression
### code and data for 1000 BYxRM haploid yeast segregant eQTL mapping by Bloom and Albert

see [eQTL_BYxRM1000_stranded.R](code/eQTL_BYxRM1000_stranded.R) for main analysis script
___

genotype and raw transcript count data can be found [here](RData/counts.RData) as an RData object
```r
#After cloning git repo, cd to local directory, load R, and then
load('RData/counts.RData')
# Structure of phenotype data (counts$pheno)
# 6,713 transcripts (rows) X 1012 haploid segregants (columns)
#
# Structure of genotype data (counts$gdata)
# 1012 haploid segregants (rows) X 42,052 markers (columns)
# -1 indicates BY allele, +1 indicates RM allele
# column names indicate chromosome:position_BYvariant/RMvariant
# postitions are based on the S.Cerevisiae SacCer3 genome build
R> str(counts)
List of 2
 $ pheno: num [1:6713, 1:1012] 550 209 141 442 35 640 89 110 307 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:6713] "YHR055C" "YPR161C" "YOL138C" "YDR395W" ...
  .. ..$ : chr [1:1012] "A01_01-A01-A1-BYxRM_eQTL_10-H6" "A01_02-A01-A2-BYxRM_eQTL_11-F3" "A01_03-A01-A3-BYxRM_eQTL_11-C6" "A01_04-A01-A4-BYxRM_eQTL_03-G2" ...
 $ gdata: num [1:1012, 1:42052] 1 -1 1 1 -1 -1 -1 -1 -1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:1012] "A01_01" "A01_02" "A01_03" "A01_04" ...
  .. ..$ : chr [1:42052] "chrI:33040_A/G" "chrI:33048_A/C" "chrI:33070_A/T" "chrI:33077_G/A" ...
```
processed and filtered tpm values per transcript can be found [here] (RData/log2_t.tpm.matrix.RData) as an RData object

```r
load('RData/log2_t.tpm.matrix.RData')
R> str(t.tpm.matrix)
 num [1:1012, 1:5720] 4.87 6.09 5.64 1.88 6.27 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1012] "A01_01-A01-A1-BYxRM_eQTL_10-H6" "A01_02-A01-A2-BYxRM_eQTL_11-F3" "A01_03-A01-A3-BYxRM_eQTL_11-C6" "A01_04-A01-A4-BYxRM_eQTL_03-G2" ...
  ..$ : chr [1:5720] "YAL062W" "YAL061W" "YAL060W" "YAL059W" ...
```
