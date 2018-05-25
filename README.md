# Genetics of trans-regulatory variation in gene expression
## code and data for 1000 BYxRM haploid yeast segregant eQTL mapping by Bloom and Albert

see [eQTL_BYxRM1000_stranded.R](code/eQTL_BYxRM1000_stranded.R) for main analysis script
___

Genotype data can be found here as an RData object
```r
#After cloning git repo, cd to local directory, load R, and then
load('genotypes/gdata_42k.RData')
# Structure of genotype data
# 1012 haploid segregants (rows) X 42,052 markers (columns)
# -1 indicates BY allele, +1 indicates RM allele
# column names indicate chromosome:position_BYvariant/RMvaiant
# postitions are based on the S.Cerevisiae SacCer3 genome build
R> str(gdata)
 num [1:1012, 1:42052] 1 -1 1 1 -1 -1 -1 -1 -1 1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1012] "A01_01" "A01_02" "A01_03" "A01_04" ...
  ..$ : chr [1:42052] "chrI:33040_A/G" "chrI:33048_A/C" "chrI:33070_A/T" "chrI:33077_G/A" ... 
```
Count data can be found here as an RData object


