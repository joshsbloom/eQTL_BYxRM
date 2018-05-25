#library(prada)
#library(gputools)
#library(gmatrix)
library(gdata)
library(plyr)
library(seqinr)
library(grofit)
library(abind)
library(intervals)
library(GenomicRanges)
library(RcppArmadillo)
library(Rcpp)
library(inline)
library(Matrix)
library(rrBLUP)
library(abind)
library(EBglmnet)
library(caret)
library(qvalue)
library(MASS)
library(WriteXLS)
library(leaps)
library(qtlDesign)
library(foreach)
library(doMC)
registerDoMC(cores=70)
X11.options(type='cairo')
# FRANK, modify base.dir as appropriate
base.dir='/data/eQTL/'

#source(paste0(base.dir, 'mmA.R'))
source(paste0(base.dir, 'code/eQTL_BYxRM1000_fx_stranded.R'))

#unique chromosomes in order
unique.chrs=c(paste0('chr', as.roman(1:16)), 'chrMito')

# load genotype data ---------------------------------------------------------------------------------------------------
#genotype data -1 =BY, +1 =RM
load(paste0(base.dir, 'genotypes/gdata_42k.RData')) #1000BYxRM_with_names.RData'))
BYxRM_orig=gdata 
#rownames are missing 4000BYxRM name

#str(BYxRM_orig)
# num [1:1040, 1:28220] 1 -1 1 1 1 -1 1 -1 1 -1 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:1040] "A01_01:01_01" "A01_02:01_09" "A01_03:01_17" "A01_13:01_21" ...
#  ..$ : chr [1:28220] "33070_chrI_33070_A_T" "33147_chrI_33147_G_T" "33152_chrI_33152_T_C" "33200_chrI_33200_C_T" ...
#-----------------------------------------------------------------------------------------------------------------------

# process growth rate data  --------------------------------------------------------------------------------------------
#growth.fit=getGrowthStats(base.dir)
#save(growth.fit, file=paste0(base.dir, 'RData/liquid_growth.RData'))
load(file=paste0(base.dir, 'RData/liquid_growth.RData'))
mOD=lapply(growth.fit, function(x) sapply(x, function(y) y$maxOD-y$minOD))
#mOD=lapply(mOD, function(x) x-mean(x, na.rm=T)        )
OD.cov=unlist(mOD, use.names=F)
names(OD.cov)=unlist(sapply(mOD, names), use.names=F)
names(OD.cov)=paste(gsub('_1$', '', names(OD.cov)), rep(paste0('BYxRM_eQTL_', sprintf("%02d", 1:13)), each=96), sep='-')
#-----------------------------------------------------------------------------------------------------------------------

#ocd=read.fasta('~/orf_coding.fasta')

#load kallisto output -------------------------------------------------------------------------------------------
# load transcript annotations -----------------------------------------------------------------------------------
    transcript.fa.file=paste0(base.dir, 'reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa')
    transcripts.fa=read.fasta(transcript.fa.file)
    transcript.info=sapply(transcripts.fa, function(x) attr(x,'Annot'))

    dir.in=paste0(base.dir, 'kallisto_out_merged_stranded_boot/')
    counts.list=list()
    for(sub.dir in list.files(dir.in)) {
        print(sub.dir)
        counts.list[[sub.dir]]=read.delim(paste0(dir.in, sub.dir, '/abundance.tsv'), header=T, sep='\t') 
    }
#---------------------------------------------------------------------------------------------------------------------

# create table of read counts per gene  ------------------------------------------------------------------------------
raw.count.matrix = do.call('cbind', lapply(counts.list, function(x) (x$est_counts)))
rownames(raw.count.matrix)=as.character(counts.list[[1]]$target_id)
#---------------------------------------------------------------------------------------------------------------------


#Extract information about transcripts--------------------------------------------------------------------------------
gene.annot.df=buildGeneAnnotationDF(transcript.info,raw.count.matrix,counts.list)
#---------------------------------------------------------------------------------------------------------------------

# reorder count matrix by chromosome and position of genes
count.matrix=raw.count.matrix[order(match(rownames(raw.count.matrix), gene.annot.df$name)),]
#--------------------------------------------------------------------------------------------------------------------

# Filter samples on total reads per sample  -------------------------------------------------------------------------
sample.sum=apply(count.matrix, 2, sum)
d.cutoff=1e6  # adjusted from 3e5 to 5e5 to 1e6
#png(file=paste0('/data/eQTL/plots/reads_per_sample.png'), width=512, height=512)
hist(sample.sum, breaks=150, xlab='reads per sample', main='', sub=paste('median =', round(median(sample.sum))))
# depth cutoff
abline(v=d.cutoff, col='red')
#dev.off()

#for now hard cutoff at 1e6 and 1e7
# also includes high count
low.count.samples=which(sample.sum<d.cutoff )#  | sample.sum>1e7)
#png(file=paste0('/data/eQTL/plots/reads_per_sample_per_batch.png'), width=1024, height=768)
boxplot( sample.sum~sapply(strsplit(names(sample.sum), '-'), function(x) x[4]), main='reads per sample per batch')
#dev.off()
print('total read')
print(length(counts.list))

# Remove low count samples here!
counts.list=counts.list[-low.count.samples]
print('total used')
print(length(counts.list))
count.matrix=count.matrix[,-low.count.samples]
sample.sum = sample.sum[-low.count.samples]
#--------------------------------------------------------------------------------------------------------------------------

# Match genotypes and phenotypes as early as possible ---------------------------------------------------------------------
match_pheno_and_geno=function(pheno_input_matrix, geno_matrix) {
    BYxRM_strain_name=rownames(geno_matrix) #do.call('rbind',strsplit(rownames(geno_matrix), ':'))[,1]
    sname=(do.call('rbind', strsplit(colnames(pheno_input_matrix), '-'))[,1])
    pheno_input_matrix=pheno_input_matrix[,-which(is.na(match(sname,BYxRM_strain_name)))]
    sname=(do.call('rbind', strsplit(colnames(pheno_input_matrix), '-'))[,1])
    gdata=geno_matrix[match(sname,BYxRM_strain_name),]
    return(list(pheno=pheno_input_matrix, gdata=gdata))
}

counts=match_pheno_and_geno(count.matrix, BYxRM_orig)
#save(counts, file='/data/eQTL/RData/counts.RData')
#load('/data/eQTL/RData/counts.RData')
gdata=counts$gdata

# Construct LD for full marker set -----------------------------------------------------------------------------------------
#cvec=c(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])
#unique.chrs=paste0('chr', as.roman(1:16))
#gdata.by.chr=list()
#gdata.s.by.chr=list()
#for(cc in unique.chrs) {   
        #gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
        #gdata.s.by.chr[[cc]]=scale(gdata.by.chr[[cc]])
        
#}
#marker.LD=lapply(gdata.by.chr, cor)
#save(marker.LD, file='/data/eQTL/RData/marker.LD.RData')
#---------------------------------------------------------------------------------------------------------------------------

# remove completely correlated markers
gdata=gdata[,!duplicated(gdata, MARGIN=2)]

# Reconstruct chromosome level lists of markers, remove perfectly correlated markers
cvec=c(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])
unique.chrs=paste0('chr', as.roman(1:16))
gdata.by.chr=list()
gdata.s.by.chr=list()
for(cc in unique.chrs) {   
        gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
        gdata.s.by.chr[[cc]]=scale(gdata.by.chr[[cc]])
        
}
#build genetic map -----------------
genetic.map=buildGeneticMap(gdata.by.chr)
genetic.map.c=as.vector(unlist(genetic.map))
names(genetic.map.c)=as.vector(unlist(sapply(genetic.map, names)))

# remove perfectly correlated markers 
gdata.downsampled=downsampleMarkers(gdata.by.chr, gdata.s.by.chr)
#save(gdata.downsampled, file='/data/eQTL/RData/gdata.downsampled.RData')


# useful for accelerated mapping function
gdata.scaled=scale(gdata)
# rebuild count.matrix for sample names properly matched with genotype data
count.matrix=counts$pheno
#--------------------------------------------------------------------------------------------------------------------------


# code for analysis by batch -----------------------------------------------------------------------------------------------
#g=gdata ; #gs=gdata.scaled; #cm=count.matrix;  #gad=gene.annot.df
#gdata=g; #gdata.scaled=gs; #count.matrix=cm #gene.annot.df=gad
#ss='eQTL_12'; #sub.set=grepl(ss, colnames(count.matrix)) ;#count.matrix=count.matrix[,sub.set] ;#gdata=gdata[sub.set,] ;
#gdata.scaled=gdata.scaled[sub.set,] ;#now subset by eQTL3, 12 and 13 
#----------------------------------------------------------------------------------------------------------------------------

# Build annotation table from sample names ----------------------------------------------------------------------------------
sample.annot.df=buildSampleAnnotationDF(count.matrix)
gbatch.fact=as.factor(sample.annot.df$growth.batch)

mOS=match(paste(sample.annot.df$name, sample.annot.df$growth.batch, sep='-'), names(OD.cov))
OD.cov=OD.cov[mOS]

#new
##mOD.slope=lapply(growth.fit, function(x) sapply(x, function(y) y$A))
##OD.cov.slope=unlist(mOD.slope, use.names=F)
##names(OD.cov.slope)=unlist(sapply(mOD.slope, names), use.names=F)
##names(OD.cov.slope)=paste(gsub('_1$', '', names(OD.cov.slope)), rep(paste0('BYxRM_eQTL_', sprintf("%02d", 1:13)), each=96), sep='-')

#OD.cov.slope=OD.cov.slope[mOS]
##OD.cov.slope[is.na(OD.cov.slope)]=mean(OD.cov.slope, na.rm=T)
##boxplot(OD.cov.slope~gbatch.fact)
# -------------------------------------------------------------------------------------------------------------------------------

# Downsample counts -------------------------------------------------------------------------------------------------------------
#downsample.2e6.count.matrix=downsampleCounts(count.matrix, 2e6, with.replacement=TRUE)
downsample.1e6.count.matrix=downsampleCounts(count.matrix, 1e6)
downsample.5e5.count.matrix=downsampleCounts(count.matrix, 5e5)
downsample.1e5.count.matrix=downsampleCounts(count.matrix, 1e5)
downsample.5e4.count.matrix=downsampleCounts(count.matrix, 5e4)

 #d2e6'=downsample.1e6.count.matrix,

downsample.count.matrices=list('d1e6'=downsample.1e6.count.matrix,
                         'd5e5'=downsample.5e5.count.matrix,
                         'd1e5'=downsample.1e5.count.matrix,
                         'd5e4'=downsample.5e4.count.matrix)
#save(downsample.count.matrices, file='/data/eQTL/RData/downsample.count.matrices.RData') 
#load('/data/eQTL/RData/downsample.count.matrices.RData')

downsample.tpm.matrices=lapply(downsample.count.matrices, function(y) {
                         log2(apply(y,2, function(x) countToTpm(x, gene.annot.df$length))+.5)
                         })

# Calculate TPM -------------------------------------------------------------------------------------------------------------------
# filter out misbehaving transcripts 
# recalculate tpm using actual transcript lengths 
# also, go ahead and calculate log2 transform
tpm.matrix=log2(apply(count.matrix,2, function(x) countToTpm(x, gene.annot.df$length))+0.5)

#covariates.OD=model.matrix(t(tpm.matrix)[,1]~gbatch.fact+OD.cov)
#residual.pheno.OD=scale(residuals(lm(t.tpm.matrix~covariates.OD) ))
#scanoneLODS=fasterLOD(nrow(residual.pheno.OD),residual.pheno.OD,gdata.scaled, betas=TRUE)

#remove invariant transcripts-------------------------------------------------------------------------------------------------------
invariant.inds=(which(apply(tpm.matrix,1,quantile, .99)<1  | apply(tpm.matrix,1,median)<1 ) )

gene.annot.df=gene.annot.df[-invariant.inds,]

tpm.matrix=tpm.matrix[-invariant.inds,]
t.tpm.matrix=t(tpm.matrix)
count.matrix=count.matrix[-invariant.inds,]

downsample.count.matrices=lapply(downsample.count.matrices, function(x) x[-invariant.inds,])
downsample.tpm.matrices=lapply(downsample.tpm.matrices, function(x) x[-invariant.inds,])
# ----------------------------------------------------------------------------------------------------------------------------------

# To estimate H^2 and find parental differences ------------------------------------------------------------------------------------
       # parental values
       dir.in.p=paste0(base.dir, 'kallisto_out_merged_parents_stranded_boot/')
       counts.list.p=list()
       for(sub.dir in list.files(dir.in.p)) {
             print(sub.dir)
             counts.list.p[[sub.dir]]=read.delim(paste0(dir.in.p, sub.dir, '/abundance.tsv'), header=T, sep='\t') 
       }

      #create table of read counts per gene for parents---------------------------------------------------------------
      raw.count.matrix.p = do.call('cbind', lapply(counts.list.p, function(x) (x$est_counts)))
      rownames(raw.count.matrix.p)=as.character(counts.list.p[[1]]$target_id)

      raw.count.matrix.p=raw.count.matrix.p[rownames(raw.count.matrix.p) %in% colnames(t.tpm.matrix),]
      count.matrix.p=raw.count.matrix.p[order(match(rownames(raw.count.matrix.p), gene.annot.df$name)),]
      p.read.cnt=apply(count.matrix.p, 2, sum)
      count.matrix.p=count.matrix.p[,p.read.cnt>1e6]
      p.groups=sapply(strsplit(colnames(count.matrix.p), '-'), function(x) x[1])
      p.batch=sapply(strsplit(colnames(count.matrix.p), '-'), function(x) x[4])
      tpm.matrix.p=log2(apply(count.matrix.p,2, function(x) countToTpm(x, gene.annot.df$length))+.5)


      segs.and.par=cbind(tpm.matrix, tpm.matrix.p)
      lm.tpm.gbatch.OD=lm(t(segs.and.par)~as.factor(c(as.character(gbatch.fact), p.batch)))
      residuals.tpm.gbatch.OD=residuals(lm.tpm.gbatch.OD)

      corrected.tpm.matrix.p=t(residuals.tpm.gbatch.OD[grep('1879|1950', rownames(residuals.tpm.gbatch.OD)),])
      corrected.tpm.matrix=t(residuals.tpm.gbatch.OD[!grepl('1879|1950', rownames(residuals.tpm.gbatch.OD)),])

      p.var=t(apply(corrected.tpm.matrix.p, 1, function(x) sapply(split(x,p.groups), var)))
      p.cnts=rle(p.groups)$lengths
      df=p.cnts-1
      pooled.par.var=(p.var[,1]*df[1]+p.var[,2]*df[2])/sum(df)
      pooled.par.var[pooled.par.var==0]=NA
      pooled.par.var[pooled.par.var>1e5]=NA

      p.diff=apply(corrected.tpm.matrix.p, 1, function(x) wilcox.test(x~p.groups)$p.value)
      
      #qs=qvalue(p.diff[!is.na(p.diff)])

      seg.var=apply(corrected.tpm.matrix, 1, var)
      H2=(seg.var-pooled.par.var)/seg.var
      H2[H2<(-35)]=NA
#---------------------------------------------------------------------------------------------------------------


# some more global variables - ---------------------------------------------------
gene.GR=GRanges(seqnames=(gene.annot.df$chr), ranges=IRanges(start=gene.annot.df$start, end=gene.annot.df$end), 
                strand=ifelse(gene.annot.df$strand==1, '+', '-'),
                ORF=(gene.annot.df$name)
                )
ma.chr=sapply(strsplit(colnames(gdata), ':'), function(x)x[1])
ma.pos=as.numeric( sapply( strsplit(sapply(strsplit(colnames(gdata), ':'), function(x)x[2]), '_'), function(x)x[1]))

marker.GR=GRanges(seqnames=ma.chr, ranges=IRanges(start=ma.pos, end=ma.pos) )
#5720
contig.lengths=get.contig.lengths('/data/eQTL/reference/sacCer3.fasta')

gcoord.key= build.gcoord.key('/data/eQTL/reference/sacCer3.fasta')
gene.GR$gcoord=as.vector(gcoord.key[as.character(seqnames(gene.GR))] + start(gene.GR) )
marker.GR$gcoord=as.vector(gcoord.key[as.character(seqnames(marker.GR))] + start(marker.GR) )
marker.GR$mname=colnames(gdata)
gene.GR=gene.GR[order(gene.GR$gcoord),]
#save(marker.GR, file='/data/eQTL/RData/markerAnnotation.RData')
#save(gene.GR, file='/data/eQTL/RData/geneAnnotation.RData')

# Key data objects
#load('/data/eQTL/RData/markerAnnotation.RData')
#load('/data/eQTL/RData/geneAnnotation.RData')

t.tpm.matrix=t.tpm.matrix[,match(gene.GR$ORF, colnames(t.tpm.matrix))]
count.matrix=count.matrix[match(gene.GR$ORF, rownames(count.matrix)),]
downsample.tpm.matrices=lapply(downsample.tpm.matrices, function(x){   x[match(gene.GR$ORF, rownames(x)),]  })

# Key data objects
#save(t.tpm.matrix, file='/data/eQTL/RData/log2_t.tpm.matrix.RData')
#save(count.matrix, file='/data/eQTL/RData/count.matrix.RData')
#load('/data/eQTL/RData/log2_t.tpm.matrix.RData') 
#load('/data/eQTL/RData/count.matrix.RData')
####################str(#load('/data/eQTL/RData/log2_t.tpm.matrix.RData')

A=tcrossprod(gdata.scaled)/ncol(gdata.scaled)

# Additional ideas for covariance matrices
# batch.fact=model.matrix(tpm.matrix[1,]~gbatch.fact-1)
# batch.mat=batch.fact%*%t(batch.fact)
# batch.by.g=batch.mat*A
# batch.by.gg=batch.mat*A*A

# could add in OD.cov or PCs here as well    
#covariates=model.matrix(t.tpm.matrix[,1]~gbatch.fact)
#residual.pheno=residuals(lm(t.tpm.matrix~covariates) )
#pheno.scaled=(scale(residual.pheno))

#mODmod=cbind(1, model.matrix(t.tpm.matrix[,1]~gbatch.fact-1), OD.cov)
#residual.pheno.OD2=residuals(lm(t.tpm.matrix~mODmod-1) )
covariates.OD=model.matrix(t.tpm.matrix[,1]~gbatch.fact+OD.cov)
residual.pheno.OD=residuals(lm(t.tpm.matrix~covariates.OD) )
pheno.scaled.OD=(scale(residual.pheno.OD))
scanoneLODS.OD=fasterLOD(nrow(pheno.scaled.OD),pheno.scaled.OD,gdata.scaled, betas=TRUE)

#save(covariates.OD, file='/data/eQTL/RData/covariates.OD.RData')
#save(scanoneLODS.OD, file = '/data/eQTL/RData/scanoneLODS_OD_stranded.RData')

#R> str(scanoneLODS)
#List of 2
# $ r  : num [1:6290, 1:11530] -0.7073 -0.02137 0.4509 0.00975 -0.01648 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:6290] "YAL062W" "YAL061W" "YAL060W" "YAL059W" ...
#  .. ..$ : chr [1:11530] "chrI:33040_A/G" "chrI:33293_A/T" "chrI:34170_T/A" "chrI:34308_C/T" ...
# $ LOD: num [1:6290, 1:11530] 152.4408 0.1004 49.9469 0.0209 0.0597 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:6290] "YAL062W" "YAL061W" "YAL060W" "YAL059W" ...
#  .. ..$ : chr [1:11530] "chrI:33040_A/G" "chrI:33293_A/T" "chrI:34170_T/A" "chrI:34308_C/T" ...
# note: 'r' here is the beta from a linear model

# Fast code for single component random effect models------------------------------
vcA.OD=calcA(pheno.scaled.OD, A)
h2A.OD=(vcA.OD[,1]/(vcA.OD[,1]+vcA.OD[,2]))
#----------------------------------------------------------------------------------
vcA.OD.unscaled=calcA(residual.pheno.OD, A)
rownames(vcA.OD.unscaled)=colnames(residual.pheno.OD)
colnames(vcA.OD.unscaled)=c('A', 'E')
names(h2A.OD)=colnames(residual.pheno.OD)
#save(vcA.OD.unscaled, file='~/Dropbox/Public/eQTL/vcA.unscaled.RData')

#### Calculate heritability for downsampled tpms ##################################
    downsampled.vcA=lapply(downsample.tpm.matrices, function(x) {
          tpm=x                     
          tpm=tpm[match(gene.GR$ORF, rownames(tpm)),]
          pheno.corrected=scale(residuals(lm(t(tpm)~covariates.OD)))
          vc.out=calcA(pheno.corrected, A)
          return(vc.out)
        } 
    )
    h2A.downsample=sapply(downsampled.vcA, function(x) x[,1]/c(x[,1]+x[,2]))

    g.counts=apply(count.matrix, 1, sum, na.rm=T)
    bins=cut(log10(as.vector(g.counts)),c(0,4,5,5.5,6,8))
    read.bins=rle(as.character(sort(bins)))
    read.bins=paste(read.bins$values, read.bins$lengths, sep=' n=')

    #png(file=paste0('/data/eQTL/plots/h2_vs_downsampling_by_abundance_class_new.png'), width=1920, height=1024)
    boxplot(h2A.OD~bins, border='black',at=c(((1:5)*5)-4), xlim=c(1,25), names=read.bins, xlab='observed log10(reads per transcript)', ylab='h^2')
    boxplot(h2A.downsample[,1]~bins, border='blue',at=c(((1:5)*5)-3), xlim=c(1,25), add=T,names=rep('',5))
    boxplot(h2A.downsample[,2]~bins, border='purple',at=c(((1:5)*5)-2), xlim=c(1,25), add=T,names=rep('',5))
    boxplot(h2A.downsample[,3]~bins, border='green',at=c(((1:5)*5)-1), xlim=c(1,25), add=T,names=rep('',5))
    boxplot(h2A.downsample[,4]~bins, border='red',at=c(((1:5)*5)), xlim=c(1,25), add=T,names=rep('',5))
    reads=c(mean(colSums(count.matrix)),1e6,5e5,1e5,5e4)
    legend('topleft', legend=round(log10(reads),2),
           title='downsampling (log10(reads per sample))',
           col=c('black', 'blue', 'purple', 'green', 'red'),
           fill=c('black', 'blue', 'purple', 'green', 'red'))
    #dev.off()

    plot(reads, c(mean(h2A.OD), apply(h2A.downsample,2, mean)), ylim=c(0,.5), xlab='reads per individual', ylab=expression(h^2), xlim=c(0,1e7), type='b')
    yin=as.numeric(c(mean(h2A.OD), apply(h2A.downsample,2, mean)))
    xin=reads
    nfit=nls(yin~d*xin/(b+xin))
    x2=seq(1e4, 20e6,length.out=100)
    predicty=((coef(nfit)['d']*x2))/(coef(nfit)['b']+x2)
    points(x2, predicty, col='blue', type='l')
    #points(x2, predict(loess(yin~xin, surface='direct', span=.2), x2), type='l', col='blue')
    sum(h2A.downsample[,1]>(h2A.OD*.9))
##################################################################################



# CIS-only model ------------------------------------------------------------------------------------------------------------------------
    #covariates=model.matrix(t.tpm.matrix[,1]~gbatch.fact)
    #residual.pheno=residuals(lm(t.tpm.matrix~covariates.OD) )
    #pheno.scaled=(scale(residual.pheno))
    #naive local model
    closest.marker.to.transcript=nearest(gene.GR, marker.GR)
    names(closest.marker.to.transcript)=colnames(t.tpm.matrix) 

    cisModel=rep(NA, 5718)
    cisModel.effect=rep(NA, 5718)

    for(i in 1:5718) {
        print(i)
        nm=lm(t.tpm.matrix[,i]~covariates.OD)
        fm=lm(t.tpm.matrix[,i]~covariates.OD+gdata[,closest.marker.to.transcript[i]])
        cisModel[i]=anova(nm,fm)$'Pr(>F)'[2]
        cisModel.effect[i]=as.numeric(coef(fm)[16])
    }
    pointwise.cis.results=data.frame(transcript=colnames(t.tpm.matrix)[1:5718], coefficient=cisModel.effect, p.value=cisModel, stringsAsFactors=F)
    #cme=(2^(cisModel.effects*2))
    #cme[cme<1]=1+(1-cme[cme<1])
    #hist(cme, xlim=c(1,4), breaks=100, xlab='(folded) fold difference between alleles', main='eQTL cis effect sizes',    sub=paste('median=1.109', '    IQR=.21') )

    #save(pointwise.cis.results, file='/data/eQTL/RData/pointwise.cis.test.RData')
    load('/data/eQTL/RData/pointwise.cis.test.RData')
    #(names(which(genes.with.cis==1)))%in%colnames(pheno.scaled)[((cisModel[1:6616]<(.05/6616)))]
    #2036
    #names(which(cisModel[1:6616]<(.05/6616))))
    sum(cisModel<(.05/5718)) #2460 
    sum(qvalue(cisModel)$qvalue<.05) #4241
#-----------------------------------------------------------------------------------------------------------------------------------------



## mapping QTL, one dimensional scan, one transcript at a time -----------------------------------------------------------------------------
background.QTL.OD = find.background.QTL(covariates.OD, t.tpm.matrix, pheno.scaled.OD, gdata, gdata.scaled)
peakList.OD = mapQTL(covariates.OD, background.QTL.OD,
                  t.tpm.matrix, pheno.scaled.OD, gdata, gdata.scaled,
                  n.perm=1000, FDR.thresh=.05)
##save(peakList.OD, file='/data/eQTL/RData/peakList_batchOD.RData')
#load('/data/eQTL/RData/peakList_batchOD.RData')
#-------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------------
# Total peak count
sum(sapply(sapply(peakList.OD, function(x) { sapply(x, function(y) nrow(y) ) } ), sum))
#36,498
# table of peaks
#36498
all.peaks.OD=buildPeakListDF(peakList.OD,gdata,gene.GR, marker.GR)
#save(all.peaks.OD, file='/data/eQTL/RData/all_peaks_OD.RData')
#load('/data/eQTL/RData/all_peaks_OD.RData')
load('/data/eQTL/RData/R_allPeaksODPad_161213.RData')
all.peaks.OD=allPeaksODPad[,-c(14,15)]

# if there are two CIS per transcript take the one closest to the gene
peaks.per.gene=split(all.peaks.OD, all.peaks.OD$gene)
peaks.per.gene=lapply(peaks.per.gene, function(x) {
       mdist=abs(x$marker.gcoord-x$gene.gcoord)
       wmm=which.min(mdist)
       if(x$cis[wmm]) {     x$cis[-wmm]=FALSE      }
       return(x) }
)
all.peaks.OD=do.call('rbind', peaks.per.gene)

# Big MAP ---------------------------------------------------------------------------
png('/data/eQTL/plots/eQTL_map_1920x1920.png', width=1920, height=1920)
eQTL_bigPlot(all.peaks.OD, gcoord.key,marker.GR)
dev.off()
#------------------------------------------------------------------------------------

peaks.per.gene=split(all.peaks.OD, all.peaks.OD$gene)
peaks.per.chr=split(all.peaks.OD, all.peaks.OD$chr)
genes.with.cis=sapply(peaks.per.gene, function(x) sum(x$cis))      #2969 .. filter 2884 genes with local  
genes.with.trans=sapply(peaks.per.gene, function(x) sum(!x$cis))


# mapping with cross validation
#doCV(gbatch,fact, covariates.OD, t.tpm.matrix, pheno.scaled.OD, gdata)

#load results
cvVE=loadCV(gbatch.fact, covariates.OD, t.tpm.matrix, gdata)

qCV=sapply(cvVE,mean)
h2As=h2A.OD[match(names(qCV), names(h2A.OD))]
qCVs=qCV[match(names(h2As), names(qCV))]
#sapply(cvVE,mean)
plot(h2As,qCVs, xlim=c(0,1), ylim=c(0,1))
h2Acv=sort(qCVs/h2As)
table(gbatch.fact)
median(h2Acv) #[h2Acv<1.5])
#save.image('/data/eQTL/RData/041917.RData')


# Remove additive effects
#pheno.additive.removed=removeAdditiveEffects(t.tpm.matrix, peaks.per.gene, covariates.OD, gdata)
    load('/data/eQTL/RData/pheno.additive.removed.RData')
    pheno.additive.removed.scaled=scale(pheno.additive.removed)

    pars=svd(pheno.additive.removed.scaled)
    pars.null=svd(apply(pheno.additive.removed.scaled, 2, function(x) x[sample(1:length(x))] ))

    plot(pars$d^2/sum(pars$d^2), xlim=c(0,100))
    points(pars.null$d^2/sum(pars.null$d^2), xlim=c(0,100), col='blue')
    # take top 20 PCs
    #pars$u[,1:20]
    background.covariates=cbind(covariates.OD, pars$u[,1:20])

# DO HOTSPOT ANALYSIS
#source('/data/eQTL/code/eQTL_Hotspots.R')

# BUILD B
#source('/data/eQTL/code/eQTL_BYxRM1000_makeB.R')


#DO 2-LOCUS AND ADDITIONAL VC ANALYSIS 
#source('/data/eQTL/code/eQTL_BYxRM1000_2Locus.R')
load('/data/eQTL/RData/peaksModel.RData')


# local eQTL median size
cil=as.numeric(sapply(strsplit(sapply(strsplit(all.peaks.OD$CI.l, ':'), function(x)x[2]), '_'), function(x)x[1]))
cir=as.numeric(sapply(strsplit(sapply(strsplit(all.peaks.OD$CI.r, ':'), function(x)x[2]), '_'), function(x)x[1]))
median(abs(cir-cil)[all.peaks.OD$cis])

# enrichment of local eQTL with interaction
#chisq.test(rbind(c(2884*36498,36498^2), c(300,1464))) 0.07902

#chisq.test(rbind(c(2884,36498), c(300,1464))) 0.07902
#data:  rbind(c(2884, 36498), c(300, 1464))    0.2049
#X-squared = 220, df = 1, p-value <2e-16


#siglev = pchisq(2.5*4.60517,1,lower.tail=FALSE)
siglev = pchisq(3.03*4.60517,1,lower.tail=FALSE)
siglev.int = pchisq(3.85*4.60517,1,lower.tail=FALSE)
dd     = seq(0.01,10,.01)
n1000 = power.t.test(n=506,  delta=seq(0.01,10,.01),sig.level=siglev)$power
#n1000.int = power.t.test(n=506,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
pvar= prop.var('ri', dd/2,1)
plot(pvar,n1000, type='l',col='blue',lwd=2)
#points(pvar,n1000.int, type='l',col='blue',lwd=2,lty=2)

dtrans=density(unlist(sapply(peaksModel, function(x) x$var.exp.Resid[!x$cis] )))
dcis=density(unlist(sapply(peaksModel, function(x) x$var.exp.Resid[x$cis] ) ))
#all.peaks.multiple.regression$var.exp[all.peaks.multiple.regression$cis])
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(dtrans, xlab='fraction of phenotypic variance', ylab='density', xlim=c(0,0.5),  type='n', main='')
polygon(dtrans, col=rgb(0, 0, 1,0.3), border=NA)
polygon(dcis, col=rgb(1, 0, 0,0.3), border=NA)
legend('topright', legend=c('distant', 'local'), fill=c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)))
sf= max(c(dtrans$y, dcis$y))
points((pvar), n1000*sf, type='l',col='black',lwd=3, lty=1, xlab=NA, ylab=NA)
#points((pvar), n1000.int*sf, type='l',col='red',lwd=3, lty=1,  xlab=NA, ylab=NA)
axis(side=4, at=seq(0,1, .2)*sf, labels= seq(0,1, .2))
mtext(side = 4, line = 3, "power")


#sum(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp))<.02)
#19603
#df=data.frame(variance.explained=pvar, power.additive.n=n1000)[1:100,]
#write.table(df[1:100,], file='~/Desktop/eQTL_power.txt', sep='\t', quote=F, row.names=F)
x=peakModel[[4]]
mNull=lm(pheno.scaled.OD~1)
mFull=lm(pheno.scaled.OD~gdata[,x$pmarker[x$cis]])


#save.image('/data/eQTL/052417.RData')
