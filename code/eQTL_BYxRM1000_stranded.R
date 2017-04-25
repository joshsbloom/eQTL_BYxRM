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
downsample.1e6.count.matrix=downsampleCounts(count.matrix, 1e6)
downsample.5e5.count.matrix=downsampleCounts(count.matrix, 5e5)
downsample.1e5.count.matrix=downsampleCounts(count.matrix, 1e5)
downsample.5e4.count.matrix=downsampleCounts(count.matrix, 5e4)

downsample.count.matrices=list(
                         'd1e6'=downsample.1e6.count.matrix,
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
#residual.pheno.OD=scale(residuals(lm(t(tpm.matrix)~covariates.OD) ))
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

    plot(reads, c(mean(h2A.OD), apply(h2A.downsample,2, mean)), ylim=c(0,.5), xlab='reads per individual', ylab=expression(h^2), xlim=c(0,20e6), type='b')
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
load('/data/eQTL/RData/peakList_batchOD.RData')
#-------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------------
# Total peak count
sum(sapply(sapply(peakList.OD, function(x) { sapply(x, function(y) nrow(y) ) } ), sum))
#36,498
# table of peaks
#36498
all.peaks.OD=buildPeakListDF(peakList.OD,gdata,gene.GR, marker.GR)
#save(all.peaks.OD, file='/data/eQTL/RData/all_peaks_OD.RData')
load('/data/eQTL/RData/all_peaks_OD.RData')
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

# Do HOTSPOT analysis
#source('/data/eQTL/code/eQTL_Hotspots.R')

# Build B
#source('/data/eQTL/code/eQTL_BYxRM1000_makeB.R')


# local eQTL median size
cil=as.numeric(sapply(strsplit(sapply(strsplit(all.peaks.OD$CI.l, ':'), function(x)x[2]), '_'), function(x)x[1]))
cir=as.numeric(sapply(strsplit(sapply(strsplit(all.peaks.OD$CI.r, ':'), function(x)x[2]), '_'), function(x)x[1]))
median(abs(cir-cil)[all.peaks.OD$cis])


# For 2D model that is only testing interactions between additive QTL ------------------
set.seed(1)
yr=residuals(lm(t.tpm.matrix~covariates.OD))
yr.perm=replicate(10, {yr[sample(1:1012),]})
y.for.2D=abind(yr, yr.perm, along=3)
lm.2Dint = foreach(p = 1:11) %dopar% {    
    print(p)
    intlist=list()
    pb =txtProgressBar(min = 1, max = length(peaks.per.gene), style = 3)
    for(g in 1:length(peaks.per.gene)) { 
        setTxtProgressBar(pb, g)
        gg=names(peaks.per.gene)[g]
        pms=peaks.per.gene[[g]]$pmarker[!duplicated(peaks.per.gene[[g]]$pmarker)]
        if(length(pms)>1) {
             X=data.frame(gdata[,pms])
             fitme=lm(y.for.2D[,gg,p]~., data=X)
             #lm.2Dint[[as.character(p)]][[gg]]=add1(fitme, ~.^2, test='F')
             intlist[[gg]]=add1(fitme, ~.^2, test='F')
        } 
    }
    close(pb)
    return(intlist)
}
#---------------------------------------------------------------------------------------
# FDR with q-value  ------------------------------------------
ps.QTL.2D=na.omit(unlist(sapply(lm.2Dint[[1]], function(x) x$'Pr(>F)')))
#qs2D =qvalue(ps.QTL.2D)
#qthresh2d=max(qs2D$pvalues[qs2D$qvalues<.1]) 
##R> qthresh2d
##[1] 0.001285
#---------------------------------------------------------------

ps.QTL.2D.perm=sapply(2:11, function(y) na.omit(unlist(sapply(lm.2Dint[[y]], function(x) x$'Pr(>F)'))))

obscnt.targetted=sapply(seq(1.5,8,.05) , function(x) sum( -log10(ps.QTL.2D) > x ) )
names(obscnt.targetted)=seq(1.5,8,.05)
permcnt.targetted=rowMeans(apply(ps.QTL.2D.perm, 2, function(y) { sapply(seq(1.5,8,.05) , function(x) sum( -log10(y) > x ) )     }))
names(permcnt.targetted)= seq(1.5,8,.05)

permcnt.targetted/obscnt.targetted
#FDR 5%  = 0.0003548    10^(3.45)
#FDR 10% = 0.001        10^(3)
qthresh2d=.001
sum(ps.QTL.2D<qthresh2d)

library(regress)

peaks.per.gene.augmented=list()
vc_multiple_components=list()
lm.2D=list()
lm.3D.test=list()


A=tcrossprod(gdata.scaled)/ncol(gdata.scaled)
AA=A*A
#Ahot=tcrossprod(gdata.scaled[,sort(match(as.character(htabler[,'peak']), colnames(gdata)))])/nrow(htabler)

set.seed(100)
registerDoMC(cores=10)
# calculate 
sample.sets=replicate(10, sample(1:nrow(t.tpm.matrix)) )
peaksModel=foreach(g = names(peaks.per.gene)) %dopar% {                   
    peaks_model2=list()
    print(g)
    yr=t.tpm.matrix[,g]
    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    peaks_model2$df=ppg

    apeaks = match(ppg$pmarker, colnames(gdata))

    X2=data.frame(gdata[,apeaks])
    X=data.frame(covariates.OD, gdata[,apeaks])
    #qA=tcrossprod(gdata.scaled[,apeaks])/length(apeaks)
    #qAxqA=qA*qA
    #qAxA=qA*A
    #rr=regress(yr~covariates.OD-1, ~qA+A+qAxqA+qAxA+AA, verbose=T, pos=c(T,T,T,T,T,T))
    #extract.rr(rr)

    fitme=lm(yr~.-1, data=X)
    aov.a = anova(fitme)
    tssq  = sum(aov.a[-c(1:14), 2])
    a.effs=aov.a[15:(nrow(aov.a)-1),2]/tssq # (aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(fitme) 

    peaks_model2$df$var.exp=a.effs
    peaks_model2$df$lm.coeff=as.vector(coeffs)[-c(1:14)]
    #peaks_model2$fit=fitme

    #if(length(apeaks)>1) {
    #     int.space= colnames(model.matrix(~.^2-1, X2))[-c(1:length(apeaks))]
    #     intps=add1(fitme, int.space, test='F')[-1,6]
    #     names(intps)=int.space
    #
    #     ypr=predict(fitme)+apply(sample.sets, 2, function(x) residuals(fitme)[x])
    #     intps.null=(as.matrix(apply(ypr,2, function(yp) {
    #            fitme=lm(yp~.-1, data=X)
    #            add1(fitme, int.space, test='F')[-1,6]
    #     }) ))
    #     if(ncol(intps.null)==1) {intps.null=t(intps.null)}
    #     rownames(intps.null)=int.space  
    #    peaks_model2$fit$intps=intps
    #    peaks_model2$fit$intps.null=intps.null
    #
    # }
    return(peaks_model2)
}
names(peaksModel)=names(peaks.per.gene)
peakModel=lapply(peaksModel, function(x) x$df)
#save(peakModel, file='/data/eQTL/RData/peakModel.RData')
load('/data/eQTL/RData/peakModel.RData')








# split this up into separate chunks 09 16 16
pb =txtProgressBar(min = 1, max = length(peaks.per.gene), style = 3)
for(gg in 1:length(peaks.per.gene)) {
    setTxtProgressBar(pb, gg)
    g=names(peaks.per.gene)[gg]
    #in names(peaks.per.gene) ) {
    #print(g)
    yr=pheno.scaled.OD[,g]
    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    apeaks = match( ppg$pmarker, colnames(gdata))

    pms=ppg$pmarker
    max.peak=ppg$pmarker[which.max(ppg$LOD)]
    qA=tcrossprod(gdata.scaled[,pms])/length(pms)
   
    rr=regress(yr~1, ~A+AA, verbose=F, pos=c(T,T,T,T,T,T))
    vc_multiple_components[[g]][['0']]=extract.rr(rr)

    #rr=regress(yr~1, ~qA+Ahot+A, verbose=F, pos=c(T,T,T,T,T,T,T))
    #vc_multiple_components[[g]][['1']]=extract.rr(rr)
    
    if(sum(ppg$cis)>0 & sum(!ppg$cis>0) ) {
        A.local=tcrossprod(gdata.scaled[,pms[ppg$cis]])/sum(ppg$cis)
        A.distant=tcrossprod(gdata.scaled[,pms[!ppg$cis]])/sum(!ppg$cis)
        A.local_X_A=A.local*A
        A.distant_X_A=A.distant*A
        rr=regress(yr~1, ~A.local + A.distant + A + A.local_X_A + A.distant_X_A + AA, verbose=T, pos=c(T,T,T,T,T,T,T))
        vc_multiple_components[[g]][['2']]=extract.rr(rr)
    }    
   
    rr=regress(yr~1, ~qA+A+AA, verbose=T, pos=c(T,T,T,T,T))
    vc_multiple_components[[g]][['3']]=extract.rr(rr)
    X=data.frame(gdata[,pms])
    fitme=lm(yr~.-1, data=X)
    #X2=data.frame(covariates.OD, gdata[,pms])
    #fitme2=lm(t.tpm.matrix[,g]~., data=X2)

    # for additive model ---------------------------------------------
    aov.a = anova(fitme)
    tssq  = sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(fitme)  
    peaks.per.gene.augmented[[g]]=ppg
    peaks.per.gene.augmented[[g]]$var.exp=a.effs
    peaks.per.gene.augmented[[g]]$lm.coeff=as.vector(coeffs)
    #----------------------------------------------------------------
    if(length(pms)>1) {
        a2=add1(fitme, ~.^2, test='F')
        nterms=rownames(a2[which(a2[,6]< qthresh2d),])
        if(length(nterms) > 0)  {
            fit2=update(fitme, as.formula(paste0('~.+',paste(nterms, collapse=' + '))))
            lm.2D[[g]]=fit2
            g2in=model.matrix(fit2)
            g2int=scale(g2in[,grep(':', colnames(g2in))])
            qAA=tcrossprod(g2int)/ncol(g2int)
            rr=  regress(yr~1, ~qA+A+qAA+AA, verbose=T, pos=c(T,T,T,T,T))
            vc_multiple_components[[g]][['4']]=extract.rr(rr)
            ints=rownames(a2)[-1]
            lm.3D.test[[g]]=add1(fit2,  paste0(gsub(':|/', '.', max.peak), ':', ints), test='F')
        }
    }
}
close(pb)
#save(eaks.per.gene.augmented,  file='/data/eQTL/RData/peaks_per_gene_augmented.RData')
#save(vc_multiple_components, file='/data/eQTL/RData/vc_QTL.RData')
#save(lm.2D, file='/data/eQTL/RData/lm.2D.RData')
#save(lm.3D.test, file= '/data/eQTL/RData/lm.3D.test.RData')

load('/data/eQTL/RData/peaks_per_gene_augmented.RData')
load('/data/eQTL/RData/vc_QTL.RData')
load('/data/eQTL/RData/lm.2D.RData')
load('/data/eQTL/RData/lm.3D.test.RData')

#median(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp)))
#[1] 0.01865
#mean(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp)))
#[1] 0.04


d=(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp)))
dl=(unlist(sapply(peaks.per.gene.augmented, function(x) x$LOD)))

  d$x[which.max(d$y)]

#R> mean(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[!x$cis] )))
#[1] 0.03417
#R> median(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[!x$cis] )))
#[1] 0.01786

#R> mean(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[x$cis] )))
#[1] 0.1167
#R> median(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[x$cis] )))
#[1] 0.05733


da=density(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[!x$cis] )))
di=density(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[x$cis] ) ))
#all.peaks.multiple.regression$var.exp[all.peaks.multiple.regression$cis])
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(da, xlab='fraction of phenotypic variance', ylab='density', xlim=c(0,0.5),  type='n', main='')
polygon(da, col=rgb(0, 0, 1,0.3), border=NA)
polygon(di, col=rgb(1, 0, 0,0.3), border=NA)
legend('topright', legend=c('distant', 'local'), fill=c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)))


































# new code to re-localize QTL
#tpm.q=apply(scale(t.tpm.matrix),2, function(x)  x[(x>quantile(x,.995))]  ) 

# remap QTL hotspots
# HOTSPOT ANALYSIS --------------------------------------------------------------------------------------------------
peaks.per.chr=split(all.peaks.OD, all.peaks.OD$chr)

gbatch.fact.OD=covariates.OD
#data.frame(gbatch.fact, OD.cov)

genetic.map=buildGeneticMap(gdata.by.chr)

dir.create('/data/eQTL/RData/hotspots/')
#add directory as an argument to this function
hotspots.OD=findHotspots.iteration1(t.tpm.matrix, gene.annot.df, 
                                 peaks.per.chr, genetic.map, gdata, gdata.s.by.chr, gbatch.fact.OD) 
#hotspots=findHotspots.iteration1(t.tpm.matrix, gene.annot.df, 
#                                 peaks.per.chr, genetic.map, gdata, gdata.s.by.chr, gbatch.fact) 
#save(hotspots.OD, file='/data/eQTL/RData/hotspots.RData')

load('/data/eQTL/RData/hotspots.RData')
hp.index=sort(match(unlist(do.call('rbind', lapply(hotspots.OD, function(x) do.call('rbind', x)))[,3]), colnames(gdata)))

hotspot.boots=list()
for( i in 1:86) { #length(hp.index)) {    
    print(i)
    hotspot.boots[[i]]=bootstrap.Hotspots(i)  
}

# aggregate information from hotspot analysis 

hotspot.boot.peaks=list()
hotspot.boot.intervals=matrix(NA, 86,5)
for(nn in c(1:86)) {
    # will load 'peaks'
    load(paste0('/data/eQTL/RData/hotspots/', nn, '.bootstrap.peaks'))
    r=gsub(':', '_', peaks)
    pos=as.numeric(sapply(strsplit(r, '_'), function(x)x[2]))
    hotspot.boot.peaks[[nn]]=pos
    hotspot.boot.intervals[nn,]=c(quantile(pos, .025), quantile(pos, .05),  quantile(pos, .5),  quantile(pos, .95),quantile(pos, .975))
}

# how hot to be a hotspot?
bins=seq(1,1.2e7,7500)
hp.null=cut(all.peaks.OD$marker.gcoord, bins)
lambda.pois=length(hp.null)/length(bins)
qpois(.05/length(bins), lambda.pois, lower.tail=F)

write.table(
            data.frame(cbind(as.character(do.call('rbind', 
                               lapply(hotspots.OD, function(x) do.call('rbind',x)))[,3]),hotspot.boot.intervals)),
            file='/data/eQTL/RData/hotspots_081016.txt', sep='\t' ,quote=F, row.names=FALSE)

hotspot.positions=as.character(do.call('rbind', lapply(hotspots.OD, function(x) do.call('rbind',x)))[,3])
hotspot.list=list(hotspot.positions=hotspot.positions,
                  hotspot.boot.peaks=hotspot.boot.peaks,
                  hotspot.boot.intervals=hotspot.boot.intervals)
#save(hotspot.list, file='/data/eQTL/RData/hotspot_with_boots_list.RData')
load('/data/eQTL/RData/hotspot_with_boots_list.RData')

#save.image('/data/eQTL/RData/090116.RData')
load('/data/eQTL/RData/090116.RData')
#load('/home/jbloom/Dropbox/Public/eQTL/hotspot_with_boots_list_062316.RData')
#--------------------------------------------------------------------------------------------------------------------

#



#here's a logical vector. TRUE => remove cis effects, FALSE => do not 
#remove cis effect

#there are 1694 FALSE, i.e. "protected" cis effects. This is a lot 
#because they include all genes that 1) physically overlap with a 95% 
#hotspot interval, and 2) that have a cis eQTL that itself overlaps with 
#a hotspot, even if the gene itself is outside of the hotspot. Thought 
#for the latter group is that these might be genes right next to a 
#hotspot where localization uncertainty might lead us astray. "cis" peak 
#is defined using padded genes (-1000, +200) and extended perfect LD eQTL 
#markers.

#So this is very generous in terms of how many genes it protects. Should 
#ensure that any cis-mediated hotspots do not get corrected away. Clearly 

obscnt.targetted=sapply(seq(1.5,8,.05) , function(x) sum( -log10(ps.QTL.2D) > x ) )
names(obscnt.targetted)=seq(1.5,8,.05)
permcnt.targetted=rowMeans(apply(ps.QTL.2D.perm, 2, function(y) { sapply(seq(1.5,8,.05) , function(x) sum( -log10(y) > x ) )     }))
names(permcnt.targetted)= seq(1.5,8,.05)

permcnt.targetted/obscnt.targetted
#FDR 5%  = 0.0003548    10^(3.45)
#FDR 10% = 0.001        10^(3)
qthresh2d=.001
sum(ps.QTL.2D<qthresh2d)

library(regress)

peaks.per.gene.augmented=list()
vc_multiple_components=list()
lm.2D=list()
lm.3D.test=list()


A=tcrossprod(gdata.scaled)/ncol(gdata.scaled)
AA=A*A
#Ahot=tcrossprod(gdata.scaled[,sort(match(as.character(htabler[,'peak']), colnames(gdata)))])/nrow(htabler)

set.seed(100)
registerDoMC(cores=70)
sample.sets=replicate(10, sample(1:nrow(t.tpm.matrix)) )
peaksModel=foreach(g = names(peaks.per.gene)) %dopar% {                   
    peaks_model2=list()
    print(g)
    yr=t.tpm.matrix[,g]
    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    peaks_model2$df=ppg

    apeaks = match(ppg$pmarker, colnames(gdata))

    X2=data.frame(gdata[,apeaks])
    X=data.frame(covariates.OD, gdata[,apeaks])
    #qA=tcrossprod(gdata.scaled[,apeaks])/length(apeaks)
    #qAxqA=qA*qA
    #qAxA=qA*A
    #rr=regress(yr~covariates.OD-1, ~qA+A+qAxqA+qAxA+AA, verbose=T, pos=c(T,T,T,T,T,T))
    #extract.rr(rr)

    fitme=lm(yr~.-1, data=X)
    aov.a = anova(fitme)
    tssq  = sum(aov.a[-c(1:14), 2])
    a.effs=aov.a[15:(nrow(aov.a)-1),2]/tssq # (aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(fitme) 

    peaks_model2$df$var.exp=a.effs
    peaks_model2$df$lm.coeff=as.vector(coeffs)[-c(1:14)]

    peaks_model2$fit=fitme
     if(length(apeaks)>1) {
         int.space= colnames(model.matrix(~.^2-1, X2))[-c(1:length(apeaks))]
         intps=add1(fitme, int.space, test='F')[-1,6]
         names(intps)=int.space

         ypr=predict(fitme)+apply(sample.sets, 2, function(x) residuals(fitme)[x])
         intps.null=(as.matrix(apply(ypr,2, function(yp) {
                fitme=lm(yp~.-1, data=X)
                add1(fitme, int.space, test='F')[-1,6]
         }) ))
         if(ncol(intps.null)==1) {intps.null=t(intps.null)}
         rownames(intps.null)=int.space  
        peaks_model2$fit$intps=intps
        peaks_model2$fit$intps.null=intps.null

     }
       return(peaks_model2)
}
names(peaksModel)=names(peaks.per.gene)
peakModel=lapply(peaksModel, function(x) x$df)
#save(peakModel, file='/data/eQTL/RData/peakModel.RData')
load('/data/eQTL/RData/peakModel.RData')



















# split this up into separate chunks 09 16 16
pb =txtProgressBar(min = 1, max = length(peaks.per.gene), style = 3)
for(gg in 1:length(peaks.per.gene)) {
    setTxtProgressBar(pb, gg)
    g=names(peaks.per.gene)[gg]
    #in names(peaks.per.gene) ) {
    #print(g)
    yr=pheno.scaled.OD[,g]
    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    apeaks = match( ppg$pmarker, colnames(gdata))

    pms=ppg$pmarker
    max.peak=ppg$pmarker[which.max(ppg$LOD)]
    qA=tcrossprod(gdata.scaled[,pms])/length(pms)
   
    rr=regress(yr~1, ~A+AA, verbose=F, pos=c(T,T,T,T,T,T))
    vc_multiple_components[[g]][['0']]=extract.rr(rr)

    #rr=regress(yr~1, ~qA+Ahot+A, verbose=F, pos=c(T,T,T,T,T,T,T))
    #vc_multiple_components[[g]][['1']]=extract.rr(rr)
    
    if(sum(ppg$cis)>0 & sum(!ppg$cis>0) ) {
        A.local=tcrossprod(gdata.scaled[,pms[ppg$cis]])/sum(ppg$cis)
        A.distant=tcrossprod(gdata.scaled[,pms[!ppg$cis]])/sum(!ppg$cis)
        A.local_X_A=A.local*A
        A.distant_X_A=A.distant*A
        rr=regress(yr~1, ~A.local + A.distant + A + A.local_X_A + A.distant_X_A + AA, verbose=T, pos=c(T,T,T,T,T,T,T))
        vc_multiple_components[[g]][['2']]=extract.rr(rr)
    }    
   
    rr=regress(yr~1, ~qA+A+AA, verbose=T, pos=c(T,T,T,T,T))
    vc_multiple_components[[g]][['3']]=extract.rr(rr)
    X=data.frame(gdata[,pms])
    fitme=lm(yr~.-1, data=X)
    #X2=data.frame(covariates.OD, gdata[,pms])
    #fitme2=lm(t.tpm.matrix[,g]~., data=X2)

    # for additive model ---------------------------------------------
    aov.a = anova(fitme)
    tssq  = sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(fitme)  
    peaks.per.gene.augmented[[g]]=ppg
    peaks.per.gene.augmented[[g]]$var.exp=a.effs
    peaks.per.gene.augmented[[g]]$lm.coeff=as.vector(coeffs)
    #----------------------------------------------------------------
    if(length(pms)>1) {
        a2=add1(fitme, ~.^2, test='F')
        nterms=rownames(a2[which(a2[,6]< qthresh2d),])
        if(length(nterms) > 0)  {
            fit2=update(fitme, as.formula(paste0('~.+',paste(nterms, collapse=' + '))))
            lm.2D[[g]]=fit2
            g2in=model.matrix(fit2)
            g2int=scale(g2in[,grep(':', colnames(g2in))])
            qAA=tcrossprod(g2int)/ncol(g2int)
            rr=  regress(yr~1, ~qA+A+qAA+AA, verbose=T, pos=c(T,T,T,T,T))
            vc_multiple_components[[g]][['4']]=extract.rr(rr)
            ints=rownames(a2)[-1]
            lm.3D.test[[g]]=add1(fit2,  paste0(gsub(':|/', '.', max.peak), ':', ints), test='F')
        }
    }
}
close(pb)
#save(eaks.per.gene.augmented,  file='/data/eQTL/RData/peaks_per_gene_augmented.RData')
#save(vc_multiple_components, file='/data/eQTL/RData/vc_QTL.RData')
#save(lm.2D, file='/data/eQTL/RData/lm.2D.RData')
#save(lm.3D.test, file= '/data/eQTL/RData/lm.3D.test.RData')

load('/data/eQTL/RData/peaks_per_gene_augmented.RData')
load('/data/eQTL/RData/vc_QTL.RData')
load('/data/eQTL/RData/lm.2D.RData')
load('/data/eQTL/RData/lm.3D.test.RData')

#median(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp)))
#[1] 0.01865
#mean(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp)))
#[1] 0.04


d=(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp)))
dl=(unlist(sapply(peaks.per.gene.augmented, function(x) x$LOD)))

  d$x[which.max(d$y)]

#R> mean(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[!x$cis] )))
#[1] 0.03417
#R> median(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[!x$cis] )))
#[1] 0.01786

#R> mean(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[x$cis] )))
#[1] 0.1167
#R> median(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[x$cis] )))
#[1] 0.05733


da=density(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[!x$cis] )))
di=density(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[x$cis] ) ))
#all.peaks.multiple.regression$var.exp[all.peaks.multiple.regression$cis])
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(da, xlab='fraction of phenotypic variance', ylab='density', xlim=c(0,0.5),  type='n', main='')
polygon(da, col=rgb(0, 0, 1,0.3), border=NA)
polygon(di, col=rgb(1, 0, 0,0.3), border=NA)
legend('topright', legend=c('distant', 'local'), fill=c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)))


































# new code to re-localize QTL
#tpm.q=apply(scale(t.tpm.matrix),2, function(x)  x[(x>quantile(x,.995))]  ) 

# remap QTL hotspots
# HOTSPOT ANALYSIS --------------------------------------------------------------------------------------------------
peaks.per.chr=split(all.peaks.OD, all.peaks.OD$chr)

gbatch.fact.OD=covariates.OD
#data.frame(gbatch.fact, OD.cov)

genetic.map=buildGeneticMap(gdata.by.chr)

dir.create('/data/eQTL/RData/hotspots/')
#add directory as an argument to this function
hotspots.OD=findHotspots.iteration1(t.tpm.matrix, gene.annot.df, 
                                 peaks.per.chr, genetic.map, gdata, gdata.s.by.chr, gbatch.fact.OD) 
#hotspots=findHotspots.iteration1(t.tpm.matrix, gene.annot.df, 
#                                 peaks.per.chr, genetic.map, gdata, gdata.s.by.chr, gbatch.fact) 
#save(hotspots.OD, file='/data/eQTL/RData/hotspots.RData')

load('/data/eQTL/RData/hotspots.RData')
hp.index=sort(match(unlist(do.call('rbind', lapply(hotspots.OD, function(x) do.call('rbind', x)))[,3]), colnames(gdata)))

hotspot.boots=list()
for( i in 1:86) { #length(hp.index)) {    
    print(i)
    hotspot.boots[[i]]=bootstrap.Hotspots(i)  
}

# aggregate information from hotspot analysis 

hotspot.boot.peaks=list()
hotspot.boot.intervals=matrix(NA, 86,5)
for(nn in c(1:86)) {
    # will load 'peaks'
    load(paste0('/data/eQTL/RData/hotspots/', nn, '.bootstrap.peaks'))
    r=gsub(':', '_', peaks)
    pos=as.numeric(sapply(strsplit(r, '_'), function(x)x[2]))
    hotspot.boot.peaks[[nn]]=pos
    hotspot.boot.intervals[nn,]=c(quantile(pos, .025), quantile(pos, .05),  quantile(pos, .5),  quantile(pos, .95),quantile(pos, .975))
}

# how hot to be a hotspot?
bins=seq(1,1.2e7,7500)
hp.null=cut(all.peaks.OD$marker.gcoord, bins)
lambda.pois=length(hp.null)/length(bins)
qpois(.05/length(bins), lambda.pois, lower.tail=F)

write.table(
            data.frame(cbind(as.character(do.call('rbind', 
                               lapply(hotspots.OD, function(x) do.call('rbind',x)))[,3]),hotspot.boot.intervals)),
            file='/data/eQTL/RData/hotspots_081016.txt', sep='\t' ,quote=F, row.names=FALSE)

hotspot.positions=as.character(do.call('rbind', lapply(hotspots.OD, function(x) do.call('rbind',x)))[,3])
hotspot.list=list(hotspot.positions=hotspot.positions,
                  hotspot.boot.peaks=hotspot.boot.peaks,
                  hotspot.boot.intervals=hotspot.boot.intervals)
#save(hotspot.list, file='/data/eQTL/RData/hotspot_with_boots_list.RData')
load('/data/eQTL/RData/hotspot_with_boots_list.RData')

#save.image('/data/eQTL/RData/090116.RData')
load('/data/eQTL/RData/090116.RData')
#load('/home/jbloom/Dropbox/Public/eQTL/hotspot_with_boots_list_062316.RData')
#--------------------------------------------------------------------------------------------------------------------

#



#here's a logical vector. TRUE => remove cis effects, FALSE => do not 
#remove cis effect

#there are 1694 FALSE, i.e. "protected" cis effects. This is a lot 
#because they include all genes that 1) physically overlap with a 95% 
#hotspot interval, and 2) that have a cis eQTL that itself overlaps with 
#a hotspot, even if the gene itself is outside of the hotspot. Thought 
#for the latter group is that these might be genes right next to a 
#hotspot where localization uncertainty might lead us astray. "cis" peak 
#is defined using padded genes (-1000, +200) and extended perfect LD eQTL 
#markers.

#So this is very generous in terms of how many genes it protects. Should 
#ensure that any cis-mediated hotspots do not get corrected away. Clearly 
#there are choices here, but let's run with this for now. We have the 
#fully cis-corrected version in the bag after all.

load('/data/eQTL/RData/R_removeCis_withPaddedCisGeneDefinitions_160901.RData')
# removeCIS yes = remove_cis

# regress out local effect
bc.resid=matrix(NA, 1012, ncol(t.tpm.matrix))
for(i in 1:ncol(t.tpm.matrix)) {
    print(i)
    if(removeCis[i]) {
        bc.resid[,i]=residuals(lm(t.tpm.matrix[,i]~covariates.OD+gdata[,closest.marker.to.transcript[i]]))
    }else{
        bc.resid[,i]=residuals(lm(t.tpm.matrix[,i]~covariates.OD))
    }
}
bc.resid=scale(bc.resid)
colnames(bc.resid)=colnames(t.tpm.matrix)[1:ncol(t.tpm.matrix)]
save(bc.resid, file='/data/eQTL/RData/bcresid.RData')

#ols.hotspot
ols.hotspots=list()
for(i in 1:ncol(bc.resid)) {
    print(i)
    ols.hotspots[[colnames(bc.resid)[i]]]=
                coef(lm(scale(t.tpm.matrix[,i])~covariates.OD+gdata[,closest.marker.to.transcript[i]]+gdata[,hp.index]-1))
                                               #tcv.EBglmnet(gdata.scaled[,as.character(htabler$peak)], bc.resid[,i])
}
ols.matrix=do.call('rbind', ols.hotspots)

da=density(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[!x$cis] )))
di=density(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp[x$cis] ) ))
#all.peaks.multiple.regression$var.exp[all.peaks.multiple.regression$cis])
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(da, xlab='fraction of phenotypic variance', ylab='density', xlim=c(0,0.5),  type='n', main='')
polygon(da, col=rgb(0, 0, 1,0.3), border=NA)
polygon(di, col=rgb(1, 0, 0,0.3), border=NA)
legend('topright', legend=c('distant', 'local'), fill=c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)))

names(h2A)=colnames(t.tpm.matrix)


vc3=sapply(vc_multiple_components, function(x) x[['3']]$sigma)
plot(vc3[1,]+vc3[2,], vc3[1,])

#mean(vc3[1,]/(vc3[1,]+vc3[2,]))
# 0.7814
#median(vc3[1,]/(vc3[1,]+vc3[2,]))
#[1] 0.8041

vc6=do.call('rbind', sapply(vc_multiple_components, function(x) x[['2']]$sigma))
plot(c(1:6), vc6[1,-7], ylim=c(0,1))

for(i in 2:2576){
points(c(1:6), vc6[i,-7], col='#00000011', type='b', lty=2)
}

vchot=sapply(vc_multiple_components, function(x) x[['1']]$sigma)

mean(c(vchot[1,]+vchot[2,])/c(vchot[1,]+vchot[2,]+vchot[3,]))

plot(vc3[1,]+vc3[2,], vc3[1,])

lm.2D.count=sapply(lm.2D, function(x) {
#   print(length(genes.near.hotspot))
#   print(length(sig.at.hotspot))
#   cis.mediation.LOD=array(NA, c(length(sig.at.hotspot), length(genes.near.hotspot), 9) )
#   dimnames(cis.mediation.LOD)[[1]]=names(sig.at.hotspot)
#   dimnames(cis.mediation.LOD)[[2]]=genes.near.hotspot$ORF
#   dimnames(cis.mediation.LOD)[[3]]=c('L0', 'L1', 'L2', 'L3', 'L32', 'L10', 'L12', 'L13', 'mLOD')
#   for(cg in 1:length(genes.near.hotspot)) {
#        #print(cg)
#        cisg=genes.near.hotspot[cg]$ORF
#        
#        sya=t.tpm.matrix[,sig.at.hotspot]
#        rmcis.ind1=match(cisg, colnames(sya))
#        if(!is.na(rmcis.ind1)) {
#            sy=sya[,-rmcis.ind1]
#        } else {sy=sya} 
#
#        rmcis.ind=match(cisg, colnames(t.tpm.matrix))
#        cis.e=(t.tpm.matrix[,rmcis.ind])
#        #testing if growth (OD)  mediates QTL effects on expression levels  -----------------------
#        #sy=scale(t.tpm.matrix[,gind])
#        #sOD=scale(OD.cov)
#        
#save(ols.matrix, file='/data/eQTL/RData/ols.matrix.hotspots.RData')
load('/data/eQTL/RData/ols.matrix.hotspots.RData')

lasso.hotspots=list()
for(i in 1:ncol(bc.resid)) {
    print(i)
    lasso.hotspots[[colnames(bc.resid)[i]]]=cv.EBglmnet(gdata.scaled[,hp.index], bc.resid[,i])
}

lasso.matrix=matrix(0, length(hp.index), ncol(bc.resid) ) 
rownames(lasso.matrix)=as.character(colnames(gdata)[hp.index]) #htabler$peak) #colnames(gdata)[as.character(htabler$peak)]
colnames(lasso.matrix)=colnames(bc.resid)

for(i in 1:ncol(bc.resid)) {
    print(i)
    lh=lasso.hotspots[[colnames(bc.resid)[i]]]
    lasso.matrix[lh$fit[,1],i]=lh$fit[,3]
}
#save(lasso.matrix, file='/data/eQTL/RData/lasso.matrix.hotspots.selectiveCisRemoval.RData')
load('/data/eQTL/RData/lasso.matrix.hotspots.selectiveCisRemoval.RData')

#save(lasso.matrix, file='/data/eQTL/RData/lasso.matrix.hotspots.0630116.RData')
#save(lasso.matrix, file='/data/eQTL/RData/lasso.matrix.hotspots.072216.selectiveCisRemoval.RData')
#save(lasso.matrix, file='/home/jbloom/Dropbox/Public/eQTL/lasso.matrix.hotspots.072216.selectiveCisRemoval.RData')


htabler=do.call('rbind', lapply(hotspots.OD, function(x) do.call('rbind', x)))
#htabler=do.call('rbind', lapply(hotspots.refined, function(x) do.call('rbind', x)))

#abline(v=marker.GR$gcoord[match(htable$peak, marker.GR$mname)])
#abline(v=marker.GR$gcoord[match(htabler$peak, marker.GR$mname)])
#abline(v=match(sapply(hotspots[[cc]], function(x) x[,'peak']), colnames(gdata)))

#naive local model
cisModel=rep(NA, 5718)
for(i in 1:5718) {
    print(i)
    cisModel[i]=anova(lm(t.tpm.matrix[,i]~covariates.OD)
                      ,lm(t.tpm.matrix[,i]~covariates.OD+gdata[,closest.marker.to.transcript[i]]) )$'Pr(>F)'[2]
}
cisModel.effects=rep(NA, 5718)
for(i in 1:5718) {
    print(i)
    cisModel.effects[i]=
    as.numeric(coefficients(lm((t.tpm.matrix[,i])~covariates.OD+gdata[,closest.marker.to.transcript[i]]-1) ))[15]
}
cme=(2^(cisModel.effects*2))
cme[cme<1]=1+(1-cme[cme<1])
hist(cme, xlim=c(1,4), breaks=100, xlab='(folded) fold difference between alleles', main='eQTL cis effect sizes',
     sub=paste('median=1.109', '    IQR=.21') )
#save(cisModel.effects, file='~/Desktop/cisEffect.RData')

which(as.character(seqnames(marker.GR)) %in% unique.chrs[1])

eQTL_bigPlot(all.peaks, gcoord.key, marker.GR, which(as.character(seqnames(marker.GR)) %in% unique.chrs[1]))
x11()
par(xaxs='i')
plot(gene.GR$gcoord[1:6288], smooth(abs(cisModel.effects),twiceit=T))

vc_multiple_components=list()
lm.2Dupdate=list()
lm.3Dint2=list()
#[3415:length(names(peaks.per.gene))]


# test for 2D interactions at peak marker locations from 1D scan
######################################################

names(h2A)=colnames(t.tpm.matrix)

library(qtlDesign)
#siglev = pchisq(2.5*4.60517,1,lower.tail=FALSE)
siglev = pchisq(3.03*4.60517,1,lower.tail=FALSE)

siglev.int = pchisq(3.85*4.60517,1,lower.tail=FALSE)
dd     = seq(0.01,10,.01)
n1000 = power.t.test(n=506,  delta=seq(0.01,10,.01),sig.level=siglev)$power
#n1000.int = power.t.test(n=506,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
pvar= prop.var('ri', dd/2,1)
plot(pvar,n1000, type='l',col='blue',lwd=2)
points(pvar,n1000.int, type='l',col='blue',lwd=2,lty=2)

sum(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp))<.02)
19603


df=data.frame(variance.explained=pvar, power.additive.n=n1000)[1:100,]
write.table(df[1:100,], file='~/Desktop/eQTL_power.txt', sep='\t', quote=F, row.names=F)


, power.interactions.n=n1000.int)






#85% power
approxfun(n1000,pvar)(.85)

vc0=sapply(vc_multiple_components, function(x) x[['0']]$sigma)
plot(vc3[1,]+vc3[2,], vc3[1,])

#mean(vc3[1,]/(vc3[1,]+vc3[2,]))
# 0.7814
#median(vc3[1,]/(vc3[1,]+vc3[2,]))
#[1] 0.8041

#local.mediation=foreach(h=1:length(hotspot.peaks)) %dopar%  {
#   hp=hotspot.peaks[h]
#   print(hp)
#   mint=marker.GR[marker.GR$mname==hp]
#   # all genes within 10 kb of peak marker
#   #make this informed by the intervals ###
#   genes.near.hotspot.int=subsetByOverlaps(gene.GR, mint+10000)
#   genes.near.hotspot95=gene.GR[na.omit(match(as.vector(genesInHotspotList[[h]][[1]]), gene.GR$ORF)),]
#    
#   if(length(genes.near.hotspot95)> length(genes.near.hotspot.int) ) { 
#       genes.near.hotspot=genes.near.hotspot95 } else{
#       genes.near.hotspot=genes.near.hotspot.int       
#   }
#   null= lm(t.tpm.matrix~covariates.OD)
#   full= lm(t.tpm.matrix~covariates.OD+gdata[,hp])
#   rn=residuals(null)
#   fn=residuals(full) 
#   rss1=colSums(rn^2)
#   rss2=colSums(fn^2)
#   Fstat=(rss1-rss2)/(rss2/(1012-14))
#   hotspot.model=pf(Fstat, 1,1012-14, lower.tail=F)
#   qhot=qvalue(hotspot.model)
#   sig.at.hotspot=which(qhot$qvalue<.05)
#   print(length(genes.near.hotspot))
#   print(length(sig.at.hotspot))
#   cis.mediation.LOD=array(NA, c(length(sig.at.hotspot), length(genes.near.hotspot), 9) )
#   dimnames(cis.mediation.LOD)[[1]]=names(sig.at.hotspot)
#   dimnames(cis.mediation.LOD)[[2]]=genes.near.hotspot$ORF
#   dimnames(cis.mediation.LOD)[[3]]=c('L0', 'L1', 'L2', 'L3', 'L32', 'L10', 'L12', 'L13', 'mLOD')
#   for(cg in 1:length(genes.near.hotspot)) {
#        #print(cg)
#        cisg=genes.near.hotspot[cg]$ORF
#        
#        sya=t.tpm.matrix[,sig.at.hotspot]
#        rmcis.ind1=match(cisg, colnames(sya))
#        if(!is.na(rmcis.ind1)) {
#            sy=sya[,-rmcis.ind1]
#        } else {sy=sya} 
#
#        rmcis.ind=match(cisg, colnames(t.tpm.matrix))
#        cis.e=(t.tpm.matrix[,rmcis.ind])
#        #testing if growth (OD)  mediates QTL effects on expression levels  -----------------------
#        #sy=scale(t.tpm.matrix[,gind])
#        #sOD=scale(OD.cov)
#        
#save(ols.matrix, file='/data/eQTL/RData/ols.matrix.hotspots.RData')
load('/data/eQTL/RData/ols.matrix.hotspots.RData')

lasso.hotspots=list()
for(i in 1:ncol(bc.resid)) {
    print(i)
    lasso.hotspots[[colnames(bc.resid)[i]]]=cv.EBglmnet(gdata.scaled[,hp.index], bc.resid[,i])
}

lasso.matrix=matrix(0, length(hp.index), ncol(bc.resid) ) 
rownames(lasso.matrix)=as.character(colnames(gdata)[hp.index]) #htabler$peak) #colnames(gdata)[as.character(htabler$peak)]
colnames(lasso.matrix)=colnames(bc.resid)

for(i in 1:ncol(bc.resid)) {
    print(i)
    lh=lasso.hotspots[[colnames(bc.resid)[i]]]
    lasso.matrix[lh$fit[,1],i]=lh$fit[,3]
}
#save(lasso.matrix, file='/data/eQTL/RData/lasso.matrix.hotspots.selectiveCisRemoval.RData')
load('/data/eQTL/RData/lasso.matrix.hotspots.selectiveCisRemoval.RData')

#save(lasso.matrix, file='/data/eQTL/RData/lasso.matrix.hotspots.0630116.RData')
#save(lasso.matrix, file='/data/eQTL/RData/lasso.matrix.hotspots.072216.selectiveCisRemoval.RData')
#save(lasso.matrix, file='/home/jbloom/Dropbox/Public/eQTL/lasso.matrix.hotspots.072216.selectiveCisRemoval.RData')


htabler=do.call('rbind', lapply(hotspots.OD, function(x) do.call('rbind', x)))
#htabler=do.call('rbind', lapply(hotspots.refined, function(x) do.call('rbind', x)))

#abline(v=marker.GR$gcoord[match(htable$peak, marker.GR$mname)])
#abline(v=marker.GR$gcoord[match(htabler$peak, marker.GR$mname)])
#abline(v=match(sapply(hotspots[[cc]], function(x) x[,'peak']), colnames(gdata)))

#naive local model
cisModel=rep(NA, 5718)
for(i in 1:5718) {
    print(i)
    cisModel[i]=anova(lm(t.tpm.matrix[,i]~covariates.OD)
                      ,lm(t.tpm.matrix[,i]~covariates.OD+gdata[,closest.marker.to.transcript[i]]) )$'Pr(>F)'[2]
}
cisModel.effects=rep(NA, 5718)
for(i in 1:5718) {
    print(i)
    cisModel.effects[i]=
    as.numeric(coefficients(lm((t.tpm.matrix[,i])~covariates.OD+gdata[,closest.marker.to.transcript[i]]-1) ))[15]
}
cme=(2^(cisModel.effects*2))
cme[cme<1]=1+(1-cme[cme<1])
hist(cme, xlim=c(1,4), breaks=100, xlab='(folded) fold difference between alleles', main='eQTL cis effect sizes',
     sub=paste('median=1.109', '    IQR=.21') )
#save(cisModel.effects, file='~/Desktop/cisEffect.RData')

which(as.character(seqnames(marker.GR)) %in% unique.chrs[1])

eQTL_bigPlot(all.peaks, gcoord.key, marker.GR, which(as.character(seqnames(marker.GR)) %in% unique.chrs[1]))
x11()
par(xaxs='i')
plot(gene.GR$gcoord[1:6288], smooth(abs(cisModel.effects),twiceit=T))

vc_multiple_components=list()
lm.2Dupdate=list()
lm.3Dint2=list()
#[3415:length(names(peaks.per.gene))]


# test for 2D interactions at peak marker locations from 1D scan
######################################################

names(h2A)=colnames(t.tpm.matrix)

library(qtlDesign)
#siglev = pchisq(2.5*4.60517,1,lower.tail=FALSE)
siglev = pchisq(3.03*4.60517,1,lower.tail=FALSE)

siglev.int = pchisq(3.85*4.60517,1,lower.tail=FALSE)
dd     = seq(0.01,10,.01)
n1000 = power.t.test(n=506,  delta=seq(0.01,10,.01),sig.level=siglev)$power
#n1000.int = power.t.test(n=506,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
pvar= prop.var('ri', dd/2,1)
plot(pvar,n1000, type='l',col='blue',lwd=2)
points(pvar,n1000.int, type='l',col='blue',lwd=2,lty=2)

sum(unlist(sapply(peaks.per.gene.augmented, function(x) x$var.exp))<.02)
19603


df=data.frame(variance.explained=pvar, power.additive.n=n1000)[1:100,]
write.table(df[1:100,], file='~/Desktop/eQTL_power.txt', sep='\t', quote=F, row.names=F)


, power.interactions.n=n1000.int)






#85% power
approxfun(n1000,pvar)(.85)

vc0=sapply(vc_multiple_components, function(x) x[['0']]$sigma)
plot(vc3[1,]+vc3[2,], vc3[1,])

#mean(vc3[1,]/(vc3[1,]+vc3[2,]))
# 0.7814
#median(vc3[1,]/(vc3[1,]+vc3[2,]))
#[1] 0.8041

vc6=do.call('rbind', sapply(vc_multiple_components, function(x) x[['2']]$sigma))
plot(c(1:6), vc6[1,-7], ylim=c(0,1))

for(i in 2:2576){
points(c(1:6), vc6[i,-7], col='#00000011', type='b', lty=2)
}

vchot=sapply(vc_multiple_components, function(x) x[['1']]$sigma)

mean(c(vchot[1,]+vchot[2,])/c(vchot[1,]+vchot[2,]+vchot[3,]))

plot(vc3[1,]+vc3[2,], vc3[1,])

lm.2D.count=sapply(lm.2D, function(x) {
    sum(grepl(':', names(coef(x))))
     })

lm.2D.which=sapply(lm.2D, function(x) {
    sum(grepl(':', names(coef(x))))
     })





#mean(h2A[names(lm.2D)])
#mean(rowSums(count.matrix)[names(lm.2D)])
#[1] 580796

#significantly higher heritability
t.test(h2A, h2A[names(lm.2D)])
# not significantly higher total abundance
t.test(rowSums(count.matrix), rowSums(count.matrix)[names(lm.2D)])
geneList= names(lm.2D) 
for( thisOntology in c("BP", "MF", "CC") ) {

    testGenes= factor(0+(colnames(t.tpm.matrix) %in% geneList))
    names(testGenes)=colnames(t.tpm.matrix)
    GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
    GOresult = runTest(GOData, algorithm="classic", statistic="fisher")
    gt=GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult)))
    genesinterms=genesInTerm(GOData, gt[,1])
    genes.enriched.list=lapply(genesinterms, function(x) x[x%in%names(testGenes[testGenes==1])])
    genes.enriched.list.simple=lapply(genes.enriched.list, function(x) as.character(SYS2ORF.key[x]))
    gt$Genes=as.vector(sapply(genes.enriched.list.simple, paste, collapse=','))
    gt$GenesSystematic=   as.vector(sapply(genes.enriched.list, paste, collapse=','))
    
     pdf(file=paste0('/data/eQTL/RData/GO/', 'epistatic_transcripts', '_', thisOntology, '.pdf'), width=25, height=25)
         plotGOToTree(GOData, GOresult, sigThres = 0.00005)
      dev.off()

    write.table(gt, file=paste0('/data/eQTL/RData/GO/', 'epistatic_transcripts', '_', thisOntology, '.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')

}







#colMeans(vc6)
#colMeans(vc6)
#      A.local     A.distant             A   A.local_X_A A.distant_X_A            AA            In 
#     0.115295      0.176857      0.051674      0.007336      0.021050      0.020913      0.610759 
##sum(vc6[1:3])
#R> sum(vc6[1:3])
#[1] 0.233
#mean((vc6[,1]+vc6[,2])/(vc6[,1]+vc6[,2]+vc6[,3]))
#[1] 0.8178

barplot(colMeans((do.call('rbind', sapply(vc_multiple_components, function(x) x[['4']]$sigma)))))
# model without qAA 
qA_A_AA=sapply(vc_multiple_components, function(x) x[['3']]$sigma)
qA_A_qAA_AA=do.call('cbind', sapply(vc_multiple_components, function(x) x[['4']]$sigma))
 hist(qA_A_AA['qA',]/(qA_A_AA['qA',]+ qA_A_AA['A',]))
hist(qA_A_AA['qA',]/(qA_A_AA['qA',]+ qA_A_AA['A',]), breaks=20, xlim=c(0,1))
 hist(qA_A_qAA_AA['qAA',]/(qA_A_qAA_AA['qAA',]+ qA_A_qAA_AA['AA',]))
qA_A_qAA_AA[3,]/
R> mean(qA_A_AA['qA',]/(qA_A_AA['qA',]+ qA_A_AA['A',]))
[1] 0.7814
R> median(qA_A_AA['qA',]/(qA_A_AA['qA',]+ qA_A_AA['A',]))
[1] 0.8041

# ---------------------------------------------------------------------------------------------------------------------------------------------

apdgc=split(all.peaks.DS, all.peaks.DS$chr)
app=lapply(apdgc, function(x) split(x, x$gene))

for(g in names(peaks.per.gene)){
    print(g)
    peaks.per.gene[[g]]= peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    apeaks = match( peaks.per.gene[[g]]$pmarker, colnames(gdata))

    #newMM.peaks[[i]]
    ax= paste('gdata[,', apeaks,']', sep='')
    aq=paste(ax, collapse= ' + ')
    i=match(g, colnames(pheno.scaled))
    # where pheno.scaled= t.tpm.matrix with batch as a fixed effect
    #am=lm(paste('pheno.scaled[,' , i,']' , ' ~ ', (aq), '-1'))

    am=lm(paste('scale(t.tpm.matrix[,' , i,'])' , ' ~ gbatch.fact + ', (aq), '-1'))
    #plot(LODS[i,])
    #abline(v=apeaks)

    aov.a = anova(am)
    tssq=sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(am)  
    peaks.per.gene[[g]]$peak.ind=apeaks
    peaks.per.gene[[g]]$var.exp=a.effs[-1]
    peaks.per.gene[[g]]$lm.coeff=as.vector(coeffs)[-c(1:13)]
    #va.a = (tssq-aov.a[nrow(aov.a),2])/(tssq)
    #lms[[g]]=summary(am)$adj.r.squared
    #lms[[g]]=scale(residuals(am))
    #summary(am)$adj.r.squared
    #coefficients(am)
}  
da=density(unlist(sapply(peaks.per.gene, function(x) x$var.exp[!x$cis] )))
di=density(unlist(sapply(peaks.per.gene, function(x) x$var.exp[x$cis] ) ))
#all.peaks.multiple.regression$var.exp[all.peaks.multiple.regression$cis])
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(da, xlab='fraction of phenotypic variance', ylab='density', xlim=c(0,0.5),  type='n', main='')
polygon(da, col=rgb(0, 0, 1,0.3), border=NA)
polygon(di, col=rgb(1, 0, 0,0.3), border=NA)
legend('topright', legend=c('distant', 'local'), fill=c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)))


all.peaks.multiple.regression=do.call('rbind', peaks.per.gene)


vexp=lm(scale(t.tpm.matrix)~covariates.OD+pc.ars[,1:7]-1)

vexp.residuals=scale(residuals(vexp))

vAresid=calcA(vexp.residuals, A)
covvexp=lm(scale(t.tpm.matrix)~covariates.OD-1)

pcvexp=sapply(summary(vexp), function(x) x$adj.r.squared)
covexp=sapply(summary(covvexp), function(x) x$adj.r.squared)

plot(vcA.OD[,1]/c(vcA.OD[,1]+vcA.OD[,2]), vAresid[,1]/c(vAres     rr=  regress(yr~1, ~qA+A+qAA+AA, verbose=T, pos=c(T,T,T,T,T))
id[,1]+vAresid[,2]), xlim=c(0,1), ylim=c(0,1),
     xlab='covariate and growth corrected h^2', ylab='covariate and growth and pc20 sv corrected h^2'   )

test=lm(t.tpm.matrix[,'YCL009C']~covariates.OD+pc.ars[,1:20]+gdata[,peaks.per.gene[['YCL009C']]$pmarker])


testprcomp=prcomp(pheno.additive.removed.scaled)
pc.ars=scale(testprcomp$x)
pcLOD=fasterLOD(1012, pc.ars, gdata.scaled)

#8th
pcAA=list()
for(i in 1:20) {
    print(i)
    pcAA[[as.character(i)]]=extract.rr(regress(pc.ars[,i]~1, ~A+AA, verbose=T, pos=c(T,T,T,T)))
}

A2.fit=list()
X=data.frame(gdata[,unlist(htabler[,'peak'])])
for( i in 9:20) {
 fitme=lm(pc.ars[,i]~.-1, data=X)
 A2.fit[[as.character(i)]]=add1(fitme, ~.^2, test='F')
 print(i)
 hist(A2.fit[[as.character(i)]][,i], breaks=1000)
}


# divide the genome until 1500 bins 
bins=seq(1,1.2e7,7500)
hp.null=cut(all.peaks.OD$marker.gcoord, bins)
lambda.pois=length(hp.null)/length(bins)
qpois(.05/length(bins), lambda.pois, lower.tail=F)
plot(qvalue(ppois(c(rle(sort(as.numeric(hp.null)))$lengths), lambda.pois, lower.tail=F)))

rle(as.vector(hp.null))
sort(unlist(do.call('rbind', lapply(hotspots.OD, function(x) do.call('rbind', x)))[,2]))

marker.GR[hp.index]$gcoord
plot(hp.null, ylim=c(0,200))
hist(all.peaks.OD$marker.gcoord, breaks=600, ylim=c(0,200))
abline(h=44)
#abline(v=cut(marker.GR[hp.index]$gcoord, bins))
#abline(v=match(cut(marker.GR[hp.index]$gcoord,bins), cut(bins,bins)))
abline(v=marker.GR[hp.index]$gcoord, col='red')

#load('/data/eQTL/RData/hotspots_042616.RData')
#Iteration 2 
#hotspots.refined=findHotspots.iteration2(hotspots.OD, t.tpm.matrix, gene.annot.df, 
#                                         peaks.per.chr, genetic.map, gdata, gdata.s.by.chr, gbatch.fact.OD ,do.save=T) 
#save(hotspots.refined, file='/data/eQTL/RData/hotspots_refined_060116.RData')
#load('/data/eQTL/RData/hotspots_refined_060116.RData')
#save(hotspots.refined, file='/data/eQTL/RData/hotspots_refined_042616.RData')
#load('/data/eQTL/RData/hotspots_refined_042616.RData')

cc=unique.chrs[[6]]
par(mfrow=c(1,1))
hist(peaks.per.chr[[cc]]$pcind, breaks=100, main=cc)
po=as.character(sapply(hotspots.OD[[cc]], function(x) x$peak))
#pr=as.character(sapply(hotspots.refined[[cc]], function(x) x$peak))

abline(v=match(po, colnames(gdata.by.chr[[cc]])), lty=2, col='blue')
text(match(po, colnames(gdata.by.chr[[cc]])), 0,
           1:length(po), col='red')
#hist(peaks.per.chr[[cc]]$pcind, breaks=100, main=paste(cc, 'refined'))
#abline(v=match(pr, colnames(gdata.by.chr[[cc]])), lty=2, col='blue')
#text(match(pr, colnames(gdata.by.chr[[cc]])), 0,
 #          1:length(pr), col='red')


# on Juno
#hotspot.boots=list()
#for( i in 43:1) {    #85:1 # 64:1 43:1
#    print(i)
#    hotspot.boots[[i]]=bootstrap.Hotspots(i, onjuno=TRUE)  
#}

#sapply(sapply(split(marker.GR, seqnames(marker.GR)), function(x) (ranges(x[length(x),]))), function(x) x[[1]][1])
#   chrI   chrII  chrIII   chrIV    chrV   chrVI  chrVII chrVIII   chrIX    chrX   chrXI  chrXII chrXIII  chrXIV   chrXV  chrXVI 
# 202825  798782  302979 1521369  563829  266145 1068261  519219  424866  724256  643662 1058607  914575  765330 1067520  928107


phname=cbind(as.character(do.call('rbind', 
                               lapply(hotspots.refined, function(x) do.call('rbind',x)))[,3]))
pdf(file='/data/eQTL/RData/hotspots_ints.pdf', width=10, height=10)
    for(nn in c(1:85)[-44]) {    hist(hotspot.boot.peaks[[nn]], breaks=100, main=phname[nn], sub=hotspot.boot.intervals[nn,4]-hotspot.boot.intervals[nn,2])
      abline(v=hotspot.boot.intervals[nn,2], col='blue', lwd=2, lty=2)
      abline(v=hotspot.boot.intervals[nn,4], col='blue', lwd=2, lty=2)
    }
dev.off()

hist(peaks.per.chr[[12]]$pcind, breaks=1000)
do.call('rbind', hotspots.OD[[12]])$peak

#----------------------------------------------------------------------------------------------------------
big.plot1='/home/jbloom/Dropbox/Public/eQTL/eQTL_map_1024x1024.png'
png(file=big.plot1, width=1024, height=1024)
eQTL_bigPlot(all.peaks, gcoord.key, marker.GR)
dev.off()

big.plot2='/home/jbloom/Dropbox/Public/eQTL/eQTL_map_1920x1920.png'
png(file=big.plot2, width=1920, height=1920)
eQTL_bigPlot(all.peaks.OD, gcoord.key, marker.GR)
dev.off()

big.plot3='/home/jbloom/Dropbox/Public/eQTL/eQTL_map_10x10.pdf'
pdf(file=big.plot3, width=10, height=10)
eQTL_bigPlot(all.peaks, gcoord.key, marker.GR)
dev.off()

for(chr in unique.chrs[-17]) {
big.plot.by.chr=paste0('/home/jbloom/Dropbox/Public/eQTL/eQTL_map_', chr, '.pdf')
pdf(file=big.plot.by.chr, width=10, height=10)
eQTL_bigPlot(all.peaks, gcoord.key, marker.GR, xlim.ind= which(as.character(seqnames(marker.GR)) %in% chr)   )
dev.off()

}
eQTL_bigPlot(all.peaks, gcoord.key, marker.GR, xlim.ind= which(as.character(seqnames(marker.GR)) %in% chr)   )
xlim.ind=which(as.character(seqnames(marker.GR)) %in% 'chrVIII'
htable=do.call('rbind', lapply(hotspots, function(x) do.call('rbind', x)))

htabler=do.call('rbind', lapply(hotspots.refined, function(x) do.call('rbind', x)))
png(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/030816/eQTL_hotspots_crunch.png', width=1024, height=256)
eQTL_bigPlot(all.peaks, gcoord.key)
abline(v=marker.GR$gcoord[match(htabler$peak, marker.GR$mname)])
dev.off()

big.plot.by.chr=paste0('/home/jbloom/Dropbox/Public/eQTL/eQTL_map_chrVIII_fixed_axis.pdf')
pdf(file=big.plot.by.chr, width=10, height=10)

dev.off()

peaks.per.gene=split(all.peaks.OD, all.peaks.OD$gene)
marker.GR.hot=marker.GR[marker.GR$mname %in% hotspot.peaks,]
par(mfrow=c(4,4)) 
for(i in unique.chrs[1:16]) {
plot(start(marker.GR[seqnames(marker.GR)==i]), genetic.map[[i]], xlab='physical position', ylab='genetic map position', main=i)
abline(v=start(marker.GR.hot[seqnames(marker.GR.hot)==i]), col='red' )
} 



# HOTSPOT ANALYSIS --------------------------------------------------------------------------------------------------
#glm(counts[i,]~offset(log(tot.counts)), family='poisson')
#note 'LL' = log10(likelihood)
#Mediation_LOD = LOD32 - LOD10

 # Frank's contrasts fl
#LL3   = (distant.transcript ~ covariates + QTL  )
#LL2   = (distant.transcript ~ covariates        )
#LL1   = (distant.transcript ~ covariates + QTL + local.transcript )
#LL0   = (distant.transcript ~ covariates + local.transcript       ) 

#For Mediation
#LOD32 = LL3 - LL2
#LOD10 = LL1 - LL0

#LOD13  = LL1 - LL3
#LOD12  = LL1 - LL2

##case 1 is LOD12 >>> LOD13 (#if theres' a very big effect of the QTL)
#case2 is   LOD12   >  LOD13  (if the local effect is larger)

# Mediation = QTL effect smaller if you add in local transcript as covariate in model 
#1 vs 0
library(intermediate)
library(foreach)
library(doMC)
registerDoMC(70)

load('/data/eQTL/RData/R_genesInHotspotList_95PadWLD_160902.RData')
#loads genesInHotspotList

#local.mediation=list()
hotspot.peaks=do.call('c', htabler[,'peak']) #as.character(htabler$peak)
local.mediation=foreach(h=1:length(hotspot.peaks)) %dopar%  {
   #test with mat hotspot
 #  h=8
   hp=hotspot.peaks[h]
   print(hp)
   mint=marker.GR[marker.GR$mname==hp]
   # all genes within 10 kb of peak marker
   #make this informed by the intervals ###
   genes.near.hotspot.int=subsetByOverlaps(gene.GR, mint+10000)
   genes.near.hotspot95=gene.GR[na.omit(match(as.vector(genesInHotspotList[[h]][[1]]), gene.GR$ORF)),]
    
   #  if(length(genes.near.hotspot95)> length(genes.near.hotspot.int) ) { 
   #      genes.near.hotspot=genes.near.hotspot95 } else{
   #      genes.near.hotspot=genes.near.hotspot.int       
   #  }
   genes.near.hotspot=unique(c(genes.near.hotspot.int, genes.near.hotspot95)) 

   null= lm(t.tpm.matrix~covariates.OD)
   full= lm(t.tpm.matrix~covariates.OD+gdata[,hp])
   rn=residuals(null)
   fn=residuals(full) 
   rss1=colSums(rn^2)
   rss2=colSums(fn^2)
   Fstat=(rss1-rss2)/(rss2/(1012-14))
   hotspot.model=pf(Fstat, 1,1012-14, lower.tail=F)
   qhot=qvalue(hotspot.model)
   sig.at.hotspot=which(qhot$qvalue<.05)
   
   all.hotspot.genes=gene.GR[match(names(sig.at.hotspot), gene.GR$ORF)]
   
   print(length(genes.near.hotspot))
   print(length(sig.at.hotspot))
   
   # creat new variable called hotspot.genes 
   hotspot.genes=genes.near.hotspot
   hotspot.genes=unique(c(all.hotspot.genes, genes.near.hotspot))
   hotspot.genes=sort(hotspot.genes, by= ~gcoord)

   cis.mediation.LOD=array(NA, c(length(sig.at.hotspot), length(hotspot.genes), 10) )
   dimnames(cis.mediation.LOD)[[1]]=names(sig.at.hotspot)
   dimnames(cis.mediation.LOD)[[2]]=hotspot.genes$ORF
   dimnames(cis.mediation.LOD)[[3]]=c('L0', 'L1', 'L2', 'L3', 'L32', 'L10', 'L12', 'L13', 'L41', 'mLOD')
   for(cg in 1:length(hotspot.genes)) {
        #print(cg)
        cisg=hotspot.genes[cg]$ORF
        
        sya=t.tpm.matrix[,sig.at.hotspot]
        rmcis.ind1=match(cisg, colnames(sya))
        if(!is.na(rmcis.ind1)) {
            sy=sya[,-rmcis.ind1]
        } else {sy=sya} 

        rmcis.ind=match(cisg, colnames(t.tpm.matrix))
        cis.e=(t.tpm.matrix[,rmcis.ind])
        #gdist=gdata[,closest.marker.to.transcript[cisg]]

        #testing if growth (OD)  mediates QTL effects on expression levels  -----------------------
        #sy=scale(t.tpm.matrix[,gind])
        #sOD=scale(OD.cov)
        #rmt4=residuals(lm(sy[,111]  ~ covariates.OD))
        #M4.t =lm(rmt4  ~ cis.e*gdata[,hp])
        #plot3d( cis.e, gdata[,hp], rmt4)
        
        #grid.lines = 1000
        #x.pred <- seq(min(cis.e), max(cis.e), length.out = grid.lines)
        #z.pred <- seq(min(gdata[,hp]), max(gdata[,hp]), length.out = grid.lines)
        #xz <-data.frame(x=x.pred, y=z.pred)
        
        #expand.grid( x = x.pred, z = z.pred)
        #y.pred <- matrix(predict(M4.t, newdata = xz), 
        #         nrow = grid.lines, ncol = grid.lines)
        #points3d(cis.e, gdata[,hp],(predict(M4.t)), color = "steelblue"), 
        #        alpha = 0.5, lit = FALSE)  
        # add in gdist???
        #if( (cisg %in% genes.near.hotspot$ORF) | 
        # if( is.na(closest.marker.to.transcript[cisg]) ) {
        
        # if local linker then recenter hotspot to be marker closest to the gene
        hp.new=hp
        if( cisg %in% genes.near.hotspot$ORF ) {hp.new=colnames(gdata)[closest.marker.to.transcript[cisg]] } 
    
             print(paste(cisg, hp.new))
            M4 =lm(sy  ~ covariates.OD+cis.e*gdata[,hp.new])
            M3= lm(sy  ~ covariates.OD+gdata[,hp.new])
            M2= lm(sy  ~ covariates.OD)
           
            M1= lm(sy  ~ covariates.OD+cis.e+gdata[,hp.new])
            M0= lm(sy  ~ covariates.OD+cis.e)
            
            L4= -(nrow(sy)/2)*log10(colSums(residuals(M4)^2))
            L3= -(nrow(sy)/2)*log10(colSums(residuals(M3)^2))
            L2= -(nrow(sy)/2)*log10(colSums(residuals(M2)^2))
            L1= -(nrow(sy)/2)*log10(colSums(residuals(M1)^2))
            L0= -(nrow(sy)/2)*log10(colSums(residuals(M0)^2))
            
            L32= L3-L2
            L10 = L1-L0
            L13 = L1-L3
            L12 = L1-L2
            L41 = L4-L1

            mLOD=(L32)-(L10)
            m=cbind(L0, L1, L2, L3, L32, L10, L12, L13, L41, mLOD)
      # reverse gives the same darn lod score
      #  M3r= lm(cis.e  ~ covariates.OD+gdata[,hp])
      #  M2r= lm(cis.e  ~ covariates.OD)
       
      #  M1r= lm(cis.e  ~ covariates.OD+sy[,1]+gdata[,hp])
      #  M0r= lm(cis.e  ~ covariates.OD+sy[,1])

      #  L3r= -(nrow(sy)/2)*log10(sum(residuals(M3r)^2))
      #  L2r= -(nrow(sy)/2)*log10(sum(residuals(M2r)^2))
      #  L1r= -(nrow(sy)/2)*log10(sum(residuals(M1r)^2))
      #  L0r= -(nrow(sy)/2)*log10(sum(residuals(M0r)^2))
      #  L32r= L3r-L2r
      #  L10r = L1r-L0r
      #  L13r = L1r-L3r
      #  L12r = L1r-L2r
      # mrLOD=(L32r)-(L10r)

        #M4.list=list()
        ## from Schadt 2009
        #for(DiTr in colnames(sy)) {
        ## iterate through for each distant gene and calculate (2)
        #    M4 = lm(cis.e~covariates.OD+gdata[,hp]+sy[,DiTr])
        if(!is.na(rmcis.ind1)) {
            cis.mediation.LOD[-rmcis.ind1,cg,]=m
        } else {cis.mediation.LOD[,cg,]=m} 
    }
    attr(cis.mediation.LOD, 'local.linker')=colnames(cis.mediation.LOD) %in% genes.near.hotspot$ORF
    return(cis.mediation.LOD)
}
names(local.mediation)=hotspot.peaks
#save(local.mediation, file='/data/eQTL/RData/local.mediation.11302016.RData')
save(local.mediation, file='/data/eQTL/RData/local.mediation.all.linking3.RData')
load('/data/eQTL/RData/local.mediation.RData')

iloc=sapply(local.mediation, function(x) x[,,'L41'])
lodToPval <-  function(x){  pchisq(x*(2*log(10)),df=1,lower.tail=FALSE)/2}
sapply(iloc, function(x) which(x>10, arr.ind=T))[[37]]
lapply(local.mediation, function(x) colSums(apply(x, 2, function(y) {y[,'mLOD']>10  & (y[,'L12'] - y[,'L13'])<20 }), na.rm=T) )


library(fields) # cfor two.colors()
library(gplots)
pdf(file='~/Desktop/mLOD_vis5.pdf', width=10, height=10)
for(hpoi in names(local.mediation)) {
par(oma=c(4,2,2,2))
jj=local.mediation[[hpoi]]

    #jj=cis.mediation.LOD
    #hpoi=names(local.mediation)[8]
    lm13=jj[,,'L13']
    lmh =jj[,,'mLOD']
    lm14=jj[,,'L41']
    lm10=jj[,,'L10']

    #lmh[lmh<0]=0
    #lmh[lm13<10]=0
    #lmh[lm14<10]=0
    lmh[lmh < 5] <- 0
    lmh[lm13 < 5] <- 0
    lmh[lm10 > 3] <- 0
    #lm14[lm14<3]=0
    # x is mediated (linked)
    # y is mediating
    lmh[cbind(1:nrow(lmh), match(rownames(lmh), colnames(lmh)))]=NA

    lmh2=lmh[,colSums(lmh,na.rm=T)>50]
    llstat=attr(jj, 'local.linker')[colSums(lmh,na.rm=T)>50]
    
tryCatch({
    heatmap.2((t(lmh2))^(1/3),col= two.colors(n=256, start="white", end="red", middle="pink", alpha=1), 
          scale='none', trace='none', main=hpoi , key.xlab='(mLOD)^(1/3)', colRow=ifelse(llstat, 'red', 'black'))
}, error=function(e) {} )
#readline()
}
dev.off()
for(hpoi in names(local.mediation)) {
    hpoi=names(local.mediation)[1]
    lmed.dim=dim(local.mediation[[hpoi]])
    #t.space=t(combn(1:lmed.dim[2],2)) #
    t.space=expand.grid(1:lmed.dim[2], 1:lmed.dim[2])
    t.space=t.space[t.space[,1]!=t.space[,2],]
    p.ks=apply(t.space, 1, function(x) { ks.test(local.mediation[[hpoi]][,x[1],'mLOD'], local.mediation[[hpoi]][,x[2],'mLOD'], alternative='greater')$p.value} ) 
    p.ks.df=cbind(t.space, p.ks)
    ecdfs=list()
    for (n in 1:lmed.dim[2]){
        ecdfs[[as.character(n)]]=ecdf(local.mediation[[hpoi]][,n,'mLOD'])
    } 
    limits=quantile(local.mediation[[hpoi]][,,'mLOD'], c(.001, .999), na.rm=T)
    plot(ecdfs[[1]],verticals=T, do.points=FALSE, xlim=limits)
    for (n in 1:lmed.dim[2]){
          plot(ecdfs[[n]],verticals=T, do.points=FALSE, add=T)
    }
}    
    




plot(ecdf(local.mediation[[1]][,1,'mLOD']))
points(ecdf(local.mediation[[1]][,1,'mLOD']))


lmLod=lapply(local.mediation, function(x) apply(x, 2, function(y) {y[,'mLOD'] }) )

l12_l13=lapply(local.mediation, function(x) apply(x, 2, function(y) {y[,'L12'] - y[,'L13'] }) )
l12=lapply(local.mediation, function(x) apply(x, 2, function(y) {y[,'L12']  }) )
l13=lapply(local.mediation, function(x) apply(x, 2, function(y) {y[,'L13'] }) )

plot(unlist(l12), unlist(l13), col=(unlist(lmLod)>5)+1, xlim=c(0,1000), ylim=c(0,1000) )

which(local.mediation[['chrVIII:116378_C/A']][,'YHR005C','mLOD']>5)
t.select=names(which(local.mediation[['chrVIII:116378_C/A']][,'YHR005C','mLOD']>5))
ytemp=t.tpm.matrix[,t.select[9]]
gtemp=gdata[,'chrVIII:116378_C/A']
mtemp=t.tpm.matrix[,'YHR005C']

plot(jitter(gtemp), ytemp)
plot((mtemp), ytemp)

      , 2, function(y) {which(y[,'mLOD']>5)})
local.mediation[['chrVIII:116378_C/A']][,'YHR005C',]
YHR005C
loc.med
g='chrVIII:116378_C/A'


library(foreach)
library(doMC)
registerDoMC(60)

### REBOOT

#hp.med=list()
#load('/data/eQTL/RData/121216.RData')

sy=t.tpm.matrix
M2= lm(sy  ~ covariates.OD)
L2= -(nrow(sy)/2)*log10(colSums(residuals(M2)^2))
# lmh[lmh < 5] <- 0
#    lmh[lm13 < 5] <- 0
#    lmh[lm10 > 3] <- 0

hp.med=foreach (cisg =colnames(t.tpm.matrix)[1:5718] ) %dopar% {
    #cisg='YGL169W'
    #cisg='YCR040W'
    print(cisg)
    cis.e=t.tpm.matrix[,cisg]
    hp.new=closest.marker.to.transcript[cisg]
    M4 =lm(sy  ~ covariates.OD+cis.e*gdata[,hp.new])
    M3= lm(sy  ~ covariates.OD+gdata[,hp.new])
    M1= lm(sy  ~ covariates.OD+cis.e+gdata[,hp.new])
   #M2= lm(sy  ~ covariates.OD)
    M0= lm(sy  ~ covariates.OD+cis.e)

    L4= -(nrow(sy)/2)*log10(colSums(residuals(M4)^2))
    L3= -(nrow(sy)/2)*log10(colSums(residuals(M3)^2))
    L1= -(nrow(sy)/2)*log10(colSums(residuals(M1)^2))
    L0= -(nrow(sy)/2)*log10(colSums(residuals(M0)^2))
    L32 = L3-L2
    L10 = L1-L0
    L13 = L1-L3
    L12 = L1-L2
    L41 = L4-L1
    mLOD=(L32)-(L10)
    return( cbind(L10, L13, L32, L41, mLOD))
    gc()
}
names(hp.med)=colnames(t.tpm.matrix)[1:5718] 


l10=unlist(lapply(hp.med, function(x) x[,'L10'] ))
l13=unlist(lapply(hp.med, function(x) x[,'L13'] ))
ml=unlist(lapply(hp.med, function(x) x[,'mLOD'] ))
l41=unlist(lapply(hp.med, function(x) x[,'L41'] ))
med.tot =sapply(hp.med, function(x)   sum(x[,'mLOD']>4 & x[,'L13']>4 & x[,'L10']<4 & x[,'L13']<10000) )
hp.reg  =sapply(hp.med, function(x)   sum(x[,'L32']>4)  )
mod.tot =sapply(hp.med, function(x)   sum(x[,'L41']>4 ) ) 

#modu=lapply(hp.med, function(x) x[,'L41'])

#& x[,'L13']>4 & x[,'L10']<3) )

plot(gene.GR$gcoord[1:5718], med.tot, type='h', ylim=c(-2000,900))
points(gene.GR$gcoord[1:5718], -(hp.reg), type='h', lty=2, lwd=.5)
segments(marker.GR$gcoord[marker.GR$mname %in%  hotspot.peaks], -50, 
         marker.GR$gcoord[marker.GR$mname %in%  hotspot.peaks], 50, col='red')
abline(v=gcoord.key, col='blue')

gmed=gene.GR[1:5718]
#cisModel[i]=anova(nm,fm)$'Pr(>F)'[2]

#cisModel.effect[i]=as.numeric(coef(fm)[16])

gmed$cis.p=cisModel
gmed$cis.effect=cisModel.effect
gmed$hp.tot=hp.reg
gmed$med.tot=med.tot
gmed$mod.tot=mod.tot

write.table(gmed, file='~/Desktop/gmed.txt', sep='\t', quote=F, row.names=F)
chr='chrI'

pdf(file='~/Desktop/hotspots_vs_mediation.pdf', width=11, height=6)
for(chr in unique.chrs[-17]) {
par(yaxs='i')
logicme=as.vector(seqnames(gmed)==chr)
#par(mfrow=c(1,2))
#yr=c(-max(gmed$hp.tot[logicme]), max(gmed$hp.tot[logicme]))
yr=c(0, max(gmed$hp.tot[logicme]))
plot( start(gmed)[logicme], gmed$hp.tot[logicme],
    # c(start(gmed)[logicme], start(gmed)[logicme]),
    # c(-gmed$hp.tot[logicme], gmed$hp.tot[logicme]), 
     xlab='position', ylab='linkages', main=paste(chr, 'mediation'),
      type='h', ylim=yr)
#abline(h=0)
hotpos=marker.GR[marker.GR$mname %in%  hotspot.peaks & seqnames(marker.GR)==chr]
#points(start(gmed)[logicme],  -gmed$mod.tot[logicme], col='blue',      type='h', lwd=3)
points(start(gmed)[logicme][gmed$med.tot[logicme]>0],   gmed$med.tot[logicme][gmed$med.tot[logicme]>0], col='red',      type='h', lwd=3)
text(start(hotpos),0, '*', col='blue', cex=5)
#points(start(gmed)[logicme], log10(gmed$cis.p)[logicme], col='purple', type='h', lwd=3)
#points(start(gmed)[logicme][gmed$mod.tot[logicme]>10],  -gmed$mod.tot[logicme][gmed$mod.tot[logicme]>10], col='blue',  type='h', lwd=2)
#readline()
}
dev.off()

points(start(gmed)[logicme], -gmed$mod.tot[logicme], col='blue', type='h')

abline(h=0)



g1=doGO(colnames(t.tpm.matrix), names(which(mod.tot>100)))


g1=doGO(colnames(t.tpm.matrix), names(which(med.tot>100)))




library(topGO)
geneList= names(which(hp.med[['YNL117W']][,'L41']>4)) #names(lm.2D) 
geneList= names(which(hp.med[['YNL177C']][,'L41']>4)) #names(lm.2D) 
geneList= names(which(hp.med[['YNL250W']][,'L41']>4)) #names(lm.2D) 
geneList= names(which(hp.med[['YNL301C']][,'L41']>4)) #names(lm.2D) 
geneList= names(which(hp.med[['YNL301C']][,'L41']>4)) #names(lm.2D) 

goi='YMR029C'
goi='YLR360W'

goi='YLL008W'
goi='YLL007C'

goi='YGL169W'

goi='YDR265W'
goi='YDL132W'
goi='YDR117C'
goi='YHR178W'


goi='YKL107W'
goi='YKL109W'
goi='YJR060W'
goi='YGL193C'



# cis gene explains some variation and (cis gene +  genotype explains little compared to just cis gene)
goi='YAL054C'

geneList=names(which((hp.med[[goi]][,'mLOD']>4 & hp.med[[goi]][,'L13']>4 & hp.med[[goi]][,'L10']<4)))
geneList= names(which(hp.med[[goi]][,'L41']>4)) #names(lm.2D) 


g1=doGO(colnames(t.tpm.matrix), geneList)

plot(-log10(gmed$cis.p), gmed$med.tot)
#cor.test(gmed$cis.p, gmed$med.tot, use='pairwise.complete.obs', method='spearman')
#rho=-0.87




local.no.hotspot=read.delim('/data/eQTL/RData/localPeaksNoHotspotAnnotated_161117.txt', header=T, sep='\t', stringsAsFactors=F)
lnh=local.no.hotspot[is.na(local.no.hotspot$paralog),]
dfg=data.frame(DataFrame(gmed))

l1=dfg$med.tot[dfg$X.ORF %in% lnh$gene]
l2=dfg$med.tot[!(dfg$X.ORF %in% lnh$gene)]


mlnh=merge(lnh, dfg, by.x='gene', by.y='X.ORF', sort=F, all.x=T)

head(mlnh,50)$med.tot

smnull=replicate(1000, {sample(dfg$med.tot,50) })







, med.tot, type='h', ylim=c(0,200))

mL=(L3-L2)-(L1-L0)
l13=L1-L3
l10=L1-L0
sum(mL>5 & l13>5 & l10>3 ) 


 lmh[lmh < 5] <- 0
    lmh[lm13 < 5] <- 0
    lmh[lm10 > 3] <- 0



SS4=crossprod(scale(residuals(M4)))
SS3=crossprod(scale(residuals(M3)))
SS2=crossprod(scale(residuals(M2)))
SS1=crossprod(scale(residuals(M1)))
SS0=crossprod(scale(residuals(M0)))
 
#t4=abs(eigen(SS4)$values)

 t3=abs(eigen(SS3)$values)
 pc.to.retain3=which(cumsum(t3/sum(t3))<.999)
 t2=abs(eigen(SS2)$values)
 pc.to.retain2=which(cumsum(t2/sum(t2))<.999)
 t1=abs(eigen(SS1)$values)
 pc.to.retain1=which(cumsum(t1/sum(t1))<.999)
 t0=abs(eigen(SS0)$values)
 pc.to.retain0=which(cumsum(t0/sum(t0))<.999)

 l3=length(pc.to.retain3)
 l2=length(pc.to.retain2)
 l1=length(pc.to.retain1)
 l0=length(pc.to.retain0)
 wml=which.min(c(l3,l2,l1,l0))
 if(wml==1) {minpc=pc.to.retain3}
 if(wml==2) {minpc=pc.to.retain2}
 if(wml==3) {minpc=pc.to.retain1}
 if(wml==4) {minpc=pc.to.retain0}


sum(log(t3[minpc]))- sum(log(t2[minpc]))
sum(log(t1[minpc]))- sum(log(t0[minpc]))

sum(log(t4[minpc]))-sum(log(t1[minpc]))


# small model over full model
# error variation / (error variation + explained variation)

(sum(log(t3[minpc]))-sum(log(t2[minpc])))-(sum(log(t1[minpc]))-sum(log(t0[minpc])))

(sum(log(t2[minpc]))/sum(log(t3[minpc])))
(sum(log(t0[minpc]))/sum(log(t1[minpc])))





         


































# sanity checking that this acting as a LOD score
        #n.pheno=1012
        #ptest= scale(sy[,1])
        #gtest = scale(gdata[,hp])
        #r=crossprod(ptest, gtest)/(n.pheno-1)
        #LOD=(-n.pheno*log(1-r^2))/(2*log(10))
        #58.13
        #MF= -(1012/2)*log10(colSums(residuals(lm(sy  ~ gdata[,hp]))^2))
        #MN= -(1012/2)*log10(colSums(residuals(lm(sy  ~ 1))^2))
        # all good 

        # well this sucks ------------------------------------------
        #rmLOD = rep(NA, length(colnames(sy)))
        #for(k  in 1:length(colnames(sy))) {
        #    rM3= lm(cis.e  ~ covariates.OD+gdata[,hp])
        #    rM2= lm(cis.e  ~ covariates.OD)
        #    rM1= lm(cis.e  ~ covariates.OD+sy[,k]+gdata[,hp])
        #    rM0= lm(cis.e  ~ covariates.OD+sy[,k])

         #   rL3= -(nrow(sy)/2)*log10(sum(residuals(rM3)^2))
         #   rL2= -(nrow(sy)/2)*log10(sum(residuals(rM2)^2))
         #   rL1= -(nrow(sy)/2)*log10(sum(residuals(rM1)^2))
         #   rL0= -(nrow(sy)/2)*log10(sum(residuals(rM0)^2))
         #   rmLOD[k]=(rL3-rL2)-(rL1-rL0)
        #}
        # --------------------------------------------------------
        # tau = hotspot effect in model without mediator
        # tau.prime = hotspot effect in model with mediator











hchrom=4
g1=gdata.by.chr[[hchrom]]
ld2=cor(g1)^2
a1=all.peaks[all.peaks$chr==unique.chrs[hchrom],]
ra1=rle(sort(a1$marker.gcoord))
crp=colorRampPalette(rev(c('white', 'blue', 'green', 'yellow', 'orange', 'red')))
c.coords=which(as.character(seqnames(marker.GR)) %in% unique.chrs[hchrom])
c.gcoord=marker.GR$gcoord[c.coords]
cgc=cbind(c.gcoord, 0)
cgc[match(ra1$values, cgc[,1]),2]=ra1$lengths
l=LDheatmap(ld2, flip=TRUE, color=crp(100), add.map=T,genetic.distances=cgc[,1]) #matlab.like2(30), add.map=F)
cvec=as.vector(cgc[,2])
l2=LDheatmap.addScatterplot(l, cvec , height=.75)
hr.peaks=as.character(sapply(hotspots.refined[[hchrom]], function(x)x$peak))
#hs.gc=match(hr.peaks, marker.GR$mname)# marker.GR$gcoord[match(hr.peaks, marker.GR$mname)]
#hs.gc= marker.GR$gcoord[match(hr.peaks, marker.GR$mname)]
l3=LDheatmap.marks(l2, hs.gc, hs.gc, col='black')

x11()
plot(cgc[,1], cgc[,2])
abline(v=marker.GR$gcoord[match(hr.peaks, marker.GR$mname)])

eQTL_bigPlot(all.peaks, gcoord.key, marker.GR, which(as.character(seqnames(marker.GR)) %in% unique.chrs[1]))
abline(v=marker.GR$gcoord[match(hr.peaks, marker.GR$mname)])




































#save(hotspots, file='/data/eQTL/RData/hotspots_112516.RData')
for(chr.inc in unique.chrs[-17]) {

png(file=paste0('/data/eQTL/plots/hotspots/hotspots_', chr.inc, '.png'), width=1080, height=1080)
#chr.inc=unique.chrs[-17]
#chr.inc='chrXVI'
#par(xaxs='i', yaxs='i')
plot(genetic.map.c[mgi[all.peaks$chr %in% chr.inc]],pgi[all.peaks$chr %in% chr.inc], pch=19, cex=all.peaks$LOD[all.peaks$chr %in% chr.inc], type='n', xlab='genetic distance (cm)', ylab='transcript index', main=chr.inc) 
#, xlim=c(3700,4500))#, col='grey') ##, xlim=c(3750,5500))
#segments(mgiL,pgi, mgiR,pgi, lwd=.5, col='#00000055' )
#, xlab='marker index', ylab='transcript index')#, col='grey') ##, xlim=c(3750,5500))

purp='#80008044'
orng='#FFA50044'
rd='#FF000044'

segments(genetic.map.c[mgiL[all.peaks$r>0 & all.peaks$chr %in% chr.inc]],pgi[all.peaks$r>0 & all.peaks$chr %in% chr.inc],
         genetic.map.c[mgiR[all.peaks$r>0  & all.peaks$chr %in% chr.inc]],pgi[all.peaks$r>0  & all.peaks$chr %in% chr.inc], lwd=.5, col='#80008055')
segments(genetic.map.c[mgiL[all.peaks$r<0 & all.peaks$chr %in% chr.inc ]],pgi[all.peaks$r<0  & all.peaks$chr %in% chr.inc ],
         genetic.map.c[mgiR[all.peaks$r<0 & all.peaks$chr %in% chr.inc]],pgi[all.peaks$r<0  & all.peaks$chr %in% chr.inc], lwd=.5, col='#FFA50055')
points(genetic.map.c[mgi[all.peaks$r>0 & all.peaks$chr %in% chr.inc ]],pgi[all.peaks$r>0 & all.peaks$chr %in% chr.inc ], pch=20, cex=log(all.peaks$LOD[all.peaks$chr %in% chr.inc], base=4.5), col=purp)
points(genetic.map.c[mgi[all.peaks$r<0 & all.peaks$chr %in% chr.inc ]],pgi[all.peaks$r<0 & all.peaks$chr %in% chr.inc ], pch=20, cex=log(all.peaks$LOD[all.peaks$chr %in% chr.inc],base=4.5), col=orng)
points(genetic.map.c[mgi[all.peaks$cis & all.peaks$chr %in% chr.inc ]],pgi[all.peaks$cis & all.peaks$chr %in% chr.inc ], pch=20, cex=log(all.peaks$LOD[all.peaks$cis & all.peaks$chr %in% chr.inc],base=4.5), col=rd)
#abline(v=hpeak.vec, lty=3, col='blue')
#abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='grey')
#abline(h=cumsum(rle(gene.annot.df$chr)$lengths), col='grey')
abline(v=genetic.map.c[match(sapply(hotspots[[chr.inc]], function(x)x$peak), colnames(gdata))])
#abline(v=genetic.map.c[match(sapply(hotspots[[chr.inc]], function(x)x$cI.l), colnames(gdata))])
gmcL=genetic.map.c[match(sapply(hotspots[[chr.inc]], function(x)x$cI.l), colnames(gdata))]
gmcR=genetic.map.c[match(sapply(hotspots[[chr.inc]], function(x)x$cI.r), colnames(gdata))]
text(gmcL-1, rep(5, length(gmcL)), 1:length(gmcL))
abline(v=gmcL, col='darkgrey', lty=2)
abline(v=gmcR, col='darkgrey', lty=2)
dev.off()
}

for(cc in unique.chrs[-17]) {

png(file=paste0('/data/eQTL/plots/hotspots/hotspots_marginal_', cc, '.png'), width=1080, height=512)

g=gdata.s.by.chr[[cc]]
XXr=crossprod(g)/(ncol(tpm.matrix)-1)
#XYr=crossprod(g,tmm)/(ncol(tpm.matrix)-1)
#scanone.all=(-1012*log(1-XYr^2))/(2*log(10))
fCX=sort(findCorrelation(XXr, cutoff=.99, exact=F))
g.s=g[,-fCX]

find=match(colnames(g.s), names(genetic.map[[cc]]))

gind=match(names(bl[[cc]]), gene.annot.df[,1])
q.bl.c=lapply(bl[[cc]], function(x) nrow(x$fit))
q.bl.c=do.call('c', q.bl.c)
greps=rep(gind,as.vector(q.bl.c))
bl.loc=sapply(bl[[cc]], function(x) x$fit[,1] )
bl.loc=sapply(bl[[cc]], function(x) x$fit[,1] )
bl.beta=sapply(bl[[cc]], function(x) x$fit[,'beta'] )
bl.var=sapply(bl[[cc]], function(x) x$fit[,'posterior variance'] )
bl.t=sapply(bl[[cc]], function(x) x$fit[,'t-value'] )
bl.p=sapply(bl[[cc]], function(x) x$fit[,'p-value'] )
#
plot(genetic.map[[cc]][find[unlist(as.vector(bl.loc))]], unlist(as.vector(bl.beta)), xlab='genetic distance (cm)', ylab='beta', main=cc)
abline(v=genetic.map.c[match(sapply(hotspots[[cc]], function(x)x$peak), colnames(gdata))], col='blue')
#abline(v=genetic.map.c[match(sapply(hotspots[[chr.inc]], function(x)x$cI.l), colnames(gdata))])
gmcL=genetic.map.c[match(sapply(hotspots[[cc]], function(x)x$cI.l), colnames(gdata))]
gmcR=genetic.map.c[match(sapply(hotspots[[cc]], function(x)x$cI.r), colnames(gdata))]
text(gmcL-1, rep(0, length(gmcL)), 1:length(gmcL))
abline(v=gmcL, col='darkgrey', lty=2)
abline(v=gmcR, col='darkgrey', lty=2)
dev.off()
}

hotspot.vector=sort(match(as.character(unlist( sapply(hotspots, function(x) sapply(x, function(y) as.character(y$peak))))), colnames(gdata)))


hotspot.lms=list()
for(i in 1:ncol(t.tpm.matrix)) {
    print(i)
    transcript=colnames(t.tpm.matrix)[i]
    cQTL=gdata[, closest.marker.to.transcript[i]]
     #yin=t.tpm.matrix[,mid.i]
    if(sum(is.na(cQTL))>0) {    cQTL[is.na(cQTL)]=1  }
    hotspot.lms[[transcript]]=lm(scale(t.tpm.matrix[,i])~gbatch.fact+cQTL+gdata[,hotspot.vector]-1)
}
ols.coefficient.matrix=sapply(hotspot.lms, coefficients)

#batch and cis de-regressed phenotypes
bc.resid=list()
for(i in 1:ncol(t.tpm.matrix)) {
    print(i)
    transcript=colnames(t.tpm.matrix)[i]
    cQTL=gdata[, closest.marker.to.transcript[i]]
     #yin=t.tpm.matrix[,mid.i]
    if(sum(is.na(cQTL))>0) {    cQTL[is.na(cQTL)]=1  }
    bc.resid[[transcript]]=scale(residuals(lm(scale(t.tpm.matrix[,i])~gbatch.fact+cQTL-1)))
}

bc.resid=do.call('cbind', bc.resid)
colnames(bc.resid)=rownames(tpm.matrix)

gdata.s=scale(gdata)

lasso.hotspots=list()
for(i in 391:ncol(bc.resid)) {
    print(i)
    lasso.hotspots[[colnames(bc.resid)[i]]]=cv.EBglmnet(gdata.s[,hotspot.vector], bc.resid[,i])
}

lasso.matrix=matrix(0, length(hotspot.vector), ncol(bc.resid) ) 
rownames(lasso.matrix)=colnames(gdata)[hotspot.vector]
colnames(lasso.matrix)=colnames(bc.resid)

for(i in 1:ncol(bc.resid)) {
    print(i)
    lh=lasso.hotspots[[colnames(bc.resid)[i]]]
    lasso.matrix[lh$fit[,1],i]=lh$fit[,3]
}

load('/data/eQTL/RData/hotspots_031016.RData')
hotspot.vector=sort(match(as.character(unlist( sapply(hotspots, function(x) sapply(x, function(y) as.character(y$peak))))), colnames(gdata)))
g.h=gdata.scaled[,hotspot.vector]

H.GRM=tcrossprod(g.h)/ncol(g.h)
vcH=calcA(pheno.scaled, H.GRM)
h2H=(vcH[,1]/(vcH[,1]+vcH[,2]))

plot(h2A, h2H)
gene.GR=GRanges(seqnames=(gene.annot.df$chr), ranges=IRanges(start=gene.annot.df$start, end=gene.annot.df$end), 
                strand=ifelse(gene.annot.df$strand==1, '+', '-'),
                ORF=(gene.annot.df$name)
                )

ma.chr=sapply(strsplit(colnames(gdata), ':'), function(x)x[1])
ma.pos=as.numeric( sapply( strsplit(sapply(strsplit(colnames(gdata), ':'), function(x)x[2]), '_'), function(x)x[1]))

marker.GR=GRanges(seqnames=ma.chr, ranges=IRanges(start=ma.pos, end=ma.pos) )
closest.marker.to.transcript=nearest(gene.GR, marker.GR)

p.batch.cis=pheno.scaled[,1:6288]
for(i in 1:6288) {
    print(i)
    p.batch.cis[,i]=scale(residuals(lm(t.tpm.matrix[,i]~gbatch.fact+gdata[,closest.marker.to.transcript[i]])))
}

#H.GRM=tcrossprod(g.h)/ncol(g.h)
vccH=calcA( p.batch.cis, H.GRM)
h2cH=(vccH[,1]/(vccH[,1]+vccH[,2]))

vccA=calcA( p.batch.cis, A)
h2cA=(vccA[,1]/(vccA[,1]+vccA[,2]))

plot(h2cA, h2cH)
identify(h2cA, h2cH)
colnames(p.batch.cis)[c(964, 2519, 2561, 2571, 3129,3756)]

x11()
plot(h2H[1:6288], h2cH)


x11()
plot(h2A[1:6288], h2cH) 









lcor=cor(lasso.matrix)

dlc=dist(t(ols.coefficient.matrix[-c(1:14),]))

dlc=dist(t(lasso.matrix))
dlc=as.dist(1-lcor^2)
dclust=hclust(dlc,method="ward.D")
cg=cutree(fit, k=ceiling(ncol(stpreds)/20))

library('RColorBrewer')
colvec=c('blue',  'white',  'red')
#, 'yellow', 'orange', 'red')
cramp=colorRampPalette(colvec)
cr=cramp(40)

image.plot(lasso.matrix[,dclust$order], col=cr)
image.plot(ols.coefficient.matrix[-c(1:14),][,dclust$order], col=cr)
image.plot(lasso.matrix[big.hotspots<40,dclust$order], col=cr) 
, col=cr)


#giant table of peaks from forward scan
save(all.peaks, file='/data/eQTL/plots/hotspots/peaks.RData')
#tpm per segregant (un-corrected)
save(t.tpm.matrix, file='/data/eQTL/plots/hotspots/tpm.RData')
#genotype data
save(gdata, file='/data/eQTL/plots/hotspots/geno.RData')
#genetic map
save(genetic.map, file='/data/eQTL/plots/hotspots/genetic_map.RData')
#factor for batch effect
save(gbatch.fact, file='/data/eQTL/plots/hotspots/gbatch_fact.RData')
#OD covariate (currently unused)
save(OD.cov, file='/data/eQTL/plots/hotspots/ODcov.RData')
#list of hotspot positions
save(hotspots, file='/data/eQTL/plots/hotspots/hotspots.RData')

#coefficients for OLS model with hotspot effects
save(ols.coefficient.matrix, file='/data/eQTL/plots/hotspots/hotspot_coefficients_new.RData')
save(lasso.matrix, file='/data/eQTL/plots/hotspots/hotspot_coefficients_lasso.RData')



#plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.var)))
#plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.t)))
#plot(unlist(as.vector(bl.loc)), -log10(unlist(as.vector(bl.p))))



#dev.off()







#match('10202890_chrXV_170921_A_G', colnames(gdata.s.by.chr[['chrXV']]))
#12670
#250

cc='chrXIV'
gforsim=gdata.s.by.chr[[cc]]

match('9623949_chrXIV_376313_C_T', colnames(gforsim))
match('9516279_chrXIV_268643_G_A', colnames(gforsim))
match('9714664_chrXIV_467028_A_G', colnames(gforsim))

n.replicates=500
h2=.015
nadditive=3
#add.qtl.ind=250
add.qtl.ind=c(284,385,482)
#482
add.qtl.sign=c(1,1,1)
a.eff=rep(0,1012)
a.eff=a.eff+add.qtl.sign[1]*gforsim[,add.qtl.ind[1]]+add.qtl.sign[2]*gforsim[,add.qtl.ind[2]]+add.qtl.sign[3]*gforsim[,add.qtl.ind[3]]

a.eff=scale(a.eff)
g=sqrt(h2)*a.eff
s=replicate(n.replicates, {g+rnorm(1012,mean=0,sd=sqrt((1-h2)/h2*var(g)))})[,1,]
s=scale(s)
simL=fasterLOD(1012, s, gforsim)
wm=apply(simL, 1, which.max)
m=apply(simL, 1, max)

par(xaxs='i')
plot(wm[m>2.5], seq_along(1:(sum(m>2.5))), xlim=c(0,ncol(gforsim)), col='#00000055',cex=log(m[m>2.5],base=4.5))
plot(wm, seq_along(1:length(m)), xlim=c(0,ncol(gforsim)), ylab='simulated traits', xlab='marker position', main =paste(h2/3, 'of phenotypic variance explained per QTL'))

#     pch='.')
#cex=log(m[m>2.5],base=4.5))
abline(v=add.qtl.ind)
#cex=.2)
    # pch='.')

#cex=log(m[m>2.5],base=4.5))
#cex=.2) #,
#IRA2
12672
12685


p=s
g=gforsim
n.marker.space=20
n.perms=5


sub.t=s
sub.t=residuals(lm(s~gdata.s.by.chr[[cc]][,386]))
sub.t=residuals(lm(s~gdata.s.by.chr[[cc]][,386]+gdata.s.by.chr[[cc]][,284]  ))
sub.t=residuals(lm(s~gdata.s.by.chr[[cc]][,386]+gdata.s.by.chr[[cc]][,284]+gdata.s.by.chr[[cc]][,482] ))

sub.t=residuals(lm(s~gdata.s.by.chr[[cc]][,388]+gdata.s.by.chr[[cc]][,284]+gdata.s.by.chr[[cc]][,484] ))

nn=ncol(sub.t)
   switchdet=FALSE
   sss=rep(NA, ncol(gdata.s.by.chr[[cc]]))
   pb =txtProgressBar(min = 1, max =length(sss), style = 3)
   
   #optimize ... SST doesn't need to be recalculated during this scan!!!!! at least 2X speedup 
   SST=crossprod((sub.t))
   dSST=determinant(SST,logarithm=T)$modulus
   switchdet=TRUE
   #if(!is.finite(dSST)){
          t5=svd(SST)
          switchdet=TRUE
          #pc.to.retain=which(t5$d/sum(t5$d)>(median(t5$d/sum(t5$d))))
           #pc.to.retain=1:100
          pc.to.retain=which(cumsum(t5$d/sum(t5$d))<.999)
          print(length(pc.to.retain))
           #  http://arxiv.org/abs/1305.5870
   #}
   for(i in 1:length(sss)){
       setTxtProgressBar(pb, i)
       test=lm(sub.t~gdata.s.by.chr[[cc]][,i])
       # Wilks.test(sub.t, gdata.s.by.chr[[cc]][,i])
       t2=residuals(test)
       #t3=predict(test)
       SSE=crossprod(t2)
       #t(t2)%*%t2
       #SST=crossprod(t2+t3)
       #if(!is.finite(dSSE) | !is.finite(dSST) )  {switchdet=TRUE}
       if(switchdet){
             t4=svd(SSE)
             #t5=svd(SST)
#       n=1012
#       p=nn = 500
       #g =2
#-(n-1-(p+g)/2)*ln(lambda)  = Xhi-square (p, g-1)
        #       http://www.statistik.tuwien.ac.at/forschung/CS/CS-2007-7.pdf
             X2=-(1012-1-(nn+2)/2)*(sum(log(t4$d[pc.to.retain]))-sum(log(t5$d[pc.to.retain])))

       }else{  
           dSSE=determinant(SSE,logarithm=T)$modulus
           X2=-(1012-1-((nn+2)/2))*(dSSE-dSST)      }
        #log((prod(as.brob(t4$values))/(prod(as.brob(t5$values)))))
        sss[i]=X2/4.605
    #    sss[i]=(-1012*log((prod(as.brob(t4$values))/ (prod(as.brob(t5$values))))))/4.605
    #   print(sss[i])
    }
   close(pb)
   which.max(sss)
plot(-log10(pchisq(sss*4.605,nn,lower.tail=F)))

s1=sss
s2=sss
s3=sss
s4=sss

plot(-log10(pchisq(s1*4.605, nn, lower.tail=F)), type='l', xlab='marker position', ylab='-log10(p)')
points(-log10(pchisq(s2*4.605, nn, lower.tail=F)), col='blue', type='l')
points(-log10(pchisq(s3*4.605, nn, lower.tail=F)), col='red', type='l')


ldg=cor(gforsim)
LL=matrix(NA,2000,774)
mLOD=rep(NA,774)
#SSEf=matrix(NA,2000,774)
#SSEr=matrix(NA,2000,774)
#LLn=list()
for(i in 1:774) {
    print(i)
    r.ld=range(which(ldg[i,]>.5))
    mdiff=abs(r.ld-i)
    rends=ldg[i,][r.ld]
    XX=cbind(gforsim[,r.ld[1]], gforsim[,i],gforsim[,r.ld[2]])

   # print(r.ld)
    if(rends[1]>.7) {
          #print('left')
          r.f=residuals(lm.fit(XX[,c(2,3)], s))
          r.r=residuals(lm.fit(as.matrix(XX[,3]), s))
     }else {
         if(rends[2]>.7 ){
          #print('right')
          r.f=residuals(lm.fit(XX[,c(1,2)],s))
          r.r=residuals(lm.fit(as.matrix(XX[,1]), s))
        }
        else {
         # print('hi')
         r.f=residuals(lm.fit(XX, s))
         r.r=residuals(lm.fit(XX[,c(1,3)], s))
        }
    }
    SSEf=apply((r.f)^2, 2, sum)
    SSEr=apply((r.r)^2, 2, sum)

    LL[,i]=1012*log(SSEr/SSEf)/4.605

    nn=150
    topset=order(LL[,i], decreasing=T)[1:nn]
    SST=crossprod(r.r[,topset])
    SSE=crossprod(r.f[,topset])
    dSST=determinant(SST,logarithm=T)$modulus
    dSSE=determinant(SSE,logarithm=T)$modulus
    X2=-(1012-1-((nn+2)/2))*(dSSE-dSST)      
    mLOD[i]=X2/4.605
}   
LL[LL<0]=0


switchdet=FALSE
sss=rep(NA, ncol(gdata.s.by.chr[[cc]]))
pb =txtProgressBar(min = 1, max =length(sss), style = 3)
   
#optimize ... SST doesn't need to be recalculated during this scan!!!!! at least 2X speedup 

SST=crossprod((sub.t))
dSST=determinant(SST,logarithm=T)$modulus
if(!is.finite(dSST)){
       t5=svd(SST)
       switchdet=TRUE
}

   for(i in 1:length(sss)){
       setTxtProgressBar(pb, i)
       test=lm(sub.t~gdata.s.by.chr[[cc]][,i])
       # Wilks.test(sub.t, gdata.s.by.chr[[cc]][,i])
       t2=residuals(test)
       #t3=predict(test)
       SSE=crossprod(t2)
       #t(t2)%*%t2
       #SST=crossprod(t2+t3)
       #if(!is.finite(dSSE) | !is.finite(dSST) )  {switchdet=TRUE}
       if(switchdet){
             t4=svd(SSE)
             #t5=svd(SST)
             X2=-(1012-1-((nn+2)/2))*(sum(log(t4$d))-sum(log(t5$d)))
       }else{  
           dSSE=determinant(SSE,logarithm=T)$modulus
           X2=-(1012-1-((nn+2)/2))*(dSSE-dSST)      }
        #log((prod(as.brob(t4$values))/(prod(as.brob(t5$values)))))
        sss[i]=X2/4.605
    #    sss[i]=(-1012*log((prod(as.brob(t4$values))/ (prod(as.brob(t5$values))))))/4.605
    #   print(sss[i])
    }
   close(pb)






LS=apply(LL,2, sum)
plot(LS)


for(i in 101:102){
    print(i)
plot(simL[i,], main=i)
points(LL[,i], col='red')
abline(v=add.qtl.ind)
readline()
}
LLn=do.call(rbind, LLn)
tLL=t(LL)
tLLn=t(LLn)
wm2=apply(LL, 2, which.max)
m2=apply(LL, 2, max)
 plot(simL[5,25:768], tLL[5,])
par(xaxs='i')
plot(wm2[m2>2.5]+25, seq_along(1:(sum(m2>2.5))), xlim=c(0,ncol(gforsim)), col='#00000055',cex=log(m2[m2>2.5],base=4.5))
#     pch='.')



wind=lm(s[,1]~gforsim[,100]+gforsim[,500]-1)
wind.p=lm(s[,1]~gforsim[,100]+gforsim[,120]+gforsim[,500]-1)
wind3=(lm(gforsim[,120]~gforsim[,100]+gforsim[,500]-1))
wind4=lm(residuals(wind)~residuals(wind3)-1)
wind5=lm(residuals(wind)~gforsim[,120]-1)

#wind3=lm(s[,1]~gforsim[,120]-1)

sum(residuals(wind)^2)
sum(residuals(wind.p)^2)
sum(residuals(wind3)^2)
sum(residuals(wind4)^2)
sum(residuals(wind5)^2)

save(hotspots, file='/data/eQTL/RData/hotspots_111615_1600.RData')
   # for(i in 5:210){
  #     kt=names(which(max.obsLOD>i))
  #     length(kt)
  #     sub.t= tmm[,kt]
  #     USVt.k=svd(t(tmm)[kt,])
  #     tmp2= fasterLOD(nrow(tmm),scale(USVt.k$v),gdata.s.by.chr[[cc]])
  #  #x11()
  #     print(paste(i,length(kt),max(apply(tmp2,1,max)))) #max(tmp2[1,])))
  # }
  # i=31

   #nn=(min(300, round(length(keep.transcripts)/10)))
   #nn=min(900,round(length(keep.transcripts)/6))


htable=do.call('rbind', lapply(hotspots, function(x) do.call('rbind', x)))
htable=htable[is.finite(htable[,1]),]
htable=htable[htable$maxLOD<900,]

hpeak.vec=sort(match(as.character(htable$peak), colnames(gdata)))

png(file='/data/eQTL/plots/peaks_genome_fe_batch_new_with_hotspots_chrXIV.png', width=1080, height=1080)
chr.inc=unique.chrs[-17]
chr.inc='chrVIII'
#par(xaxs='i', yaxs='i')
plot(mgi[all.peaks$chr %in% chr.inc],pgi[all.peaks$chr %in% chr.inc], pch=19, cex=all.peaks$LOD[all.peaks$chr %in% chr.inc], type='n', xlab='marker index', ylab='transcript index') 
#, xlim=c(3700,4500))#, col='grey') ##, xlim=c(3750,5500))
#segments(mgiL,pgi, mgiR,pgi, lwd=.5, col='#00000055' )
#, xlab='marker index', ylab='transcript index')#, col='grey') ##, xlim=c(3750,5500))

purp='#80008044'
orng='#FFA50044'
rd='#FF000044'

segments(mgiL[all.peaks$r>0 & all.peaks$chr %in% chr.inc],pgi[all.peaks$r>0 & all.peaks$chr %in% chr.inc], mgiR[all.peaks$r>0  & all.peaks$chr %in% chr.inc],pgi[all.peaks$r>0  & all.peaks$chr %in% chr.inc], lwd=.5, col='#80008055')
segments(mgiL[all.peaks$r<0 & all.peaks$chr %in% chr.inc ],pgi[all.peaks$r<0  & all.peaks$chr %in% chr.inc ], mgiR[all.peaks$r<0 & all.peaks$chr %in% chr.inc],pgi[all.peaks$r<0  & all.peaks$chr %in% chr.inc], lwd=.5, col='#FFA50055')
points(mgi[all.peaks$r>0 & all.peaks$chr %in% chr.inc ],pgi[all.peaks$r>0 & all.peaks$chr %in% chr.inc ], pch=20, cex=log(all.peaks$LOD[all.peaks$chr %in% chr.inc], base=4.5), col=purp)
points(mgi[all.peaks$r<0 & all.peaks$chr %in% chr.inc ],pgi[all.peaks$r<0 & all.peaks$chr %in% chr.inc ], pch=20, cex=log(all.peaks$LOD[all.peaks$chr %in% chr.inc],base=4.5), col=orng)
points(mgi[all.peaks$cis & all.peaks$chr %in% chr.inc ],pgi[all.peaks$cis & all.peaks$chr %in% chr.inc ], pch=20, cex=log(all.peaks$LOD[all.peaks$cis & all.peaks$chr %in% chr.inc],base=4.5), col=rd)
#abline(v=hpeak.vec, lty=3, col='blue')
abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='grey')
abline(h=cumsum(rle(gene.annot.df$chr)$lengths), col='grey')
abline(v=match(sapply(hotspots[[chr.inc]], function(x)x$peak), colnames(gdata)))
dev.off()

png(file='/data/eQTL/plots/hotspot_CI.png', width=512, height=512)
hist(abs(ll-rr), breaks=5000, xlim=c(0,1e4), main='1.5 LOD Intervals for eQTL hotspots', xlab='bp')
dev.off()     


trans.CI.sizes=abs(as.numeric(sapply(strsplit(all.peaks$CI.l[all.peaks$cis==FALSE], '_'), function(x) x[3]))-as.numeric(sapply(strsplit(all.peaks$CI.r[all.peaks$cis==FALSE], '_'), function(x) x[3])))
cis.CI.sizes=abs(as.numeric(sapply(strsplit(all.peaks$CI.l[all.peaks$cis==TRUE], '_'), function(x) x[3]))-as.numeric(sapply(strsplit(all.peaks$CI.r[all.peaks$cis==TRUE], '_'), function(x) x[3])))

png(file='/data/eQTL/plots/QTL_CI.png', width=512, height=512)
par(mfrow=c(2,1))
hist(trans.CI.sizes, breaks=4000, xlim=c(1,.5e5), main='distant QTL 1.5 LOD CIs', xlab='bp')
hist(cis.CI.sizes, breaks=2000, xlim=c(1,.5e5), main='local QTL 1.5 LOD CIs', xlab='bp')
dev.off()










ints=Intervals(cbind(mgiL[all.peaks$r>0 & all.peaks$chr %in% chr.inc], mgiR[all.peaks$r>0  & all.peaks$chr %in% chr.inc]))
png(file='/data/eQTL/plots/all_peaks_as_intervals.png', width=1028, height=512)
par(xaxs='i', yaxs='i')
plot(ints, use_points=FALSE, cex=.5, col='#11000011')
abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='blue')
dev.off()

hotspot.vector=as.vector(unlist(hotspots))
hotspot.vector=hotspot.vector[order(match(hotspot.vector, colnames(gdata)))]

#save(hotspots, file='/data/eQTL/RData/hotspot_markers_new.RData')
load('/data/eQTL/RData/hotspot_markers_new.RData')

hotspot.lms=list()
for(i in 1:ncol(t.tpm.matrix)) {
    print(i)
    transcript=colnames(t.tpm.matrix)[i]
    cQTL=gdata[, closest.marker.to.transcript[i]]
     #yin=t.tpm.matrix[,mid.i]
    if(sum(is.na(cQTL))>0) {    cQTL[is.na(cQTL)]=1  }
    hotspot.lms[[transcript]]=lm(scale(t.tpm.matrix[,i])~gbatch.fact+cQTL+gdata[,hotspot.vector]-1)
}
coefficient.matrix=sapply(hotspot.lms, coefficients)
#save(coefficient.matrix, file='/data/eQTL/RData/hotspot_coefficients_new.RData')


#batch and cis de-regressed phenotypes
bc.resid=list()
for(i in 1:ncol(t.tpm.matrix)) {
    print(i)
    transcript=colnames(t.tpm.matrix)[i]
    cQTL=gdata[, closest.marker.to.transcript[i]]
     #yin=t.tpm.matrix[,mid.i]
    if(sum(is.na(cQTL))>0) {    cQTL[is.na(cQTL)]=1  }
    bc.resid[[transcript]]=scale(residuals(lm(scale(t.tpm.matrix[,i])~gbatch.fact+cQTL-1)))
}

bc.resid=do.call('cbind', bc.resid)

cor.bc.resid=cor(bc.resid)
diag(cor.bc.resid)=0
ots=which(cor.bc.resid>.99, arr.ind=T)
ots=unique(sort(c(ots[,1], ots[,2])))
bc.resid=bc.resid[,-ots]
####HERE

rownames(coefficient.matrix)=gsub('gdata\\[, hotspot.vector\\][0-9]+_chr', '', rownames(coefficient.matrix)) 

gd=hclust(dist(t(coefficient.matrix[15:nrow(coefficient.matrix),])))


par(oma=c(5,1,1,1))
heatmap.2(t(coefficient.matrix[15:nrow(coefficient.matrix),]), col=colorRampPalette(c('orange', 'white', 'purple'))(4),
          density.info='none', trace='none', dendrogram='both', scale=c('row'), keysize=.5)

image.plot(coefficient.matrix[15:nrow(coefficient.matrix),gd$order],col=colorRampPalette(c('orange', 'white', 'purple'))(4))


test=lm(scale(t.tpm.matrix)~gdata[,hotspot.vector[14]]-1)
mt=manova(test)
summary(mt, test='Wilks')
library(Brobdingnag)

pheno.scaled=(scale(residuals.tpm.gbatch.OD))

x=which(test>.99999, arr.ind=T)
prm = pheno.scaled[,-unique(sort(c(x[,1], x[,2])))]
g.hot=gdata[,hotspot.vector]

# experiment with candisc'

g.hot=scale(g.hot)
#142*1000 * 
#g.by.p=t(g.hot)%*% prm
g.by.p=cor((g.hot),(prm))


gc.ind=as.numeric(as.factor(sapply(strsplit(colnames(g.hot), '_'), function(x) x[2])))
abline(v=sapply(1:16, function(x) max(which(gc.ind==x))/141.5), lty=2, col='black')

par(oma=c(4,4,4,4))
heatmap.2((g.by.p), col=redgreen(15), density.info='none', trace='none', dendrogram='column', scale=c('column'), keysize=.2)



#cc.g.p=candisc::cancor(g.hot, prm[,-c(1002,1003,1004)])
#library(CCA)

# ADDITIVE QTL MODEL
lms=list()
for(g in names(peaks.per.gene)){
    print(g)
    peaks.per.gene[[g]]= peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    apeaks = match( peaks.per.gene[[g]]$pmarker, colnames(gdata))

    #newMM.peaks[[i]]
    ax= paste('gdata[,', apeaks,']', sep='')
    aq=paste(ax, collapse= ' + ')
    i=match(g, colnames(pheno.scaled))
    # where pheno.scaled= t.tpm.matrix with batch as a fixed effect
    #am=lm(paste('pheno.scaled[,' , i,']' , ' ~ ', (aq), '-1'))

    am=lm(paste('scale(t.tpm.matrix[,' , i,'])' , ' ~ gbatch.fact + ', (aq), '-1'))
               
                
    #plot(LODS[i,])
    #abline(v=apeaks)

    aov.a = anova(am)
    tssq=sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(am)  
    peaks.per.gene[[g]]$peak.ind=apeaks
    peaks.per.gene[[g]]$var.exp=a.effs[-1]
    peaks.per.gene[[g]]$lm.coeff=as.vector(coeffs)[-c(1:13)]
    #va.a = (tssq-aov.a[nrow(aov.a),2])/(tssq)
    #lms[[g]]=summary(am)$adj.r.squared
    lms[[g]]=scale(residuals(am))
    #summary(am)$adj.r.squared
    #coefficients(am)
}  

#6477
#regress out residual A
A=tcrossprod(gdata)/ncol(gdata)
residual_additiveQTL_pheno=do.call('cbind', lms)
colnames(residual_additiveQTL_pheno)=names(peaks.per.gene)

vc.rA=cbind(rep(NA, ncol(residual_additiveQTL_pheno)), rep(NA, ncol(residual_additiveQTL_pheno)))
rownames(vc.rA)=colnames(residual_additiveQTL_pheno)
eig.rA=doEigenA_forMM(residual_additiveQTL_pheno,A)
for(i in 1:ncol(residual_additiveQTL_pheno)){
        if(is.na(sd(residual_additiveQTL_pheno[,i]))) {
        next;
        }
        vc.rA[i,]=m.S(residual_additiveQTL_pheno[,i], K=A,  theta=eig.rA$theta, Q=eig.rA$Q)
        print(i)
}

b.resid=matrix(NA, nrow(residual_additiveQTL_pheno), ncol(residual_additiveQTL_pheno))
rownames(b.resid)=rownames(residual_additiveQTL_pheno)
colnames(b.resid)=colnames(residual_additiveQTL_pheno)

# b.resid is the phenotype with batch, significant additive QTL, and residual additive variance regressed out
for(i in 1:ncol(residual_additiveQTL_pheno)){
     if(is.na(sd(residual_additiveQTL_pheno[,i]))) {
        next;
        }
    G=(vc.rA[i,1]*A)
    E=vc.rA[i,2]*diag(1012)
    VinV=solve(G+E)
    b.resid[,i]=scale(residual_additiveQTL_pheno[,i]-G%*%VinV%*%residual_additiveQTL_pheno[,i])
    print(i)
}

save.image('/data/eQTL/110915.RData')
colnames(residual_additiveQTL_pheno)=names(peaks.per.gene)

vc.rA2=cbind(rep(NA, ncol( b.resid)), rep(NA, ncol(b.resid)))
rownames(vc.rA2)=colnames(b.resid)
eig.rA2=doEigenA_forMM(b.resid, A)
for(i in 1:ncol(b.resid)){
        if(is.na(sd(b.resid[,i]))) {
        next;
        }
        vc.rA2[i,]=m.S(b.resid[,i], K=A,  theta=eig.rA2$theta, Q=eig.rA2$Q)
        print(i)
}

g.hot=gdata[,hotspot.vector]

g.hot.combos=matrix(NA,1012, 10011)
g.hot.combos.names=rep(NA, 10011)
g.hot.combos.names.ind=rep(NA, 10011)

n=0
for(i in 1:141) {
    print(i)
    for(j in (i+1):142) {
        n=n+1
        g.hot.combos[,n]=g.hot[,i]*g.hot[,j]
        g.hot.combos.names[n]=paste(colnames(g.hot)[i], colnames(g.hot)[j], sep=':')
        g.hot.combos.names.ind[n]=paste(i, j, sep=':')

}}

g.hot.combos.scaled=scale(g.hot.combos)


LODS.2D.hot=fasterLOD(nrow(b.resid),b.resid,g.hot.combos.scaled)
LODS.2D.hot.p1=fasterLOD(nrow(b.resid),b.resid[sample(1:1012),],g.hot.combos.scaled)
LODS.2D.hot.p2=fasterLOD(nrow(b.resid),b.resid[sample(1:1012),],g.hot.combos.scaled)
LODS.2D.hot.p3=fasterLOD(nrow(b.resid),b.resid[sample(1:1012),],g.hot.combos.scaled)
LODS.2D.hot.p4=fasterLOD(nrow(b.resid),b.resid[sample(1:1012),],g.hot.combos.scaled)
LODS.2D.hot.p5=fasterLOD(nrow(b.resid),b.resid[sample(1:1012),],g.hot.combos.scaled)

g.hot.combos=
max.obs=apply(LODS.2D.hot,1,max)
obsPcnt = sapply(seq(2, 7, .01), function(thresh) { sum(max.obs>thresh) }   )
names(obsPcnt) = seq(2, 7, .01)
       

exp.2d.peaks=cbind(apply(LODS.2D.hot.p1,1,max),
                   apply(LODS.2D.hot.p2,1,max),
                   apply(LODS.2D.hot.p3,1,max),
                   apply(LODS.2D.hot.p4,1,max),
                   apply(LODS.2D.hot.p5,1,max))


               
   # expected number of QTL peaks with LOD greater than threshold
   expPcnt = sapply(seq(2, 7, .01),  
                         function(thresh) { 
                                #print(thresh); 
                                mean(apply(exp.2d.peaks, 2, function(ll) {sum(ll>thresh) }) )
                            } )
   names(expPcnt) = seq(2, 7, .01)
   pFDR = expPcnt/obsPcnt
   fdrFX=approxfun(pFDR, seq(2, 7, .01))
   thresh=fdrFX(.1)

   sig.2d=which(LODS.2D.hot>5.61, arr.ind=T)
# first pass, reduce interaction space to 142 hotspot markers
s2h=do.call('rbind', (strsplit(g.hot.combos.names.ind[sig.2d[,2]], ':')))
s2h=cbind(s2h, g.hot.combos.names.ind[sig.2d[,2]])
s2h=cbind(s2h,  names(sig.2d[,1]))
s2h[,1]=as.numeric(s2h[,1])
s2h[,2]=as.numeric(s2h[,2])


s2h.c=rle(sort(s2h[,3]))
s2h.n=do.call('rbind', strsplit(s2h.c$values, ':'))
s2h.n[,1]=as.numeric(s2h.n[,1])
s2h.n[,2]=as.numeric(s2h.n[,2])
png(file='/data/eQTL/plots/2D_142_hotspots.png', width=1080, height=1080)
par(xaxs='i', yaxs='i')
plot(s2h.n[,1], s2h.n[,2], cex=s2h.c$lengths, xlim=c(0,142), ylim=c(0,142), xlab='', ylab='')
abline(v=sapply(1:16, function(x) max(which(gc.ind==x))))
abline(h=sapply(1:16, function(x) max(which(gc.ind==x))))
abline(0,1)
dev.off()

igsig.2d.hotspots=do.call('rbind', (strsplit(g.hot.combos.names.ind[sig.2d[,2]], ':')))












    h2=(vcA[,1]/(vcA[,1]+vcA[,2]))




calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%t(Z)%*%Vinv%*%(y- X%*%B)     }
res



V=residual_additiveQTL_pheno[,1] %*% t(residual_additiveQTL_pheno[,1])




YFR008W

varQTL=list()
#g='YAR029W'

genes.with.cis=sapply(peaks.per.gene, function(x) sum(x$cis))
genes.with.trans=sapply(peaks.per.gene, function(x) sum(!x$cis)>0)
genes.with.trans=sapply(peaks.per.gene, function(x) sum(!x$cis))

genes.with.both=names(which((genes.with.cis+genes.with.trans)>1))
#YDR077W
for(g in genes.with.both[420:length(genes.with.both)]) {
    flag=TRUE
   # names(peaks.per.gene)[1979:length(names(peaks.per.gene))] ){
    print(g)
    peaks.per.gene[[g]]= peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    apeaks = match( peaks.per.gene[[g]]$pmarker, colnames(gdata))
    
    tpeaks = match( peaks.per.gene[[g]]$pmarker[peaks.per.gene[[g]]$cis==FALSE], colnames(gdata))
    cpeaks = match( peaks.per.gene[[g]]$pmarker[peaks.per.gene[[g]]$cis==TRUE], colnames(gdata))

    Qt=(gdata[,tpeaks]%*%t(gdata[,tpeaks])/length(tpeaks))
    Qc=(gdata[,cpeaks]%*%t(gdata[,cpeaks])/length(cpeaks))
    i=match(g, colnames(pheno.scaled))
    
    QtQt=Qt*Qt
    QcQt=Qc*Qt
    QcA=Qc*A
    QtA=Qt*A
    tryCatch( {rQ=regress(pheno.scaled[,i]~1, ~Qc+Qt+QcQt+QtQt, verbose=T)}, error=function(e) {flag=FALSE;})
    if(flag==FALSE) {next;}
    varQTL[[g]]$sigma=rQ$sigma
    varQTL[[g]]$sigma.cov=sqrt(diag(rQ$sigma.cov))
}

vQQ=sapply(varQTL, function(x) x$sigma)
#save(varQTL, file='/data/eQTL/RData/varQTL_AA.RData')
load('/data/eQTL/RData/varQTL_AA.RData')
vQQ.s=vQQ[,-which(vQQ[4,]>vQQ[5,]-.05)]
plot((vQQ.s[,1]),type='l', ylim=c(0,1))
for(i in 2:2512)
points((vQQ.s[,i]),type='l', ylim=c(0,1), col='#00000011')
png(file='/data/eQTL/plots/QTL_VC_model.png', width=1080, height=512)
boxplot(t(vQQ.s)[,-5], ylab='fraction of phenotypic variance', ylim=c(0,1))
dev.off()
#ubpoints(vQQ[4,][vQQ[4,]>vQQ[5,]-.05], vQQ[5,][vQQ[4,]>vQQ[5,]-.05], col='red')

#restart QQ here: YDR370C
AAvar=list()
for(g in names(peaks.per.gene) ){
    print(g)
   
    i=match(g, colnames(pheno.scaled))
    peaks.per.gene[[g]]= peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    apeaks = match( peaks.per.gene[[g]]$pmarker, colnames(gdata))

    Q=(gdata[,apeaks]%*%t(gdata[,apeaks])/length(apeaks))
    QA=Q*A
    QQ=Q*Q
    QQQ=Q*Q*Q
    QQA=QQ*A

    #newMM.peaks[[i]]
    ax= paste('gdata[,', apeaks,']', sep='')
    aq=paste(ax, collapse= ' + ')
        
    am=lm(paste('pheno.scaled[,' , i,']' , ' ~ ', (aq), '-1'))
    im=lm(paste('pheno.scaled[,' , i,']' , ' ~ ', (aq), '-1'))

    #rQA=regress(pheno.scaled[,i]~1, ~A+QA, verbose=T)
    #rAA=regress(pheno.scaled[,i]~1, ~A+AA, verbose=T)
    rQQ=regress(pheno.scaled[,i]~gdata[,apeaks], ~A+QQ, verbose=T)
    #rQQA=regress(pheno.scaled[,i]~1, ~Q+A+QA+QQQ, verbose=T)

    #rQQQ=
    #plot(LODS[i,])
    #abline(v=apeaks)

    #aov.a = anova(am)
    #tssq=sum(aov.a[,2])
    #a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    #coeffs=coefficients(am)  
    #AAvar[[g]]$AA$sigma=rAA$sigma
    #AAvar[[g]]$AA$sigma.cov=sqrt(diag(rAA$sigma.cov))
    #AAvar[[g]]$QA$sigma=rQA$sigma
    #AAvar[[g]]$QA$sigma.cov=sqrt(diag(rQA$sigma.cov))
    
    AAvar[[g]]$QQ$sigma=rQQ$sigma
    AAvar[[g]]$QQ$sigma.cov=sqrt(diag(rQQ$sigma.cov))
 
   # AAvar[[g]]$QQA$sigma=rQQA$sigma
   # AAvar[[g]]$QQA$sigma.cov=sqrt(diag(rQQA$sigma.cov))

   # peaks.per.gene[[g]]$peak.ind=apeaks
   # peaks.per.gene[[g]]$var.exp=a.effs
   # peaks.per.gene[[g]]$lm.coeff=as.vector(coeffs)
    #va.a = (tssq-aov.a[nrow(aov.a),2])/(tssq)
    #lms[[g]]=summary(am)$adj.r.squared
    #coefficients(am)
}    
#ppg=do.call('rbind', peaks.per.gene)

save(AAvar, file='/data/eQTL/RData/AAvar2.RData')

load('/data/eQTL/RData/AAvar2.RData')

jA.A=sapply(AAvar, function(x) x$AA$sigma[1])
jQ.A=sapply(AAvar, function(x) x$QA$sigma[1])
jA.AA=sapply(AAvar, function(x) x$AA$sigma[2])
jQ.AA=sapply(AAvar[names(peaks.per.gene)], function(x) x$QA$sigma[2])
jQ.QQ=as.vector(unlist(sapply(AAvar, function(x) x$QQ$sigma[2])))

png(file='/data/eQTL/plots/QTL_VC_model_AA.png', width=1080, height=1080)
smoothScatter(jA.AA, jQ.AA, xlim=c(-.25,.5), ylim=c(-.25,.5), xlab='AAvariance (genome)', ylab='AAvariance (qtl x genome)', nbin=300, nrpoints=7000)
abline(0,1)
dev.off()
png(file='/data/eQTL/plots/QTL_VC_model_AA_QQ.png', width=1080, height=1080)

smoothScatter(jQ.AA, jQ.QQ, xlim=c(-.1,.5), ylim=c(-.1,.5), xlab='AAvariance (genome x qtl)', ylab='AAvariance (qtl x qtl)', nbin=300, nrpoints=7000)
abline(0,1)
dev.off()


x=do.call('rbind', sapply(AAvar, function(x) c(x$AA$sigma[1], x$QQ$sigma[1])))
png(file='/data/eQTL/plots/A_explained_by_QTL.png', width=1080, height=1080)
par(xaxs='i', yaxs='i')
smoothScatter(x[,1], x[,1]-x[,2], xlim=c(0,1), ylim=c(0,1), xlab='A variance (genome)', ylab='A variance (qtl)', nbin=300, nrpoints=7000)
abline(0,1)
dev.off()

jQQ.AA=sapply(AAvar, function(x) x$QQ$sigma[2])
jQQ.AAA=sapply(AAvar, function(x) x$QQ$sigma[3])


vcE=sapply(lms, function(x) x)
#summary(x)$adj.r.squared)
vcA.with.qtl=vcA.downsample[[1]][match(names(vcE), names(vcA.downsample[[1]]))]

/(vcA[,1][match(names(vcE), rownames(vcA))]+vcA[,2][match(names(vcE), rownames(vcA))])

smoothScatter(vcE, vcA.with.qtl, xlim=c(0,1), ylim=c(0,1), xlab='h^2 (genome)', ylab='h^2 (qtl)', nbin=300)




























load('/data/eQTL/RData/hotspot_coefficients.RData')

hm2=heatmap.2(coefficient.matrix[-c(1:11),])

my_palette=colorRampPalette(c("blue", "white", "red"))(n = 50)

hm2=heatmap.2(t(coefficient.matrix[-c(1:11),]), col=my_palette, 
        density.info="none", trace="none", 
            dendrogram=c("none"), symm=F,symkey=F,symbreaks=T, scale=c("row"), margins=c(20,20), keysize=.5, key=T)



(residuals(lm(t.tpm.matrix[,mid.i]~gbatch.fact+cQTL+gdata[,pmarker])))



library(car)
test=lm(scale(t.tpm.matrix)~gdata[,hotspot.vector[14]]-1)
mt=manova(test)
summary(mt, test='Wilks')
library(Brobdingnag)



t2=residuals(test)
t3=predict(test)
SSE=t(t2)%*%t2
SST=t(t2+t3) %*% (t2+t3)
t4=(eigen(SSE))
t5=eigen(SST)
prod(as.brob(t4$values))/prod(as.brob(t5$values))

#determinant(mt$SSPE, logarithm=TRUE)
#$$modulus
#[1] 2891

y=scale(t.tpm.matrix[,10:600])
determinant(t(y) %*% y, logarithm=TRUE)


2919




2891-2919

2891-



   plot(apply(tmp, 2, sum))
   points(apply(tmp, 2, sum), col='red')
   points(apply(tmp, 2, sum), col='blue')
   points(apply(tmp, 2, sum), col='green')
    
   peak= which.max(apply(tmp, 2, sum))
   maxLOD=rep(NA,100)
   for(i in 1:100) {
    tmp2= fasterLOD(nrow(tmm),tmm[sample(1:nrow(tmm)),],gdata.s.by.chr[[cc]])
    maxLOD[i]=(max(apply(tmp2, 2, sum)))
    print(i)
   }
   thresh=quantile(maxLOD, .95)
    
   pmarker=c(pmarker, names(peak))




}


tmp2= fasterLOD(nrow(tmm),tmm[sample(1:nrow(tmm)),],gdata.s.by.chr[[cc]])




test=Manova(lm(tmm~gdata[,100]))


 


       else {
     tmm[,i]=scale(residuals(lm(t.tpm.matrix[,i]~gbatch.fact+gdata[,bQTL]+gdata[,cQTL])))
     }                    
                              
                              background.QTL[[cc]][[i]]]+gdata[,found.q.covs[[i]]])))
}

trans3=matri
scale(residuals(lm(t.tpm.matrix[,i]~gbatch.fact+gdata[,background.QTL[[cc]][[i]]]+gdata[,found.q.covs[[i]]]))



offCnt=colSums(count.matrix)
XF=covariates.OD
XQ=gdata[,peaks.per.gene[['YAL062W']]$pmarker]
lresid=residuals(lm(t.tpm.matrix[,'YAL062W']~XF))
summary(lm(t.tpm.matrix[,'YAL062W']~XQ*XF-1))

covariates.OD+gdata[,peaks.per.gene[['YAL062W']]$pmarker]))

summary(glm(count.matrix['YAL062W',]~-1+XF*XQ+offset(log(offCnt)), family='poisson'))

plot(predict(glm(count.matrix['YAL062W',]~-1+XF*XQ+offset(log(offCnt)), family='poisson')), count.matrix['YAL062W',])

XQ2=XQ*XQ
plot(predict(lm(t.tpm.matrix[,'YAL062W']~XQ+XQ2+XF-1)), t.tpm.matrix[,'YAL062W'])




summary(glmer(count.matrix['YAL062W',]~-1+XF+XQ+offset(log(offCnt)), family='poisson'))






plot(mgi[all.peaks$cis==FALSE], (all.peaks$LOD[all.peaks$cis==FALSE]), cex=.5)

# random effect model with batch correction, without OD correction
#FDR of 5% genome wide
#2,259 cis
#21,677  trans
# smith
# 3809 transe
#1200 cis

# fixed effect model with OD correction and batch correction
#FDR of 5% genome wide
#2262 cis
#19622 trans

# fixed effect model (without OD correction) and with batch correction
#FDR of 5% genome wide
#2192 cis
#sum(!all.peaks$cis)
#[1] 21583












pcis.cnt=sapply(peaks.per.gene, function(x){sum(x$cis)})
ptrans.cnt=sapply(peaks.per.gene, function(x){sum(!x$cis)})
qtl.cnt=sapply(peaks.per.gene, function(x){nrow(x)})

png(file='/data/eQTL/plots/QTL_cnt.png', width=768, height=768)
hist(qtl.cnt, breaks=70, xlab='QTL count per transcript' ,xaxt='n', ylab='transcripts', main='QTL count')
axis(side = 1, at = c(0,1,2,3,4,5,6,7,8,9,10,15,20), labels=c(0,1,2,3,4,5,6,7,8,9, 10,15,20) )
dev.off()


#plot(density(ppg$var.exp[!ppg$cis]), col='blue')
#points(density(ppg$var.exp[ppg$cis]), col='red', type='l')
png(file='/data/eQTL/plots/eff_size_dist_new.png', width=768, height=768)
da=density(ppg$var.exp[!ppg$cis])
di=density(ppg$var.exp[ppg$cis])
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(da, xlab='fraction of phenotypic variance', ylab='density', xlim=c(0,0.5),  type='n', main='')
polygon(da, col=rgb(0, 0, 1,0.3), border=NA)
polygon(di, col=rgb(1, 0, 0,0.3), border=NA)
legend('topright', legend=c('distant', 'local'), fill=c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)))
dev.off()


hist(ppg$var.exp[!ppg$cis], col=rgb(0,0,1,0.5), border=rgb(0,0,1,0.5),  xlim=c(0,1), ylim=c(0,2000), breaks=50)
hist(ppg$var.exp[ppg$cis], col=rgb(1,0,0,0.5),border=rgb(1,0,0,0.5), add=T, breaks=250)


hist(ppg$var.exp[!ppg$cis], col=rgb(0,0,1,0.5), border=rgb(0,0,1,0.5),  xlim=c(0,1), ylim=c(0,2000), breaks=500, add=T)


   plot(fasterLOD(797, pheno.scaled[,g], gdata.scaled)[1,], cex=.5)

   scale(residuals(lm(pheno.scaled[,g]~gdata[,gind[-1][-1]])))
   points(fasterLOD(797, scale(residuals(lm(pheno.scaled[,g]~gdata[,gind[-1][-1]]))), gdata.scaled)[1,], col='red', cex=.5)


   points( fasterLOD(797, pheno.scaled[,g], gdata.scaled)[1,]


    lms[[g]]=anova(lm(pheno.scaled[,g]~gdata[,gind]-1))
}

A=tcrossprod(gdata)/ncol(gdata)
vcA=cbind(rep(NA, ncol(pheno.scaled)), rep(NA, ncol(pheno.scaled)))
rownames(vcA)=colnames(pheno.scaled)
eigA=doEigenA_forMM(pheno.scaled,A)
# calculate mixed model, one term for additive variance  -------------------------------------------
for(i in 1:ncol(pheno.scaled)){
        if(is.na(sd(pheno.scaled[,i]))) {
        next;
        }
        vcA[i,]=m.S(pheno.scaled[,i], K=A,  theta=eigA$theta, Q=eigA$Q)
     print(i)
}
   # save(vcA, file=paste0(base.dir, 'RData/vcA.', pset, '.RData'))
#------------------------------------------------------------------------------------------------------





vcE=sapply(lms, function(x) summary(x)$adj.r.squared)
vcA.with.qtl=vcA[,1][match(names(vcE), rownames(vcA))]/(vcA[,1][match(names(vcE), rownames(vcA))]+vcA[,2][match(names(vcE), rownames(vcA))])
png(file='/data/eQTL/plots/so_much_missing_h2.png',width=500,height=500)
smoothScatter(vcE, vcA.with.qtl, xlim=c(0,1), ylim=c(0,1), xlab='h^2 (genome)', ylab='h^2 (qtl)', nbin=300)
points(vcE, vcA.with.qtl,pch='.')
abline(0,1)
dev.off()



str(pheno.scaled)

t(pheno.scaled)





    LODSb=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled, betas=T, pheno=pheno.scaled)
   # png(file='~/Desktop/peff.png' ,width=4096,height=4096)
   # k=pheatmap(LODSb$beta, cluster_rows=TRUE, cluster_cols=FALSE)
   # dev.off()
   # x=c(2,10,50,100,150,200,250,300)
   # y=c(207458, 157525, 130440,118162,111160,106021,101437,97754)
    kmean=list()
    for(k in c(2,10,50,100,150,200,250,300) ){
    print(k)
    k2=kmeans(LODSb$beta,centers=k)
    n=as.character(k)
    kmean[[n]]=k2
    print((k2$tot.withinss))
    }
    plot(kmean[['300']]$centers[2,])
        abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='grey')

    #str(LODSb$beta)
    #library(fastcluster)
    #dmat=dist(LODSb$beta)
    #dmat[is.na(dmat)]=0
    #  x=(hclust(dmat))

    #save(LODS, file=paste0(base.dir, 'RData/LODS.', pset, '.RData'))

    wl3=which(LODS>3, arr.ind=T)
    wl10=which(LODS>10, arr.ind=T)
    wl100=which(LODS>100, arr.ind=T)

    #rcolorbrewer red to blue with white in middle
    # Visualization -------------------------------------------------------------------------------------------------------------------------------
    png(file=paste0(base.dir, 'plots/BYxRM_eQTL_', pset, '.png'), width=1080, height=1080)
    par(xaxs='i', yaxs='i')
    plot(wl3[,2],wl3[,1], xlab='marker position (marker genome index)', ylab='gene position (gene genome index)', pch=19, col='#00000022', cex=.5,
         main=paste('n=', nrow(pheno.scaled), '   BYxRM eQTL'))
    points(wl10[,2],wl10[,1], pch=19, col='#0000FF22', cex=.5)
    points(wl100[,2],wl100[,1], pch=19, col='#FF000022', cex=.5)
    abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='grey')
    abline(h=cumsum(rle(gene.annot.df$chr)$lengths), col='grey')
    dev.off()
    # --------------------------------------------------------------------------------------------------------------------------------------------

    # crude hotspot plot ------------------------------------------------------------------------------
    png(file=paste0(base.dir, 'plots/BYxRM_eQTL_hotspots_', pset, '.png'), width=1080, height=1080)
    hist(wl3[,2], breaks=100000, main=pset)
    abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='grey')
    dev.off()
    #--------------------------------------------------------------------------------------------------



    png(file=paste0(base.dir, 'plots/hist_vcA_', pset, '.png'), width=1080, height=1080)
    hist(vcA, breaks=100, xlab='h2', ylab='Transcripts', 
         main=paste('histogram of h^2  for', pset), sub=paste('median h^2 =', round(median(vcA, na.rm=T) ,digits=3) ))
    dev.off()
    
    png(file=paste0(base.dir, 'plots/vcA_vs_abundance_', pset, '.png'), width=1080, height=1080)
    plot(log10(apply(count.matrix, 1, sum)), vcA, col='#00000066', xlab='log10 (total reads per gene across samples)', ylab='h^2')
    dev.off()
}

pset.vec=c(
 'downsample.1e4.tpm.bc',
 'downsample.5e4.tpm.bc',
 'downsample.1e5.tpm.bc',
 'downsample.5e5.tpm.bc',
 'tpm.bc'
)

vcDF=list()
for(pset in pset.vec) {
    load(file=paste0(base.dir, 'RData/vcA.', pset, '.RData'))
    print(i)
    cisModel[i]=coef(lm(pheno[,i]~gdata[,closest.marker.to.transcript[i]]))[2]
}


6685

cisLODS=LODS[cbind(which(!is.na(closest.marker.to.transcript)),closest.marker.to.transcript[which(!is.na(closest.marker.to.transcript))])]
sum(cisLODS>2.5, na.rm=T)                  

gene.GR2=gene.GR[1:6685]
gene.GR2$cisLOD=cisLODS
plot(start(gene.GR2[seqnames(gene.GR2)=='chrVII']         ), gene.GR2[seqnames(gene.GR2)=='chrVII']$cisLOD)
sig7=which(gene.GR2[seqnames(gene.GR2)=='chrVII']$cisLOD>100)
gene.GR2[seqnames(gene.GR2)=='chrVII'][sig7]


plot(log10(apply(count.matrix, 1, sum))[which(!is.na(closest.marker.to.transcript))], cisLODS)
cor.test(log10(apply(count.matrix, 1, sum))[which(!is.na(closest.marker.to.transcript))], cisLODS)
sum(cisLODS>2.0, na.rm=T)








#library("DESeq2")
#count.ints=(apply(count.matrix, 2, as.integer))
#rownames(count.ints)=rownames(count.matrix)
#ddsMat=DESeqDataSetFromMatrix(countData=count.ints, colData=sample.annot.df, design=~1)
#rld=rlog(ddsMat)
#tpm.matrix=apply(count.matrix,2, function(x) countToTpm(x, gene.annot.df$eff_length))
#se = SummarizedExperiment(scale(tpm.matrix),
#                                colData=colData(ddsMat))
#plotPCA(DESeqTransform( se ), intgroup='growth.batch', ntop=500)
#pout=plotPCA(DESeqTransform( se ), intgroup='growth.batch', ntop=500, returnData=T)
#plot(pout$PC1, sample.sum)
#growth.batch' )




