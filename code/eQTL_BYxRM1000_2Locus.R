# code to compute exhaustive 2 locus scan (and marginal scan)
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
library(foreach)
library(doMC)
registerDoMC(cores=70)
base.dir='/data/eQTL/'

# load funcitons 
source(paste0(base.dir, 'code/eQTL_BYxRM1000_fx_stranded.R'))
unique.chrs=c(paste0('chr', as.roman(1:16)), 'chrMito')

# load genotype data
load('/data/eQTL/RData/gdata.RData')
gdata.scaled=scale(gdata)

# load expression data 
load('/data/eQTL/RData/log2_t.tpm.matrix.RData')

# load expression covariates 
load('/data/eQTL/RData/covariates.OD.RData')

# load marker and transcript annotation data
load('/data/eQTL/RData/markerAnnotation.RData')
load('/data/eQTL/RData/geneAnnotation.RData')

cvec=(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])
gdata.s.by.chr=list()
gdata.by.chr=list()
for(cc in unique.chrs[1:16]) {   
   gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
   gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]  
}

genetic.map=buildGeneticMap(gdata.by.chr)
genetic.map.c=as.vector(unlist(genetic.map))
names(genetic.map.c)=as.vector(unlist(sapply(genetic.map, names)))

# do 2 locus scan --------------begin by pruning marker data by LD  ---------------------------------------------------------------------------------------------------------
gdata.downsampled=downsampleMarkers(gdata.by.chr, gdata.s.by.chr)

#save(gdata, file='/data/eQTL/RData/gdata.RData')
#load('/data/eQTL/RData/2D/gdata.RData')
#save(gdata.downsampled, file='/data/eQTL/RData/gdata.downsampled.RData')

# residualize trait (remove additive effects) (note that output here is not scaled)
pheno.additive.removed=removeAdditiveEffects(t.tpm.matrix, peaks.per.gene, covariates.OD, gdata)
#save(pheno.additive.removed, file='/data/eQTL/RData/pheno.additive.removed.RData')

load('/data/eQTL/RData/2D/pheno.additive.removed.RData')

pheno.additive.removed.scaled=scale(pheno.additive.removed)
# sanity check that additive signal is gone
#scanoneLODS.ar=fasterLOD(nrow(pheno.additive.removed.scaled),pheno.additive.removed.scaled,gdata.scaled, betas=TRUE)

marker.GR.downsampled=marker.GR[match(as.vector(unlist(sapply(gdata.downsampled, colnames))), marker.GR$mname),]
all.peaks.DS=all.peaks.OD
all.peaks.DS$marker.ds.ind=findInterval(all.peaks.OD$marker.gcoord, marker.GR.downsampled$gcoord)
all.peaks.DS$pmarker.ds=marker.GR.downsampled$mname[all.peaks.DS$marker.ds.ind]

#save(all.peaks.DS, file='/data/eQTL/RData/all.peaks.DS.RData')
load('/data/eQTL/RData/all.peaks.DS.RData')

obs2D=lapply(do2locusScan.reduced_output(pheno.additive.removed.scaled, 
                                         gdata.downsampled, all.peaks.DS, marker.gap=20), collapse2D)
set.seed(1)
perm2D=list()
for(i in 1:5) {
    perm2D[[as.character(i)]]=lapply(do2locusScan.reduced_output(pheno.additive.removed.scaled[sample(1:nrow(pheno.additive.removed.scaled)),], 
                                                                 gdata.downsampled, all.peaks.DS, marker.gap=20), collapse2D)
}
saveRDS(obs2D , file= '/data/eQTL/RData/obs2D.peaks.RDS')
saveRDS(perm2D, file='/data/eQTL/RData/perm2D.peaks.RDS')

obscnt=lapply(obs2D, function(o){ 
                  v=sapply(seq(2,8,.05), function(x) sum(o$lod>x)) 
                  names(v)=seq(2,8,.05)
                  return(v)    
                      } )
permcnt=lapply(perm2D, function(y){ 
          lapply(y, function(o){ 
                  v=sapply(seq(2,8,.05), function(x) sum(o$lod>x)) 
                  names(v)=seq(2,8,.05)
                  return(v)    
                      } ) })


permcnt.full=rowMeans(sapply(permcnt, function(x) x[[1]]))
permcnt.marginal=rowMeans(sapply(permcnt, function(x) x[[2]]))

permcnt.full/obscnt$full.scan.peaks
# 7.5 5%
# 7.15 10%

permcnt.marginal/obscnt$marginal.scan.peaks
#5.35 5% 
#4.9  10%

full2d.table=obs2D$full.scan.peaks[obs2D$full.scan.peaks$lod>7.15,]
write.table(full2d.table, file='/data/eQTL/RData/full2D_table.txt', sep='\t', quote=F, row.names=F)

marginal2d.table=obs2D$marginal.scan.peaks[obs2D$marginal.scan.peaks$lod>4.9,]
write.table(marginal2d.table, file='/data/eQTL/RData/marginal2D_table.txt', sep='\t', quote=F, row.names=F)

#Test this ... fast 2D scan for epistasis chromsome by chromosome
set.seed(10)
for(pp in 1:6) {
    VC_2Dp=list()
    for(cc1 in 1:16) {
            for(cc2 in cc1:16) {
                print(paste(cc1, cc2))
                c1=gdata.s.by.chr[[unique.chrs[cc1]]]
                c2=gdata.s.by.chr[[unique.chrs[cc2]]]
                A1=A.mat(c1)/2
                A2=A.mat(c2)/2
                VC2=A1*A2
                if(pp==1) {  s2vc=calcA(pheno.additive.removed.scaled, VC2, do.print=F) } else{
                    s2vc=calcA(pheno.additive.removed.scaled[sample(1:1012),], VC2, do.print=F)
                }
                VC_2Dp[[paste0(cc1, '_', cc2) ]]=s2vc
            }
    }
    #save(VC_2D, file= '/data/eQTL/RData/stranded/VC_2D.RData')
    saveRDS(VC_2Dp, file= paste0('/data/eQTL/RData/VC_2D.RDS.', pp))
}
# process this output 



v2=sapply(VC_2D, function(x) sum(x[,1]>.03))
v3=do.call('rbind', strsplit(names(v2), '_'))
v4=data.frame(v3,v2, stringsAsFactors=F)
plot(v4[,1], v4[,2], type='n', ylim=c(1,16), xlim=c(1,16))
text(as.numeric(v4[,1]), as.numeric(v4[,2]),v4[,3], cex=log10(v4[,3]/10))
#------------------------------------------------------------------------------------------------------------------------


# Visualizing results of 2 locus scan 
mm2d1=marginal2d.table$m1[marginal2d.table$lod>4.9]
mm2d2=marginal2d.table$m2[marginal2d.table$lod>4.9]
mm2d1.ind=marker.GR$gcoord[match(mm2d1,colnames(gdata))]
mm2d2.ind=marker.GR$gcoord[match(mm2d2,colnames(gdata))]

bins2d=seq(1,1.2e7,99000)

b1=cut(mm2d1.ind, bins2d)
b2=cut(mm2d2.ind, bins2d)

b3=cbind(as.numeric(b1), as.numeric(b2))
b3=na.omit(b3)
b3=b3[order(b3[,1]),]

b4=paste(b3[,1], b3[,2], sep='_')
b5=rle(b4)

b6=data.frame( do.call('rbind', strsplit(b5$values, '_')), b5$lengths, stringsAsFactors=F)
b6[,1]=as.numeric(b6[,1])
b6[,2]=as.numeric(b6[,2])
par(xaxs='i', yaxs='i')
plot(bins2d[b6[,1]],bins2d[b6[,2]], type='n', ylim=c(0, max(gcoord.key)), xlim=c(0, max(gcoord.key)),
     xaxt='n', yaxt='n', xlab='', ylab='') 
text(bins2d[b6[,1]],bins2d[ b6[,2]], b6[,3], cex=log(b6[,3],7)+.5, col=ifelse(b6[,3]>10, 'red', 'black'))
abline(v=gcoord.key, lty=1, col='grey')
abline(h=gcoord.key, lty=1, col='grey')
abline(0,1)
#abline(v=marker.GR$gcoord[hp.index], col='blue', lty=2)
#abline(h=marker.GR$gcoord[hp.index], col='blue', lty=2)

plot(mm2d1.ind, mm2d2.ind, col='#00ff0022', pch=19 , cex=2)
points(mm2d1.ind, mm2d2.ind, col='#0000ff22', pch=19, cex=2)

plot(hexbin(mm2d1.ind, mm2d2.ind))
smoothScatter(mm2d1.ind, mm2d2.ind)
abline(v=gcoord.key, lty=1, col='black')

abline(v=hp.index, lty=2, col='blue')
abline(h=hp.index, lty=2, col='blue')
#--------------------------------------------------------------------------------------------------------------------------


# targetted testing for 2 and 3 locus interactions 

lm.2Dint=list()
for(g in names(peaks.per.gene)){
    print(g)
    yr=residuals(lm(t.tpm.matrix[,g]~covariates.OD))
    pms=peaks.per.gene[[g]]$pmarker[!duplicated(peaks.per.gene[[g]]$pmarker)]
    if(length(pms)>1) {
     X=data.frame(gdata[,pms])
       fitme=lm(yr~., data=X)

     lm.2Dint[[g]]=add1(fitme, ~.^2, test='F')
    } else {next;}
 }
ps.QTL.2D=unlist(sapply(lm.2Dint, function(x) x$'Pr(>F)'))
qs2D =qvalue(na.omit(ps.QTL.2D))
qthresh2d=max(qs2D$qvalues[qs2D$qvalues<.1]) 


#lm.3Dint=list()
#for(g in names(peaks.per.gene)[3415:length(names(peaks.per.gene))] ) {
#    print(g)
#    yr=residuals(lm(t.tpm.matrix[,g]~covariates.OD))
#    
#    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
#    pms=ppg$pmarker
#    if(length(pms)>14) {
#         pms=peaks.per.gene[[g]]$pmarker[!duplicated(peaks.per.gene[[g]]$pmarker) & abs(peaks.per.gene[[g]]$r)>.14]
#    }
#    if(length(pms)>14) {
#         pms=peaks.per.gene[[g]]$pmarker[!duplicated(peaks.per.gene[[g]]$pmarker) & abs(peaks.per.gene[[g]]$r)>.17]
#    }
#
#
#    if(length(pms)>2) {
#      X=data.frame(gdata[,pms])
#      fitme=lm(yr~., data=X)
#      a2=add1(fitme, ~.^2, test='F')
#      nterms=rownames(a2[which(a2[,6]< qthresh2d),])
#      
#      #fit2=update(fitme, as.formula(paste0('~.+',paste(nterms, collapse=' + '))))
#      fitme2=lm(yr~.^2, data=X)
#      a3=add1(fitme2, ~.^3, test='F')
#     lm.3Dint[[g]]=a3
#    } else {next;}
# }
#ps.QTL.3D=unlist(sapply(lm.3Dint, function(x) x$'Pr(>F)'))
#qs3D =qvalue(na.omit(ps.QTL.3D))
#qthresh3d=max(qs3D$qvalues[qs3D$qvalues<.1]) 

#covariates.OD=model.matrix(t.tpm.matrix[,1]~gbatch.fact+OD.cov)
#residual.pheno.OD=residuals(lm(t.tpm.matrix~covariates.OD) )
#pheno.scaled.OD=(scale(residual.pheno.OD))
ps.QTL.3D.f=unlist(sapply(lm.3Dint2, function(x) x$'Pr(>F)'))
qs3Df =qvalue(na.omit(ps.QTL.3D.f))


#almost 2X more QTL at FDR 10% for targetted 2D scan
# top candidates for 3x interaction
g='YPL277C'
g='YFL059W'
g='YOR142W-B'
g='YPL278C'
g='YDR365W-B'
g='YDR034C'
g='YJL186W'
g='YPL277C'
yr=residuals(lm(t.tpm.matrix[,g]~covariates.OD))
    
    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    pms=ppg$pmarker
    max.peak=ppg$pmarker[which.max(ppg$LOD)]
       X=data.frame(gdata[,pms])

fitme=lm(yr~.^3, data=X)
summary(fitme)

