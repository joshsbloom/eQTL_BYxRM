# code to detect QTL-QTL interactions 
# full, marginal, and additive only scans
# also 2D variance component analysis

# load expression data 
# load('/data/eQTL/RData/log2_t.tpm.matrix.RData')
# load expression covariates 
# load('/data/eQTL/RData/covariates.OD.RData')

# load marker and transcript annotation data
# load('/data/eQTL/RData/markerAnnotation.RData')
# load('/data/eQTL/RData/geneAnnotation.RData')


# do 2 locus scan --------------begin by pruning marker data by LD  ---------------------------------------------------------------------------------------------------------
# load('/data/eQTL/RData/gdata.downsampled.RData')
# load('/data/eQTL/RData/pheno.additive.removed.RData')

pheno.additive.removed.scaled=scale(pheno.additive.removed)
# sanity check that additive signal is gone
#scanoneLODS.ar=fasterLOD(nrow(pheno.additive.removed.scaled),pheno.additive.removed.scaled,gdata.scaled, betas=TRUE)

AA=(tcrossprod(gdata.scaled)/ncol(gdata.scaled))^2
## this should be about the same as the 
vcAA=calcA(pheno.additive.removed.scaled, AA)
#h2AA=vcAA[,1]/(vcAA[,1]+vcAA[,2])

set.seed(10)
yAA.perm=replicate(10, {pheno.additive.removed.scaled[sample(1:1012),]})
WG.2D.VC.perm=sapply(1:10, function(p) {
    vcAA=calcA(yAA.perm[,,p], AA)
    return(vcAA[,1]) #/(vcAA[,1]+vcAA[,2]))
})

o2dvc=sapply(seq(0,.5,.001), function(x) sum(h2AA > x ) )
names(o2dvc)=seq(0,.5,.001)
e2dvc=rowMeans(apply(WG.2D.VC.perm, 2, function(y) { sapply(seq(0,.5,.001) , function(x) sum( y > x ) )     }))
names(e2dvc)= seq(0,.5,.001)

e2dvc/o2dvc
colnames(t.tpm.matrix)[which(h2AA>.4)]
#FDR 10%
# h2AA>0.39 

#load('/data/eQTL/RData/all_peaks_OD.RData')
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
#saveRDS(obs2D , file= '/data/eQTL/RData/obs2D.peaks.RDS')
#saveRDS(perm2D, file='/data/eQTL/RData/perm2D.peaks.RDS')

obs2D=readRDS('/data/eQTL/RData/obs2D.peaks.RDS')
perm2D=readRDS('/data/eQTL/RData/perm2D.peaks.RDS')

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
#387 (306 transcripts)
#write.table(full2d.table, file='/data/eQTL/RData/full2D_table.txt', sep='\t', quote=F, row.names=F)

marginal2d.table=obs2D$marginal.scan.peaks[obs2D$marginal.scan.peaks$lod>4.9,]
#784 (590 transcripts)
#write.table(marginal2d.table, file='/data/eQTL/RData/marginal2D_table.txt', sep='\t', quote=F, row.names=F)

# For 2D model that is only testing interactions between additive QTL ------------------
set.seed(1)
registerDoMC(cores=11)
#yr=residuals(lm(t.tpm.matrix~covariates.OD))
yr=pheno.additive.removed.scaled
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
####old FDR 5%  = 0.0003548    10^(3.45)
####old  FDR 10% = 0.001        10^(3)
qthresh2d=10^(-2.9) #.001


#peaks.per.gene.augmented=list()
#Ahot=tcrossprod(gdata.scaled[,sort(match(as.character(htabler[,
# split this up into separate chunks 09 16 16

peaksModel=list()
lm.2D=list()
lm.3D=list()
pb =txtProgressBar(min = 1, max = length(peaks.per.gene), style = 3)
for(gg in 1:length(peaks.per.gene)) {
    setTxtProgressBar(pb, gg)
    g=names(peaks.per.gene)[gg]
    #in names(peaks.per.gene) ) {
    #print(g)
    yt=t.tpm.matrix[,g]
    yr=pheno.scaled.OD[,g]
    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    peaksModel[[g]]=ppg
    
    apeaks = match( ppg$pmarker, colnames(gdata))
    X=data.frame(covariates.OD, gdata[,apeaks])

    fitme=lm(yt~.-1, data=X)
    aov.a = anova(fitme)
    tssq  = sum(aov.a[-c(1:14), 2])
    a.effs=aov.a[15:(nrow(aov.a)-1),2]/tssq # (aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(fitme) 
    peaksModel[[g]]$var.exp.Raw=a.effs
    peaksModel[[g]]$lm.coeff.Raw=as.vector(coeffs)[-c(1:14)]

    X2=data.frame(gdata[,apeaks])
    fitme=lm(yr~.-1, data=X2)
    #X2=data.frame(covariates.OD, gdata[,pms])
    #fitme2=lm(t.tpm.matrix[,g]~., data=X2)
    # for additive model ---------------------------------------------
    aov.a = anova(fitme)
    tssq  = sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    coeffs=coefficients(fitme)  
    peaksModel[[g]]$var.exp.Resid=a.effs
    peaksModel[[g]]$lm.coeff.Resid=as.vector(coeffs)

     if(length(ppg$pmarker)>1) {
            a2=add1(fitme, ~.^2, test='F')
            nterms=rownames(a2[which(a2[,6]< qthresh2d),])
            if(length(nterms) > 0)  {
                fit2=update(fitme, as.formula(paste0('~.+',paste(nterms, collapse=' + '))))
                lm.2D[[g]]=fit2
                ints=rownames(a2)[-1]
                max.peak=ppg$pmarker[which.max(ppg$LOD)]
                lm.3D[[g]]=add1(fit2,  paste0(gsub(':|/', '.', max.peak), ':', ints), test='F')
                }
     }
}
close(pb)
#save(peaksModel, file='/data/eQTL/RData/peaksModel.RData')
#load('/data/eQTL/RData/peaksModel.RData')

lm.2D.count=sapply(lm.2D, function(x) {  sum(grepl(':', names(coef(x))))   })
sum(lm.2D.count)
[1] 1464

ps.3D=na.omit(unlist(sapply(lm.3D, function(x) x[,6])))
#sum(qvalue(ps.3D)$qvalue<.1)
# 16
q3d=qvalue(ps.3D)
q3dthresh=max(q3d$pvalues[q3d$qvalues<.1])

sapply(lm.3D.test, function(x) x
l3dthresh=lapply(lm.3D.test, function(x) {
                     ret=x[,6]<q3dthresh
                     if(sum(ret,na.rm=T)>0 ) {
                     return(x[ret,])
                     } else{return(NULL) } })

l3dthresh[which(sapply(l3dthresh, length)>0)]
t.tpm.matrix[,'YFL059W']

boxplot(pheno.scaled.OD[,'YFL059W']~gdata[,'chrVI:30187_G/A']+gdata[,'chrX:722130_T/C']+gdata[,'chrXIV:18838_C/T'])
#SNZ3 CIS
#SNZ2 trans (chrXIV)
# PGU1
#boxplot(pheno.scaled.OD[,'YNL333W']~gdata[,'chrVI:30187_G/A']+gdata[,'chrX:722130_T/C']+gdata[,'chrXIV:18838_C/T'])

#YPL277C CIS
#boxplot(pheno.scaled.OD[,'YPL277C']~gdata[,'chrXVI:26510_T/C']+gdata[,'chrXV:1065815_A/G']+gdata[,'chrXIII:914575_T/C'])
#chrXIII FET4
#PHR1 or FRE5

# restructure 2D interactions 
lm.2D.df=list()
for(goi in names(lm.2D)) {
   print(goi) 
    #goi='YAL039C'
    #lm.2D[[goi]]
    ppga=peaks.per.gene[[goi]]
    lm.goi=lm.2D[[goi]]
    aov.a = anova(lm.goi)
    tssq  = sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    c2d=coefficients(lm.goi)  
    c2d=c2d[grep(':', names(c2d))]
    vexp2d=a.effs[grep(':', rownames(aov.a))]
    c2d.df=do.call('rbind', strsplit(names(c2d), ':'))
    pmarker.mod=gsub('\\:|\\/', '.', ppga$pmarker)

    m1.ind=match(c2d.df[,1], pmarker.mod)
    m2.ind=match(c2d.df[,2], pmarker.mod)

    lm.2D.df[[goi]]=data.frame(
               coef2D=as.numeric(c2d),
               vexp2D=as.numeric(vexp2d),
               m1.pmarker=ppga$pmarker[m1.ind],
               m2.pmarker=ppga$pmarker[m2.ind],
               m1.gcoord=ppga$marker.gcoord[m1.ind],
               m2.gcoord=ppga$marker.gcoord[m2.ind],
               t.gcoord=ppga$gene.gcoord[m1.ind],
               partner.cis=ppga$cis[m1.ind]+ppga$cis[m2.ind],
               t.regulated=ppga$gene[m1.ind]
               )
}
lm.2D.df.long=do.call('rbind', lm.2D.df)



# for supplementary table 040818
lm.2D.2=list()
for(goi in names(lm.2D)) {
   print(goi) 
    #goi='YAL039C'
    #lm.2D[[goi]]
    ppga=peaks.per.gene[[goi]]
    lm.goi=lm.2D[[goi]]
    aov.a = anova(lm.goi)
    tssq  = sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    c2d=coefficients(lm.goi)  
    c2d=c2d[grep(':', names(c2d))]
    vexp2d=a.effs[grep(':', rownames(aov.a))]
    c2d.df=do.call('rbind', strsplit(names(c2d), ':'))
    pmarker.mod=gsub('\\:|\\/', '.', ppga$pmarker)

    m1.ind=match(c2d.df[,1], pmarker.mod)
    m2.ind=match(c2d.df[,2], pmarker.mod)
    csl=coefficients(summary(lm.goi))
    csli=grep(':', rownames(csl))
    lm.2D.2[[goi]]=data.frame(
               trait=ppga$gene[m1.ind],
               m1_pmarker=ppga$pmarker[m1.ind],
               m1_beta=coefficients(lm.goi)[m1.ind],
               m1_variance_explained=a.effs[m1.ind],
               m1_p=csl[m1.ind,4],
               m1_is_cis=ppga$cis[m1.ind],
               m2_pmarker=ppga$pmarker[m2.ind],
               m2_beta=coefficients(lm.goi)[m2.ind],
               m2_variance_explained=a.effs[m2.ind],
               m2_p=csl[m2.ind,4],
               m2_is_cis=ppga$cis[m2.ind],
               interaction_beta=as.numeric(c2d),
               interaction_variance_explained=as.numeric(vexp2d),
               interaction_p=csl[csli,4]
               )
}
lm2D_between_additive=do.call('rbind', lm.2D.2)

save(lm2D_between_additive,file='/data/eQTL/RData/lm2D_between_additive.RData')

m1=data.table::tstrsplit(full2d.table$m1, '_|:', type.convert=T)
m2=data.table::tstrsplit(full2d.table$m2, '_|:', type.convert=T)
names(m1)=c('chr','start','var')
names(m2)=c('chr','start','var')
m1$end=m1$start
m2$end=m2$start

f2d.tchr=as.character(seqnames(gene.GR)[match(full2d.table$trait, gene.GR$ORF)])
f2d.tpos=start(gene.GR)[match(full2d.table$trait, gene.GR$ORF)]
full2d.table$m1.cis=(f2d.tchr==m1$chr & abs(f2d.tpos-m1$start)<5e4)  
full2d.table$m2.cis=(f2d.tchr==m2$chr & abs(f2d.tpos-m2$start)<5e4)

f2dts=list();
#split(full2d.table, full2d.table$trait)
for(i in 1:nrow(full2d.table)) {
    #i=1
    x=full2d.table[i,]
    yr=pheno.scaled.OD[,x$trait]
    pmarkers=unique(c(x$m1, x$m2))
    apeaks=match(pmarkers, colnames(gdata))
    X=data.frame(gdata[,apeaks])
    lm.goi=lm(yr~.^2-1,data=X)
    aov.a = anova(lm.goi)
    tssq  = sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    csl=coefficients(summary(lm.goi))
   f2dts[[i]]=data.frame(trait=x$trait,
               m1_pmarker=x$m1,
               m1_beta=coef(lm.goi)[1],
               m1_variance_explained=a.effs[1],
               m1_p=csl[1,4],
               m1_is_cis=x$m1.cis,
               m2_pmarker=x$m2,
               m2_beta=coef(lm.goi)[2],
               m2_variance_explained=a.effs[2],
               m2_p=csl[2,4],
               m2_is_cis=x$m2.cis,
               interaction_beta=coef(lm.goi)[3],
               interaction_variance_explained=a.effs[3],
               interaction_p=csl[3,4])
}
lm2D_from_full=do.call('rbind', f2dts)
#save(lm2D_from_full,file='/data/eQTL/RData/lm2D_full.RData')

library(tidyr)
m1=data.table::tstrsplit(marginal2d.table$m1, '_|:', type.convert=T)
m2=data.table::tstrsplit(marginal2d.table$m2, '_|:', type.convert=T)
names(m1)=c('chr','start','var')
names(m2)=c('chr','start','var')
m1$end=m1$start
m2$end=m2$start

m2d.tchr=as.character(seqnames(gene.GR)[match(marginal2d.table$trait, gene.GR$ORF)])
m2d.tpos=start(gene.GR)[match(marginal2d.table$trait, gene.GR$ORF)]
marginal2d.table$m1.cis=(m2d.tchr==m1$chr & abs(m2d.tpos-m1$start)<5e4)  
marginal2d.table$m2.cis=(m2d.tchr==m2$chr & abs(m2d.tpos-m2$start)<5e4)


# add in code to reformat marginal table
m2dts=list();
#split(full2d.table, full2d.table$trait)
for(i in 1:nrow(marginal2d.table)) {
    #i=1
    print(i)
    x=marginal2d.table[i,]
    yr=pheno.scaled.OD[,x$trait]
    pmarkers=unique(c(x$m1, x$m2))
    apeaks=match(pmarkers, colnames(gdata))
    X=data.frame(gdata[,apeaks])
    lm.goi=lm(yr~.^2-1,data=X)
    aov.a = anova(lm.goi)
    tssq  = sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    csl=coefficients(summary(lm.goi))
   m2dts[[i]]=data.frame(trait=x$trait,
               m1_pmarker=x$m1,
               m1_beta=coef(lm.goi)[1],
               m1_variance_explained=a.effs[1],
               m1_p=csl[1,4],
               m1_is_cis=x$m1.cis,
               m2_pmarker=x$m2,
               m2_beta=coef(lm.goi)[2],
               m2_variance_explained=a.effs[2],
               m2_p=csl[2,4],
               m2_is_cis=x$m2.cis,
               interaction_beta=coef(lm.goi)[3],
               interaction_variance_explained=a.effs[3],
               interaction_p=csl[3,4])
}
lm2D_from_marginal=do.call('rbind', m2dts)
save(lm2D_from_marginal,file='/data/eQTL/RData/lm2D_marginal.RData')




cint=cbind(lm.2D.df.long$m2.gcoord,lm.2D.df.long$m1.gcoord)
cint=t(apply(cint,1, function(x) if(x[1]<x[2]) {return(x)} else { return(rev(x)) }))

plot(cint[,1], cint[,2], col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'),
     #lm.2D.df.long$partner.cis+1, 
     pch=21, xlab='QTL 1 position', ylab='QTL 2 position')
points(cint[,2], cint[,1],  col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'),
       #col= lm.2D.df.long$partner.cis+1, 
       pch=21)
abline(v=gcoord.key, lty=2, col='lightblue')
abline(h=gcoord.key, lty=2, col='lightblue')

pdf(file='~/Desktop/2Dplot_ideas.pdf', width=40, height=20)
par(xaxs='i', yaxs='i', mfrow=c(1,2))
plot(cint[,1], cint[,2], col= lm.2D.df.long$partner.cis+1, pch=21, xlab='QTL 1 position', ylab='QTL 2 position',
     xlim=c(0, max(cint[,2])), ylim=c(c(0, max(cint[,2]))))
points(cint[,2], cint[,1], col= lm.2D.df.long$partner.cis+1, pch=21)
abline(v=gcoord.key, lty=2, col='lightblue')
abline(h=gcoord.key, lty=2, col='lightblue')

full2d.table$m1.gcoord=marker.GR$gcoord[match(full2d.table$m1, marker.GR$mname)]
full2d.table$m2.gcoord=marker.GR$gcoord[match(full2d.table$m2, marker.GR$mname)]

marginal2d.table$m1.gcoord=marker.GR$gcoord[match(marginal2d.table$m1, marker.GR$mname)]
marginal2d.table$m2.gcoord=marker.GR$gcoord[match(marginal2d.table$m2, marker.GR$mname)]
pdf(file='~/Desktop/2Dscans_2.pdf',width=11, height=11)
#par(mfrow=c(1,3)) #, xaxs='i', yaxs='i')
plot(full2d.table$m1.gcoord, full2d.table$m2.gcoord, xlab='QTL 1 position', ylab='QTL 2 position', main='', xaxt='n', yaxt='n',
     xlim=c(0, max(cint[,2])), ylim=c(c(0, max(cint[,2]))), pch=20, col='#0000FF44', cex=2)
#abline(v=gcoord.key, lty=2, col='lightblue')
#abline(h=gcoord.key, lty=2, col='lightblue')

points(marginal2d.table$m1.gcoord, marginal2d.table$m2.gcoord, xlab='QTL 1 position', ylab='QTL 2 position', main='Marginal Scan',xaxt='n', yaxt='n',
     xlim=c(0, max(cint[,2])), ylim=c(c(0, max(cint[,2]))), pch=20, col='#FFFF0044', cex=2)
#abline(v=gcoord.key, lty=2, col='lightblue')
#abline(h=gcoord.key, lty=2, col='lightblue')

points(cint[,1], cint[,2], xlab='QTL 1 position', ylab='QTL 2 position', main='QTL Scan',xaxt='n', yaxt='n',
     xlim=c(0, max(cint[,2])), ylim=c(c(0, max(cint[,2]))), pch=20, col='#00000044', cex=2)
abline(v=gcoord.key, lty=2, col='lightblue')
abline(h=gcoord.key, lty=2, col='lightblue')
dev.off()

#plot(lm.2D.df.long$m1.gcoord, lm.2D.df.long$t.gcoord, col=lm.2D.df.long$partner.cis+1, pch=20, xlab='QTL position', ylab='transcript position')
#points(lm.2D.df.long$m2.gcoord, lm.2D.df.long$t.gcoord, col=lm.2D.df.long$partner.cis+1, pch=20)
#segments(lm.2D.df.long$m1.gcoord, lm.2D.df.long$t.gcoord, lm.2D.df.long$m2.gcoord, lm.2D.df.long$t.gcoord,
#         col=ifelse(lm.2D.df.long$partner.cis, '#00000055', '#ff000044'))

plot(cint[,1], lm.2D.df.long$t.gcoord, 
     col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'),
     #col=lm.2D.df.long$partner.cis+1, 
     pch=20, xlab='QTL position', ylab='transcript position',xlim=c(0, max(cint[,2])), ylim=c(c(0, max(cint[,2]))))
points(cint[,2], lm.2D.df.long$t.gcoord, 
       col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'),
       #col=lm.2D.df.long$partner.cis+1, 
       pch=20)
segments(cint[,1], lm.2D.df.long$t.gcoord, cint[,2], lm.2D.df.long$t.gcoord,
         col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'))
         #ifelse(lm.2D.df.long$partner.cis, '#00000055', '#ff000044'))
abline(v=gcoord.key, lty=2, col='lightblue')
abline(h=gcoord.key, lty=2, col='lightblue')
dev.off()

#plot(Intervals(cint)[order(cint[,2]-cint[,1])] )
##, col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'), use_point=F, use_names=F)
#library(circlize)
#dfc=data.frame(from=cint[,1], to=cint[,2])
#plot(Intervals(cint), col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'), use_point=F, use_names=F)
#plot(Intervals(cint), col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'),  use_point=F, use_names=F)


vc_multiple_components=list()
A=tcrossprod(gdata.scaled)/ncol(gdata.scaled)
AA=A*A
pb =txtProgressBar(min = 1, max = length(peaks.per.gene), style = 3)
for(gg in 1:length(peaks.per.gene)) {
    setTxtProgressBar(pb, gg)
    g=names(peaks.per.gene)[gg]
    yr=pheno.scaled.OD[,g]
    ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
    apeaks = match( ppg$pmarker, colnames(gdata))

    pms=ppg$pmarker
    max.peak=ppg$pmarker[which.max(ppg$LOD)]
    qA=tcrossprod(gdata.scaled[,pms])/length(pms)
   
    rr=regress(yr~1, ~A+AA, verbose=F, pos=c(T,T,T,T,T,T))
    vc_multiple_components[[g]][['A_AA']]=extract.rr(rr)

    ##rr=regress(yr~1, ~qA+Ahot+A, verbose=F, pos=c(T,T,T,T,T,T,T))
    ##vc_multiple_components[[g]][['1']]=extract.rr(rr)
    
    if(sum(ppg$cis)>0 & sum(!ppg$cis>0) ) {
        A.local=tcrossprod(gdata.scaled[,pms[ppg$cis]])/sum(ppg$cis)
        A.distant=tcrossprod(gdata.scaled[,pms[!ppg$cis]])/sum(!ppg$cis)
        A.local_X_A=A.local*A
        A.distant_X_A=A.distant*A
        rr=regress(yr~1, ~A.local + A.distant + A + A.local_X_A + A.distant_X_A + AA, verbose=F, pos=c(T,T,T,T,T,T,T))
        vc_multiple_components[[g]][['local_distant']]=extract.rr(rr)
    }    
 
    if(!is.null(lm.2D[[g]]) ) {
        g2in=model.matrix(fit2)
        g2int=scale(g2in[,grep(':', colnames(g2in))])
        qAA=tcrossprod(g2int)/ncol(g2int)

        rr=  regress(yr~1, ~qA+A+qAA+AA, verbose=F, pos=c(T,T,T,T,T))
        vc_multiple_components[[g]][['QTL']]=extract.rr(rr)
    }
}
close(pb)
#save(vc_multiple_components, file='/data/eQTL/RData/vc_multiple_components_042617.RData')
stripchart(data.frame(t(sapply(vc_multiple_components, function(x) x[['A_AA']]$sigma))), vertical=T, method='jitter', pch=21, col='#00000022', ylim=c(0,1), ylab='fraction of phenotypic variance', group.names=c('A','AA', 'E'))
boxplot(t(sapply(vc_multiple_components, function(x) x[['A_AA']]$sigma)), add=T, border='darkred', outline=FALSE, boxwex=.5, xaxt='n')

stripchart(data.frame((do.call('rbind', sapply(vc_multiple_components, function(x) x[['local_distant']]$sigma)))), vertical=T, method='jitter', pch=21, col='#00000022', ylim=c(0,1), ylab='fraction of phenotypic variance',
            group.names=c('local additive', 'distant additive', 'unmapped additive', 'local by genome, interactions', 'distant by genome, interactions', 'interactions between unmapped QTL', 'E')
           )
boxplot((do.call('rbind', sapply(vc_multiple_components, function(x) x[['local_distant']]$sigma))), add=T, border='darkred', outline=FALSE, boxwex=.5, xaxt='n')


stripchart(data.frame((do.call('rbind', sapply(vc_multiple_components, function(x) x[['QTL']]$sigma)))), vertical=T, method='jitter', pch=21, col='#00000022', ylim=c(0,1), ylab='fraction of phenotypic variance', 
           group.names=c('additive QTL', 'unmapped additive', 'interaction QTL', 'unmapped interactions', 'E'))
boxplot((do.call('rbind', sapply(vc_multiple_components, function(x) x[['QTL']]$sigma))),add=T, border='darkred', outline=FALSE, boxwex=.5 , xaxt='n')


lds=do.call('rbind', sapply(vc_multiple_components, function(x) x[['local_distant']]$sigma))
plot(vest[1,], vest[1,]+vest[2,], xlim=c(0,1), ylim=c(0,1))
plot(lds[1,], type='l', col='#00000009', ylim=c(0,1))
for(i in 1:2846) {
points(lds[i,], type='l', col='#00000009', ylim=c(0,1))
}

