# code to detect QTL-QTL interactions 

# load expression data 
#load('/data/eQTL/RData/log2_t.tpm.matrix.RData')
# load expression covariates 
#load('/data/eQTL/RData/covariates.OD.RData')

# load marker and transcript annotation data
#load('/data/eQTL/RData/markerAnnotation.RData')
#load('/data/eQTL/RData/geneAnnotation.RData')


# do 2 locus scan --------------begin by pruning marker data by LD  ---------------------------------------------------------------------------------------------------------
#load('/data/eQTL/RData/gdata.downsampled.RData')
#load('/data/eQTL/RData/pheno.additive.removed.RData')

pheno.additive.removed.scaled=scale(pheno.additive.removed)
# sanity check that additive signal is gone
#scanoneLODS.ar=fasterLOD(nrow(pheno.additive.removed.scaled),pheno.additive.removed.scaled,gdata.scaled, betas=TRUE)

#AA=(tcrossprod(gdata.scaled)/ncol(gdata.scaled))^2
## this should be about the same as the 
#vcAA=calcA(pheno.additive.removed.scaled, AA)
#h2AA=vcAA[,1]/(vcAA[,1]+vcAA[,2])

#set.seed(10)
#yAA.perm=replicate(10, {pheno.additive.removed.scaled[sample(1:1012),]})
#WG.2D.VC.perm=sapply(1:10, function(p) {
#    vcAA=calcA(yAA.perm[,,p], AA)
#    return(vcAA[,1]/(vcAA[,1]+vcAA[,2]))
#})

#o2dvc=sapply(seq(0,1,.005), function(x) sum(h2AA > x ) )
#names(o2dvc)=seq(0,1,.005)
#e2dvc=rowMeans(apply(WG.2D.VC.perm, 2, function(y) { sapply(seq(0,1,.005) , function(x) sum( y > x ) )     }))
#names(e2dvc)= seq(0,1,.005)

#e2dvc/o2dvc
#colnames(t.tpm.matrix)[which(h2AA>.4)]
#FDR 10%
# h2AA>0.39 

load('/data/eQTL/RData/all_peaks_OD.RData')
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

# Structure 2D interactions 

lm.2D.df=list()
for(goi in names(lm.2D)) {
   print(goi) 
    #goi='YAL039C'
    #lm.2D[[goi]]
    ppga=peaks.per.gene[[goi]]
    c2d=coef(lm.2D[[goi]])
    c2d=c2d[grep(':', names(c2d))]
    c2d.df=do.call('rbind', strsplit(names(c2d), ':'))
    pmarker.mod=gsub('\\:|\\/', '.', ppga$pmarker)

    m1.ind=match(c2d.df[,1], pmarker.mod)
    m2.ind=match(c2d.df[,2], pmarker.mod)

    lm.2D.df[[goi]]=data.frame(
               coef2D=as.numeric(c2d),
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






























i


                g2in=model.matrix(fit2)
                g2int=scale(g2in[,grep(':', colnames(g2in))])
                qAA=tcrossprod(g2int)/ncol(g2int)
                rr=  regress(yr~1, ~qA+A+qAA+AA, verbose=T, pos=c(T,T,T,T,T))
                vc_multiple_components[[g]][['4']]=extract.rr(rr)
                ints=rownames(a2)[-1]
                lm.3D.test[[g]]=add1(fit2,  paste0(gsub(':|/', '.', max.peak), ':', ints), test='F')
            }
        }







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
   }
close(pb)















































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


load('/data/eQTL/RData/peaks_per_gene_augmented.RData')
load('/data/eQTL/RData/lm.2D.RData')
#peaks.per.gene.augmented
lm.2D.df=list()
for(goi in names(lm.2D)) {
   print(goi) 
    #goi='YAL039C'
    #lm.2D[[goi]]
    ppga=peaks.per.gene.augmented[[goi]]
    c2d=coef(lm.2D[[goi]])
    c2d=c2d[grep(':', names(c2d))]
    c2d.df=do.call('rbind', strsplit(names(c2d), ':'))
    pmarker.mod=gsub('\\:|\\/', '.', ppga$pmarker)

    m1.ind=match(c2d.df[,1], pmarker.mod)
    m2.ind=match(c2d.df[,2], pmarker.mod)

    lm.2D.df[[goi]]=data.frame(
               coef2D=as.numeric(c2d),
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

plot(Intervals(cint)[order(cint[,2]-cint[,1])] )
#, col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'), use_point=F, use_names=F)
library(circlize)
dfc=data.frame(from=cint[,1], to=cint[,2])


plot(Intervals(cint), col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'), use_point=F, use_names=F)
plot(Intervals(cint), col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'), 
     use_point=F, use_names=F)

plot3d(lm.2D.df.long$m1.gcoord, lm.2D.df.long$m2.gcoord, lm.2D.df.long$t.gcoord,col=lm.2D.df.long$partner.cis+1)



segments(lm.2D.df.long$m1.gcoord, lm.2D.df.long$partner.cis+1
         
         lm.2D.df.long$t.gcoord, lm.2D.df.long$m2.gcoord, lm.2D.df.long$t.gcoord,
         col=ifelse(lm.2D.df.long$partner.cis, '#00000044', '#ff000044'))


abline(v=gcoord.key)
18.7% of sig 2D coming from one partner local


# Sanity check the mapping 
f=split(full2d.table,full2d.table$trait)
f=split(marginal2d.table,marginal2d.table$trait)

library(car)
ofl=colSums(count.matrix)
for(i in 1:length(f)) {
    gene=names(f)[i]
    print(gene)
    X=data.frame(gdata[,peaks.per.gene[[gene]]$pmarker])
    fitme=glm(t(count.matrix)[,gene]~-1+offset(log(ofl))+covariates.OD+., data=X,family=poisson(link=log))

    #fitme=lm(t.tpm.matrix[,gene]~covariates.OD+., data=X) #,family=poisson(link=log))
    
    f2=f[[gene]]
    print(f2)
    niseq=1:nrow(f2)
    interms=paste(sapply(niseq, function(n) {
    paste0(c('gdata[,f2$m1[', n, ']]*gdata[,f2$m2[',n,']]'), collapse='')
          }), collapse='+')
    uf2=update(fitme, paste0('~.+', interms))
    print(   summary( uf2) ) 

    yr=residuals(glm(t(count.matrix)[,gene]~-1+offset(log(ofl))+covariates.OD,family=poisson(link=log)))
    yr2=residuals(lm(t.tpm.matrix[,gene]~-1+covariates.OD))
    #Anova(uf2, type='III')
    par(mfrow=c(1, length(niseq)))
    boxplot(yr~gdata[,f2$m1[1]]*gdata[,f2$m2[1]])
    tryCatch({
    boxplot(yr~gdata[,f2$m1[2]]*gdata[,f2$m2[2]]) }, error=function(e) {} )
    tryCatch({
    boxplot(yr~gdata[,f2$m1[3]]*gdata[,f2$m2[3]]) }, error=function(e) {} )
    readline()
}






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


obs2DVCchr=readRDS('/data/eQTL/RData/VC_2D.RDS.1')
ochr=(sapply(obs2DVCchr, function(x) x[,1]/(x[,1]+x[,2]) ))
ochr.vec=as.vector(ochr)

echr=sapply(2:6, function(y) {
            indata=readRDS(paste0('/data/eQTL/RData/VC_2D.RDS.', y))
            return(sapply(indata, function(x) x[,1]/(x[,1]+x[,2]) ))
          })

chr.o2dvc=sapply(seq(0,1,.001), function(x) sum(ochr.vec > x ) )
names(chr.o2dvc)=seq(0,1,.001)
chr.e2dvc=rowMeans(apply(echr, 2, function(y) { sapply(seq(0,1,.001) , function(x) sum( y > x ) )     }))
names(chr.e2dvc)= seq(0,1,.001)
chr.e2dvc/chr.o2dvc
# FDR 10% 0.105
# n= 63 (60 unique traits)




v2=sapply(VC_2D, function(x) sum(x[,1]>.03))
v3=do.call('rbind', strsplit(names(v2), '_'))
v4=data.frame(v3,v2, stringsAsFactors=F)
plot(v4[,1], v4[,2], type='n', ylim=c(1,16), xlim=c(1,16))
text(as.numeric(v4[,1]), as.numeric(v4[,2]),v4[,3], cex=log10(v4[,3]/10))
#------------------------------------------------------------------------------------------------------------------------

