load('/data/eQTL/RData/detected1D.list.RData')
load('/data/eQTL/RData/compositeMV.LOD.list.RData')

hotspot.index=sort(match(as.vector(unlist(sapply(compositeMV.LOD.list, rownames))), colnames(gdata.scaled)))
hotspot.index.chr=sapply(sapply(compositeMV.LOD.list, rownames), function(x) sort(match(x, colnames(gdata.scaled))))
                                       
Xreduced=gdata.scaled[,hotspot.index]
Xr2=(gdata[,hotspot.index]+1)/2
g.s=gdata.scaled

peaks.per.chr.a=split(all.peaks.OD, all.peaks.OD$chr)

best.model.markers=list()

# get the best cis marker for each gene
# either an observed local linkage peak or the closest marker to each transcript
#CIS.marker.matrix=matrix(0, 1012, length(peaks.per.gene))
CIS.marker.chosen=rep(NA, length(peaks.per.gene))
for(ii in 1:length(peaks.per.gene) ) {
      i=names(peaks.per.gene)[ii]
      pgi=peaks.per.gene[[i]]
      cis.marker=pgi$pmarker[pgi$cis][1]
      if(is.na(cis.marker) ) { cis.marker=colnames(gdata)[closest.marker.to.transcript[i]]  }
      CIS.marker.chosen[ii]=cis.marker
}
names(CIS.marker.chosen)=names(peaks.per.gene)
#CIS.marker.chosen=CIS.marker.chosen[!is.na(CIS.marker.chosen)]
#g.CIS=gdata.scaled[,CIS.marker.chosen]


best.model.markers=list()
for(cc in unique.chrs[-c(17)] ) {
      print(cc) 
      ppc=peaks.per.chr.a[[cc]]
      keep.transcripts=unique(ppc$gene)
      keep.transcripts=keep.transcripts[keep.transcripts %in% names(CIS.marker.chosen)]
      print('calculating gene expression residual')
      tmm=model.QTL.effects(t.tpm.matrix, gdata, keep.transcripts, peaks.per.gene, background.covariates,c(),cc)
      Xrset=gdata[,hotspot.index.chr[[cc]]]
      
      n.perm=10
      sampleme=replicate(n.perm, sample(1:1012))
      tmm.s=scale(tmm)
      Xrset=scale(Xrset)
  
      term.ps=list()
      term.ps.null=list()
      print('doing forward scan with permutations')
      pb =txtProgressBar(min = 1, max =length(keep.transcripts), style = 3)
      for(ii in 1:length(keep.transcripts) ) {
            setTxtProgressBar(pb, ii)
            i=keep.transcripts[ii]
            cis.marker=CIS.marker.chosen[i]
            if( gene.annot.df$chr[gene.annot.df$name==i]==cc ) { 
                     if(cis.marker %in% colnames(Xrset)) { Xr=Xrset }
                     else {        
                         Xr=cbind(gdata[,cis.marker],Xrset)  
                         colnames(Xr)[1]=cis.marker
                    }
            }
            else { Xr=Xrset }
            l1=regsubsets(Xr,tmm.s[,i], method='forward', intercept=FALSE, nbest=1)
            term.ps[[i]]=l1
            for(k in 1:n.perm) {
                l2=regsubsets(Xr,tmm.s[sampleme[,k],i], method='forward', intercept=FALSE, nbest=1)
                term.ps.null[[i]][[k]]=l2
            }
      }
     close(pb)
     stest= seq(0.00001, 0.025, .00001)

     mstat=rep(0, length(keep.transcripts))
     names(mstat)=keep.transcripts
     mstatp=matrix(0, length(keep.transcripts), n.perm)
     rownames(mstatp)=keep.transcripts
         
     nstats=ncol(Xrset)+1
     thresh=rep(NA, nstats)
     best.model=rep(0, length(keep.transcripts))
     names(best.model)=keep.transcripts
     for(i in 1:nstats ) {
         maxr2.1=sapply(term.ps, function(x) summary(x)$rsq[i] )
         maxr2.perm=t(sapply(term.ps.null, function(y) sapply(y, function(x) summary(x)$rsq[i] ) ))
         thresh[i]=doFDR(stest, maxr2.1-mstat, maxr2.perm-mstatp, .05)
         sig.transcripts=(maxr2.1-mstat>thresh[i])
         best.model[which(sig.transcripts)]=i
         cnt.sig=sum(maxr2.1-mstat>thresh[i], na.rm=T)
         print(cnt.sig)
         if(cnt.sig==0) {break;}
         mstat=maxr2.1
         mstatp=maxr2.perm
         #if(sum(maxr2.1-mstat>thresh[i])==0) {break;}
     }

     for(i in keep.transcripts) {
        best.model.markers[[i]]=c(best.model.markers[[i]], names(which(summary(term.ps[[i]])$which[best.model[i],])))
     }
}
#save(best.model.markers, file='/data/eQTL/RData/bestModelMarkers.Forward.RData')
load('/data/eQTL/RData/bestModelMarkers.Forward.RData')

B=matrix(0, length(hotspot.index), length(best.model.markers))
rownames(B)=colnames(gdata)[hotspot.index]
colnames(B)=names(best.model.markers)

B.cis=matrix(0, ncol(gdata), length(best.model.markers))
rownames(B.cis)=colnames(gdata)
colnames(B.cis)=names(best.model.markers)

#ii=names(best.model.markers)[1]
pb =txtProgressBar(min = 1, max =length(best.model.markers), style = 3)
Yr=matrix(NA, 1012, length(best.model.markers))
colnames(Yr)=names(best.model.markers)

for(i in 1:length(best.model.markers) ) { 
    ii=names(best.model.markers)[i]
    setTxtProgressBar(pb, i)
    
    nZm=best.model.markers[[ii]]
    gdata.sub=gdata[,nZm]
    
    dall=cbind(background.covariates, gdata.sub)
    model0=lm(t.tpm.matrix[,ii]~.-1, data=data.frame(background.covariates)) 
    model1=lm(t.tpm.matrix[,ii]~.-1, data=data.frame(dall)) 
    model2=lm(scale(residuals(model0))~.-1, data=data.frame(gdata.sub)) 
    
    Yr[,i]=residuals(model2)
        
    coefm2=coef(model2)
    names(coefm2)=nZm

    hpm=nZm %in% rownames(B)
    hpm.cis=!(nZm %in% rownames(B))

    B[nZm[hpm],ii]=coefm2[hpm]
    B.cis[nZm[hpm.cis],ii]=as.vector(coefm2[hpm.cis])
}
close(pb)

#save(B, file ='/data/eQTL/RData/B.Forward.RData')
#save(B.cis, file ='/data/eQTL/RData/B.cis.Forward.RData')
load('/data/eQTL/RData/B.Forward.RData')
load('/data/eQTL/RData/B.cis.Forward.RData')

# some visualization ideas
mgc=marker.GR$gcoord[match(rownames(B), marker.GR$mname)]
pB=apply(B, 1, function(x) sum(x>0))
nB=apply(B, 1, function(x) sum(x<0))
plot(mgc,pB, type='h', col='orange', ylim=c(-3000,3500), lwd=2)
points(mgc,-nB, type='h', col='purple', lwd=2) 
abline(v=sapply(split(marker.GR$gcoord, seqnames(marker.GR)), max), lty=2, col='lightblue')
abline(v=gcoord.key)


Blist=apply(B, 1, function(x) x[x!=0])
Bcislist=apply(B.cis, 1, function(x) x[x!=0])


#pdf(file='~/Desktop/hotspot_plot_ideas.pdf', height=1920, width=1080)
mgc=marker.GR$gcoord[match(rownames(B), marker.GR$mname)]
par(mfrow=c(3,1))
plot(mgc,pB, type='h', col='orange', ylim=c(-3000,3500), lwd=2,xlim=c(1e4,12e6),xaxt='n', xlab='genome position', ylab='genes influenced')
points(mgc,-nB, type='h', col='purple', lwd=2) 
abline(v=sapply(split(marker.GR$gcoord, seqnames(marker.GR)), max), lty=2, col='lightblue')
stripchart(Blist,vertical=T, method='jitter', jitter=10000, at=mgc, xlim=c(1e4,12e6), las=2, pch=20, cex=.5, main='hotspot linkages', ylab='standardized effect size') 
abline(v=sapply(split(marker.GR$gcoord, seqnames(marker.GR)), max), lty=2, col='lightblue')
stripchart(Bcislist,vertical=T, at=marker.GR$gcoord, xlim=c(1e4,12e6), las=2, pch=20, cex=1, xaxt='n', main='local linkages', ylab='standardized effect size')
abline(v=sapply(split(marker.GR$gcoord, seqnames(marker.GR)), max), lty=2, col='lightblue')
#dev.off()
#plot(sapply(Blist, function(x) quantile(x^2, .999)), sapply(Blist, length))

