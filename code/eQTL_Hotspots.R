# NEW HOTSPOT ANALYSIS  020117

#cc='chrIII'
#hotspots=list()
#total.hotspot.n=0
#unique.chrs=paste0('chr', as.roman(1:16))
# simplify further
compositeMV.LOD.list=list()
detected1D.list=list()
bootstrap.list=list()


# constitutes after1D.hotspots list()
hotspot.pheno=list()
hotspot.geno=list()
only1D=list()


#### code to detect hotspots
#for(cc in unique.chrs) { 
for(cc in unique.chrs[-c(17)] ) {
       print(cc) 
       #for each chromosome keep the composite LOD score, the matrix of detected QTL positions, and the input (phenotype) into the multivariate mapping 
       #LODthresh=2.5
       #cc='chrXV'
       #------------------------------------------------------------- 
       g.r=gdata.by.chr[[cc]]
       g.s=gdata.s.by.chr[[cc]]
       gd=duplicated(g.s, MARGIN=2)
       if(sum(gd>0) ) {   g.s=g.s[,-gd] }

       g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
       names(g.s.pos)=colnames(g.s)
       #hs.n=1
      
       ppc=peaks.per.chr[[cc]]
       keep.transcripts=colnames(t.tpm.matrix)
      
       # remove local QTL -------------------------------------------
       keep.transcripts=unique(gene.annot.df$name[gene.annot.df$chr!=cc]  )
       ppc3=ppc[ppc$gene %in% keep.transcripts,]
      
       #check that this was executed
       keep.transcripts=unique(ppc3$gene)

       # visualize
       #par(mfrow=c(2,1))
       #plot(g.s.pos[ppc3$pcind], ppc3$LOD*sign(ppc3$r), ylim=range( ppc3$LOD*sign(ppc3$r)), cex=.7, main=cc, ylab='LOD', xlab='1D peak position') #ppc2$LOD*sign(ppc2$r)
       plot(g.s.pos[ppc3$pcind], ppc3$r, ylim=range( -1,1), cex=.7,  main=cc, ylab='r', xlab='1D peak position') #ppc2$LOD*sign(ppc2$r)

       tmm=model.QTL.effects(t.tpm.matrix, gdata, keep.transcripts, peaks.per.gene, background.covariates,c(),cc)
       #tmm=model.QTL.effects(t.tpm.matrix, gdata, names(to.scan), peaks.per.gene,add.cov, pmarker,cc)
       tmms=scale(tmm)

       hotspot.pheno[[cc]]=tmms
       hotspot.geno[[cc]]=g.s
       
       # 1/21 was here ######################3
       tmms.svd=svd(tmms)
       tmms.svd.null=svd(apply(tmms, 2, function(x) x[sample(1:length(x))] ))
       #x11()
       #plot(tmms.svd$d^2/sum(tmms.svd$d^2), xlim=c(0,200))
       #points(tmms.svd.null$d^2/sum(tmms.svd.null$d^2), xlim=c(0,200), col='blue')
       pc.null=tmms.svd.null$d^2/sum(tmms.svd.null$d^2)
       pc.to.keep=which((tmms.svd$d^2/sum(tmms.svd$d^2)  > pc.null) ) #& tmms.svd$d>pc.null[1000])
       Ysub=tmms.svd$u[,pc.to.keep]

       #Ysub=icf$K
       mnull.1=lm(Ysub~1)
       #det_AtA(residuals(mnull.1))
       RSSn.1=crossprod(residuals(mnull.1))
       ldetRSSn.1=determinant(RSSn.1, logarithm=T)$modulus
       
       detected1D=matrix(1,1012,1)
       p.detected=c()
       
       all.null=c()
       repeat{
           #scanone -------------------------------------------------------------- 
           yR=residuals(.lm.fit(detected1D, Ysub))
           ldetnull=determinant(crossprod(residuals(lm(yR~1))), logarithm=T)$modulus
           ldetRSSf.1=mvn.scanone(g.s, yR)
           mvLOD.1=1012*(ldetnull-ldetRSSf.1)/(2*log(10))
           #plot(mvLOD.1)
           print(paste('maxLOD', max(mvLOD.1),which.max(mvLOD.1) ))
           #--------------------------------------------------------------------
              
           #scanone null----------------------------------------------------------
           n.perm.m=100
           mvLOD.1null=rep(0,n.perm.m)
           pb =txtProgressBar(min = 1, max = n.perm.m, style = 3)
           for(k in 1:n.perm.m) {
               setTxtProgressBar(pb, k)
               ldetperm=mvn.scanone(g.s, yR[sample(1:nrow(yR)),]) #tmms[sample(1:nrow(tmms)),to.test.vec])
               #mvLOD.1null[k]=max(1012*(ldetRSSn.1-ldetperm)/(2*log(10)))
               #fL=fasterLOD(1012, scale(Ysub), g.s)
               mvLOD.1null[k]=max(1012*(ldetnull-ldetperm)/(2*log(10)))
           }
           close(pb)
           all.null=c(all.null, mvLOD.1null)

          #--------------------------------------------------------------------
            if(max(mvLOD.1)< (quantile(mvLOD.1null,.99) ) ) { 
                break 
            } else {
             print(paste('perm thresh', max(mvLOD.1null) ) ) # quantile(mvLOD.1null,.99)))
             p.detected=c(p.detected, which.max(mvLOD.1))
             print(p.detected)
             detected1D=cbind(detected1D, g.s[,which.max(mvLOD.1)])
            }
       }
       colnames(detected1D)=c('Intercept', colnames(g.s)[p.detected])


         #composite multivariate scan
       # ****
       compositeMV.LOD=calc.composite.MV.LOD(detected1D, Ysub, g.s) 
       # plot it  
       #plotMV.LOD(compositeMV.LOD, g.s)
       #relocate

       # error here chrVII collapsing 1D peaks !! (debugging)
       #detected1D.2=relocate.peaks(compositeMV.LOD, g.s, g.s.pos, detected1D)
       # double check that the same amount of remapping is done for chrI and chrIII??

       detected1D=relocate.peaks(compositeMV.LOD, g.s, g.s.pos, detected1D)
       
       #recompute composite LOD 
       compositeMV.LOD=calc.composite.MV.LOD(detected1D, Ysub, g.s) 
       #x11()
       plotMV.LOD(compositeMV.LOD, g.s, g.s.pos)

       #x11()
       #plotMV.LOD(compositeMV.LOD2, g.s)
       
       only1D[[cc]]$detected1D=detected1D
       only1D[[cc]]$compositeMV.LOD=compositeMV.LOD
}
after1D.hotspots=list(hotspot.pheno=hotspot.pheno, hotspot.geno=hotspot.geno, only1D=only1D)

# This is a very important data structure !!!
#save(after1D.hotspots, file= '/data/eQTL/RData/after1D.hotspots.RData') 
#file/home/jbloom/Dropbox/Public/eQTL/stranded/

load('/data/eQTL/RData/after1D.hotspots.RData')
only1D=after1D.hotspots$only1D
hotspot.pheno=after1D.hotspots$hotspot.pheno
hotspot.geno=after1D.hotspots$hotspot.geno

# now scan for ghost QTL 
compositeMV.LOD.list=list()
detected1D.list=list()
for (cc in unique.chrs[-17] ) {
    #cc='chrVII'
    detected1D= only1D[[cc]]$detected1D
    compositeMV.LOD= only1D[[cc]]$compositeMV.LOD
    tmms=hotspot.pheno[[cc]]
    
    tmms.svd=svd(tmms)
    tmms.svd.null=svd(apply(tmms, 2, function(x) x[sample(1:length(x))] ))
    pc.null=tmms.svd.null$d^2/sum(tmms.svd.null$d^2)
    pc.to.keep=which((tmms.svd$d^2/sum(tmms.svd$d^2)  > pc.null) ) #& tmms.svd$d>pc.null[1000])
    Ysub=tmms.svd$u[,pc.to.keep]

    g.s=hotspot.geno[[cc]]
    g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
    names(g.s.pos)=colnames(g.s)

    plotMV.LOD(compositeMV.LOD, g.s, g.s.pos)

    record1D=detected1D
    recordCM=compositeMV.LOD


    peakMVN=apply(compositeMV.LOD, 1, max)
    max.peak=max(peakMVN)
    wm.mp=which.max(peakMVN)
    wmax.peak=names(which.max(peakMVN))

   # potential threshold for 2D
   #max(all.null)*3
   #peak.to.fine.map=names(which(peakMVN>200))
   #mvLOD.thresh=max(all.null)*4
   
   mvLOD.thresh=200

   found2D=FALSE
   peak.blacklist=c()
   print('scanning for 2D peaks')
   while( max.peak> mvLOD.thresh ) {
            
            pname=wmax.peak
            print(pname)
            #roi=range(which(mvLOD.1>(max(mvLOD.1)-max(mvLOD.1)/2)))
            roi1=getROI(compositeMV.LOD[pname,],g.s)
            pm=match(pname, colnames(detected1D))
            pii=match(wmax.peak, colnames(g.s))
            peak.blacklist=c(peak.blacklist, c((pii-10):(pii+10)))

            #narrow scantwo
            ldetRSSn.2=determinant(crossprod(residuals(lm(Ysub~detected1D[,-pm]))), logarithm=T)$modulus
            #ldetRSSn.2=determinant(crossprod(residuals(lm(Ysub~detected1D))), logarithm=T)$modulus
            ldetRSSf.2=mvn.scantwo(g.s,  Ysub, roi1, add.cov=detected1D[,-pm]) #, g.s[,which.max(mvLOD.1.1)])
            
            Y.reduced=residuals(lm(Ysub~detected1D[,-pm]))
            #ldetRSSn.2=determinant(crossprod(residuals(lm(Y.reduced~1))), logarithm=T)$modulus
            #ldetRSSf.2=mvn.scantwo(g.s,  Y.reduced, roi1) #, g.s[,which.max(mvLOD.1.1)])

            mvLOD.2=1012*(-ldetRSSf.2+ldetRSSn.2)/(2*log(10))
            
            mv2Dpeak=arrayInd(which.max(mvLOD.2), dim(mvLOD.2))+(min(roi1)-1)
            sig2D=max(mvLOD.2, na.rm=T)-max(compositeMV.LOD[pname,]) 
            
            #scantwo null
            n.perm.m=100
            mv.scantwo.null=rep(0,n.perm.m)
            pb =txtProgressBar(min = 1, max = n.perm.m, style = 3)
            for(k in 1:n.perm.m) {
                 Y=Y.reduced[sample(1:nrow(Y.reduced)),]
                 #ldetRSSn.2null=determinant(crossprod(residuals(lm(Y~1))), logarithm=T)$modulus

                 ldetperm=mvn.scanone(g.s, Y)
                 #mvLOD.1null.all=(1012*(ldetRSSn.1-ldetperm)/(2*log(10)))
                 mvLOD.1null.all=(1012*(ldetRSSn.2-ldetperm)/(2*log(10)))
                 mvLOD.1null=max(mvLOD.1null.all)
                 roinull=getROI(mvLOD.1null.all,g.s)
                 
                 ldetRSSf.2null=mvn.scantwo(g.s, Y, roinull)
                 mvLOD.2null=1012*(-ldetRSSf.2null+ldetRSSn.2)/(2*log(10))
                 mv.scantwo.null[k]=max(mvLOD.2null, na.rm=T)-max(mvLOD.1null, na.rm=T)
                 setTxtProgressBar(pb, k)
            }
            close(pb)
            print(paste(sig2D, (quantile(mv.scantwo.null,.99)) ))

            sig.suff1= sig2D> (quantile(mv.scantwo.null,.99))+10
            sig.suff2= sum((c((pii-3):(pii+3)) %in% mv2Dpeak))==0

            if(sig.suff1 & sig.suff2) { 
                found2D=TRUE
                print( colnames(g.s)[mv2Dpeak] )
                detected1D=detected1D[,-pm]
                detected1D=cbind(detected1D, g.s[, mv2Dpeak[1]], g.s[, mv2Dpeak[2]])
                colnames(detected1D)[ncol(detected1D)-1]=colnames(g.s)[mv2Dpeak[1]]
                colnames(detected1D)[ncol(detected1D)]  =colnames(g.s)[mv2Dpeak[2]]
             }
                peakMVN=apply(compositeMV.LOD, 1, max)
                pii=match(names(peakMVN), colnames(g.s))
                peakMVN=peakMVN[!(pii %in% peak.blacklist)]
                max.peak=max(peakMVN)
                wm.mp=which.max(peakMVN)
                wmax.peak=names(which.max(peakMVN))
      }
      if(found2D) {
         
          dstore=detected1D
          cstore=compositeMV.LOD

          # recalculate composite LOD
          compositeMV.LOD=calc.composite.MV.LOD(detected1D, Ysub, g.s)
          #plotMV.LOD(compositeMV.LOD, g.s)
          # relocate peaks 
          detected1D=relocate.peaks(compositeMV.LOD, g.s, g.s.pos, detected1D)
          
          compositeMV.LOD=calc.composite.MV.LOD(detected1D, Ysub, g.s)
          plotMV.LOD(compositeMV.LOD, g.s, g.s.pos)
      }
    compositeMV.LOD.list[[cc]]=compositeMV.LOD
    #=list()
    detected1D.list[[cc]]=detected1D
}
#save(compositeMV.LOD.list,file='/data/eQTL/RData/compositeMV.LOD.list.RData')
#save(detected1D.list, file='/data/eQTL/RData/detected1D.list.RData') 


load('/data/eQTL/RData/compositeMV.LOD.list.RData')
load('/data/eQTL/RData/detected1D.list.RData') 


bootstrap.list=list()
boot.window=list()
#[[cc]]boots
#pdf('/data/eQTL/RData/mvLODplots.pdf', width=14, height=8)
for (cc in unique.chrs[-17] ) {
    #cc='chrVII'
    detected1D=detected1D.list[[cc]]
    compositeMV.LOD=compositeMV.LOD.list[[cc]]
    
    tmms=hotspot.pheno[[cc]]
    g.s=hotspot.geno[[cc]]
    g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
    names(g.s.pos)=colnames(g.s)

    tmms.svd=svd(tmms)
    tmms.svd.null=svd(apply(tmms, 2, function(x) x[sample(1:length(x))] ))
    pc.null=tmms.svd.null$d^2/sum(tmms.svd.null$d^2)
    pc.to.keep=which((tmms.svd$d^2/sum(tmms.svd$d^2)  > pc.null) ) #& tmms.svd$d>pc.null[1000])
    Ysub=tmms.svd$u[,pc.to.keep]

    g.2=gdata.by.chr[[cc]]
    g.s.pos2=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.2),':'), function(x)x[2])), '_'), function(x) x[1]))

    # 020317 ... extract exact window tested
    g.r=gdata.by.chr[[cc]]
    g.s=gdata.s.by.chr[[cc]]
    gd=duplicated(g.s, MARGIN=2)
    if(sum(gd>0) ) {   g.s=g.s[,-gd] }
    g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
    names(g.s.pos)=colnames(g.s)
    dfin=data.frame(detected1D[,-1])
   
    bet = calcIndividualEffects_afterMV( tmms, detected1D, dfin)    
      
      boots=list()
      for(p in 1:ncol(dfin)) {
        #p=10 
        peaks=colnames(detected1D)[p+1] 
        pm = match(peaks, colnames(detected1D))
        boot.window[[cc]][[peaks]]=range(g.s.pos[getROI(compositeMV.LOD[p,],g.s)])

        betr=names(which(sapply(bet, function(x) sum(names(x$coefficients) %in% colnames(dfin)[p]))>0 ))
        print(paste('calculate CI', peaks, length(betr)))
      
        svoi=svd(tmms[,betr])
        svor=svd(apply(tmms[,betr],2, function(x) sample(x)))
        #plot(svoi$d^2/sum(svoi$d^2), xlim=c(0,200))
        #points(svor$d^2/sum(svor$d^2), xlim=c(0,200), col='blue')
        pc.to.keep2=which(svoi$d^2/sum(svoi$d^2)  > svor$d^2/sum(svor$d^2) )
        Yn=svoi$u[,pc.to.keep2]
        pN=determinant(crossprod(residuals(lm(Yn~detected1D[,-pm]))), logarithm=T)$modulus
        pM=mvn.scanone(g.s, Yn, add.cov=detected1D[,pm])
        mvLOD.1=1012*(pN-pM)/(2*log(10))
      
        pm = match(peaks, colnames(detected1D))
        roi1=getROI(compositeMV.LOD[p,],g.s)
        
        bootpos=rep(NA,1000)
        pb =txtProgressBar(min = 1, max = 1000, style = 3)
          for(bb in 1:1000) {
             smpme=sample(1:1012, replace=T)
             #bnull=determinant(crossprod(residuals(lm(Yn[smpme,]~detected1D[smpme,-pm]))), logarithm=T)$modulus
             bnull=determinant(crossprod(residuals(lm(Yn[smpme,]~1))), logarithm=T)$modulus
             pM=mvn.scanone(g.s[smpme,], Yn[smpme,],add.cov=detected1D[smpme,-pm] , roi1 )
             bmvLOD.1=1012*(bnull-pM)/(2*log(10))
             #bmvLOD.2=1012*(bnull-pM2)/(2*log(10))
             #bmv2Dpeak=arrayInd(which.max(bmvLOD.2), dim(bmvLOD.2))+(min(roi1)-1)
             #print(bmv2Dpeak)
             setTxtProgressBar(pb, bb)
             #print(which.max(bmvLOD.1))
             bootpos[bb]=which.max(bmvLOD.1)
        }
        close(pb)
        print(colnames(g.s)[quantile(bootpos, c(.025, .975))])
        boots[[peaks]]=bootpos
      }
       bootstrap.list[[cc]]=boots
}

#save(bootstrap.list, file='/data/eQTL/RData/bootstrap.list.RData')

mv.bootstrap.stats=list()
mv.bootstrap.boots=list()
for(cc in unique.chrs[-17]) {
   y=bootstrap.list[[cc]]
   hg=hotspot.geno[[cc]]
   hgp=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(hg),':'), function(x)x[2])), '_'), function(x) x[1]))
   index.range=t(sapply(y, function(x) range(x)))

   index.range=cbind(index.range, (index.range[,2]-index.range[,1])<40)

   bCI=(t(sapply(y, function(x) 
                  hgp[quantile(na.omit(x), c(.025, .05, .125, .25,.5,.75, .875, .95, .975))] ) ))

   allB=t(sapply(y, function(x) hgp[sort(x)]))
   
   bCI=cbind(bCI, index.range[,3])
   colnames(bCI)=c(.025, .05, .125, .25,.5,.75, .875, .95, .975, 'localization.confidence')
   bSort=order(bCI[,4])
   bCI=bCI[bSort,] 
   allB=allB[bSort,]

   mv.bootstrap.stats[[cc]]=bCI
   mv.bootstrap.boots[[cc]]=allB
}
mv.bootstrap.table=do.call('rbind', mv.bootstrap.stats)
mv.boots=do.call('rbind', mv.bootstrap.boots)
#save(mv.boots, file='/data/eQTL/RData/mv.bootstrap.pos.RData')
#WriteXLS(data.frame(mv.bootstrap.table),'/data/eQTL/RData/MV.12417.xlsx', row.names=T)
#save(mv.bootstrap.table, file='/data/eQTL/RData/mv.bootstrap.table.RData')

#mv.bootstrap.table[,9]-mv.bootstrap.table[,8]>39
#mv.bootstrap.table[,9]-mv.bootstrap.table[,8]>39
#sort(match(rownames(mv.bootstrap.table), colnames(gdata.scaled)))








