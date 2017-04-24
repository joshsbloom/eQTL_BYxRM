# consider this in addition to or in place of cufflinks
# an alternate splicing aware aligner
# http://bioconductor.org/packages/release/bioc/html/Rsubread.html
# glms from scratch
# http://www.r-bloggers.com/poisson-regression-fitted-by-glm-maximum-likelihood-and-mcmc/


calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%crossprod(Z,Vinv)%*%(y- X%*%B)    }
#using grofit package
getGrowthStats=function(base.dir, batches=paste0('eQTL', 1:13)) {
    #batches = paste0('eQTL', 1:13)
    growth.fit=list()
    name.list=list()
    for( batch in batches) {
        print(batch)
        #batch = batches[1]
        sample.sheet.file = paste0(base.dir, 'sample_sheets/sample_sheet_', batch, '.csv')
        ss=read.delim(sample.sheet.file, header=T, sep=',', stringsAsFactors=F) #skip=20 # if miseq sample sheet
        ss_names=matrix(do.call('rbind', strsplit(ss$Sample_Name, '-'))[,1],8,12)

        growth.file = paste(base.dir, 'eQTL_growth/', batch, '.xlsx', sep='')
        scnt=sheetCount(growth.file)
        gtime=sheetNames(growth.file)
        growth=list()
        for (s in 1:scnt) {  print(s);  growth[[s]]=read.xls(growth.file, sheet=s)[,2:13]-0.079 }
        growth=abind(growth, along=3)
            for(row in 1:8) {
                for(col in 1:12) {
                    ss_name=ss_names[row,col]
                    esn=exists(ss_name, where=name.list)
                    if(esn) { ss_name=paste0(ss_name, '_', length(grep(ss_name, names(name.list)))) }
                    name.list[[ss_name]]=ss_name
                    if(dim(growth)[3]<5) {
                        g=c(min(growth[row, col, ]),  growth[row, col, ])
                        t=c(min(gtime), gtime)
                        growth.fit[[batch]][[ss_name]]=getGrowthCurve(t, g, plotme=F)
                    }
                        else{
                    growth.fit[[batch]][[ss_name]]=getGrowthCurve(growth[row, col, ], gtime, plotme=F)
                    }
                  growth.fit[[batch]][[ss_name]]$maxOD=max(growth[row, col, ])
                  growth.fit[[batch]][[ss_name]]$minOD=min(growth[row, col, ])
                #readline()
                }
            }
    }
    return(growth.fit)
}


# Two functions to calculate linkage statistics -----------------------------------------------------------------------
#get.LOD.by.COR = function(n.pheno, pheno, gdata, doGPU=F) {
#    if(doGPU) {
#    # Lynch and Walsh p. 454
#    return( (-n.pheno*log(1-gpuCor(pheno, gdata, use='pairwise.complete.obs')$coefficients^2))/(2*log(10)) )  }
#    else {
#    return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10)) )  }
#}
# this one can't handle NAs but is ridiculously fast
# input is scaled phenotypes and genotypes
fasterLOD=function(n.pheno, pheno.s,gdata.s, betas=FALSE, sdx=1, pheno=NULL){
   r=crossprod(pheno.s, gdata.s)/(n.pheno-1)
   LOD=(-n.pheno*log(1-r^2))/(2*log(10))
   if(betas==FALSE) {
       return(LOD)
   } else {
       #beta=r*apply(cbind(pheno),2, sd,na.rm=T)/sdx
       return(list(r=r, LOD=LOD))
   }
}
#----------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# functions to convert counts to other units
# from https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#countToTpm <- function(counts, effLen){
    #effLen=ifelse(effLen>.0001, effLen, 1)
    #counts=ifelse(counts>.0001, counts, 1)
#    rate <- log(counts) - log(effLen)
#    rate[!is.finite(rate)]=0
#    denom <- log(sum(exp(rate)))
#    exp(rate - denom + log(1e6)) }

countToTpm <- function(counts, effLen){
    #effLen=ifelse(effLen>.0001, effLen, 1)
    #counts=ifelse(counts>.0001, counts, 1)
    rate <- log(counts) - log(effLen)
    infFlag = !is.finite(rate)
    rate[!is.finite(rate)]=NA
    denom <- log(sum(exp(rate), na.rm=TRUE))
    ret = exp(rate - denom + log(1e6))
    ret[infFlag] <- 0
    ret
}

countToFpkm <- function(counts, effLen){
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) ) }

fpkmToTpm <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6)) }

countToEffCounts <- function(counts, len, effLen){
    counts * (len / effLen) }

# borrowed and modified from  from mixed.solve() in rrBLUP 08/02/15 ------------------------------------
# super optimized for one VC and fixed effects ~1000X speedup by precomputing eigen decomp
m.S=function (y, K = NULL, bounds = c(1e-09, 1e+09), theta=NULL, Q=NULL, X=NULL ) 
{
    n <- length(y)
    y <- matrix(y, n, 1)
    if(is.null(X) ) {  p <- 1    } else { p = ncol(X) }
    Z <- diag(n)
    m <- ncol(Z)
       
    omega <- crossprod(Q, y)
    omega.sq <- omega^2
    
    f.REML <- function(lambda, n.p, theta, omega.sq) {
        n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
    }
    soln <- optimize(f.REML, interval = bounds, n - p, theta,  omega.sq)
    lambda.opt <- soln$minimum
    
    df <- n - p
    Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
    Ve.opt <- lambda.opt * Vu.opt
    VCs=c(Vu.opt, Ve.opt)
    return(VCs)
}
#---------------------------------------------------------------------------------------------------------------

doEigenA_forMM=function(pheno.scaled,A ,X=NULL ) {
        n=nrow(pheno.scaled)
        if(is.null(X) ) {  X = matrix(rep(1, n), n, 1); p=1 } else {p=ncol(X) }
        XtX = crossprod(X, X)
        XtXinv = solve(XtX)
        S = diag(n) - tcrossprod(X %*% XtXinv, X)
        SHbS = S %*% A %*% S
        SHbS.system = eigen(SHbS, symmetric = TRUE)
        theta = SHbS.system$values[1:(n - p)] 
        Q = SHbS.system$vectors[, 1:(n - p)]
        return(list(theta=theta, Q=Q))
        }

#Extract info about transcripts----------------------------------------------------------------------------------
buildGeneAnnotationDF=function(transcript.info, raw.count.matrix, counts.list)  {
    annot.matrix=do.call('rbind', strsplit(transcript.info[rownames(raw.count.matrix)], ':'))
    name=rownames(raw.count.matrix)
    chr=paste0('chr', annot.matrix[,4])
    start=as.numeric(annot.matrix[,5])
    end=as.numeric(annot.matrix[,6])
    strand=as.numeric(gsub(' gene', '', as.character(annot.matrix[,7])))
    length=as.numeric(counts.list[[1]]$length)
    eff_length=as.numeric(counts.list[[1]]$eff_length)
    annot.df=data.frame(name,chr, start,end, strand, length, eff_length, stringsAsFactors=F) 
    return(annot.df)
}
#---------------------------------------------------------------------------------------------------------------------



# Build annotation table from sample names ----------------------------------------------------------------------------------
buildSampleAnnotationDF=function( count.matrix) {
    sinfo=colnames(count.matrix)
    annot.matrix=do.call('rbind', strsplit(sinfo, '-'))
    name=as.character(annot.matrix[,1])
    source.plate=as.character(annot.matrix[,2])
    source.well=as.character(annot.matrix[,3])
    growth.batch=as.character(annot.matrix[,4])
    growth.well=as.character(annot.matrix[,5])
    #seq.lane=as.character(annot.matrix[,6])
    #annot.df=data.frame(name, source.plate, source.well, growth.batch, growth.well, seq.lane,stringsAsFactors=F)
    annot.df=data.frame(name, source.plate, source.well, growth.batch, growth.well, stringsAsFactors=F)
    return(annot.df)
}
#---------------------------------------------------------------------------------------------------------------------

# function to downsample reads ------------------------------------------------------------------------------------------
downsampleCounts=function(count.matrix, downsample.total, with.replacement=FALSE) {
    downsample.count.matrix=matrix(0, nrow=nrow(count.matrix), ncol=ncol(count.matrix))
    rownames(downsample.count.matrix)=rownames(count.matrix)
    colnames(downsample.count.matrix)=colnames(count.matrix)
    total=ncol(count.matrix)
    # create progress bar
    pb =txtProgressBar(min = 1, max = ncol(count.matrix), style = 3)
    for(i in 1:ncol(count.matrix)){
        setTxtProgressBar(pb, i)
        
        if(with.replacement){
        dsgenes=sample(rep(rownames(count.matrix), as.integer(count.matrix[,i])), downsample.total,replace=T)
        
        }else{
        dsgenes=sample(rep(rownames(count.matrix), as.integer(count.matrix[,i])), downsample.total)
        }
        dcounts=plyr::count(dsgenes)
        downsample.count.matrix[match(as.character(dcounts$x), rownames(count.matrix)),i]=dcounts$freq
    }
    close(pb)
    return(downsample.count.matrix)
}
#--------------------------------------------------------------------------------------------------------------------------



#‘A’ = maximum growth value; 
#‘mu’, maximum slope; 
#‘lambda’ = lag duration
#‘integral’= integral under growth curve.
#time is in hours
# wrapper function to extract growth curve info
getGrowthCurve = function(x,gtime ,plotme=F) {

       temp=gcFitSpline(gtime, x)
  #     temp=gcFitModel(gtime, x)
       if(plotme==T) {plot(temp); readline() } 
       return(list(A=temp$parameters$A[[1]],
                   mu=temp$parameters$mu[[1]],
                   lambda=temp$parameters$lambda[[1]],
                   integral=temp$parameters$integral[1]))
}

# calculate relatedness matrix from genomewide marker data
# also, in this case, this is equivalent to (gdata %*% t(gdata))


calcA=function(p,A,do.print=T) {
    vcA.e=cbind(rep(NA, ncol(p)), rep(NA, ncol(p)))
    rownames(vcA.e)=colnames(p)
    eigA.e=doEigenA_forMM(p,A)
    # calculate mixed model, one term for additive variance  -------------------------------------------
    vcA.e=foreach(i=1:ncol(p), .combine='rbind') %dopar% {
           if(is.na(sd(p[,i]))) {
               return(c(NA,NA))
         # next;
          }
        #vcA.e[i,]=m.S(p[,i], K=A,  theta=eigA.e$theta, Q=eigA.e$Q)
        return(m.S(p[,i], K=A,  theta=eigA.e$theta, Q=eigA.e$Q))
        #if(do.print){     print(i)}
    }
    return(vcA.e)
}



extract.rr=function(x) {  return( list(sigma=x$sigma, sigma.cov=x$sigma.cov, llik=x$llik)  ) }

















# CLUSTER ANALYSIS ON TPMs -------------------------------------------------------------
##pearson tpm
#c=cor(tpm.matrix, method='pearson')
#d=dist(t(tpm.matrix))
#h2=hclust(d)
#hcd=c[h2$order,h2$order]
#X11()
#image.plot(hcd, main='pearson tpm')
#
##spearman tpm
#c=cor(tpm.matrix, method='spearman')
#d=dist(t(tpm.matrix))
#h2=hclust(d)
#hcd=c[h2$order,h2$order]
#X11()
#image.plot(hcd, main='spearman tpm')
#
##pearson batch
#c=cor(t(residuals.tpm.gbatch.OD), method='pearson')
#d=dist(residuals.tpm.gbatch.OD)
#h2=hclust(d)
#hcd=c[h2$order,h2$order]
#X11()
#image.plot(hcd, main='pearson batch')
#
##spearman batch
#c=cor(t(residuals.tpm.gbatch.OD), method='spearman')
#d=dist(residuals.tpm.gbatch.OD)
#h2=hclust(d)
#hcd=c[h2$order,h2$order]
#X11()
#image.plot(hcd, main='spearman batch')
#
##pearson batch scale
#c=cor(t(pheno.scaled), method='pearson')
#d=dist(pheno.scaled)
#h2=hclust(d)
#hcd=c[h2$order,h2$order]
#X11()
#image.plot(hcd,main='pearson batch scale')
#
##spearman batch scale
#c=cor(t(pheno.scaled), method='spearman')
#d=dist(pheno.scaled)
#h2=hclust(d)
#hcd=c[h2$order,h2$order]
#X11()
#image.plot(hcd,main='spearman batch scale')

#------------------------------------------------------------------------------------------





# Heritability analysis and heritability vs counts

#h2=vcA.downsample[[1]]
#par(xaxs='i', yaxs='i')
#plot(h2, H2, xlim=c(0,1), ylim=c(0,1))
##-------------------------------------------------------------------------------------------------------------------------------------------
#
#vcA.downsample=list()
#for(pheno in names(pheno.list)){
#    #pheno.scaled=scale(pheno.list[[pheno]])
#    #pheno.scaled=(scale(t.tpm.matrix))
#    # Comment out for subsetting
#    #lm.tpm.gbatch.OD=lm(pheno.list[[pheno]]~gbatch.fact+OD.cov)
#    lm.tpm.gbatch.OD=lm((pheno.list[[pheno]])~gbatch.fact)
#    residuals.tpm.gbatch.OD=residuals(lm.tpm.gbatch.OD)
#    pheno.scaled=(scale(residuals.tpm.gbatch.OD))
#    A=tcrossprod(gdata.scaled)/ncol(gdata.scaled)
#vcA=cbind(rep(NA, ncol(pheno.scaled)), rep(NA, ncol(pheno.scaled)))
#rownames(vcA)=colnames(pheno.scaled)
#eigA=doEigenA_forMM(pheno.scaled,A)
#    # calculate mixed model, one term for additive variance  -------------------------------------------
#for(i in 1:ncol(pheno.scaled)){
#        if(is.na(sd(pheno.scaled[,i]))) {
#        next;
#        }
#        vcA[i,]=m.S(pheno.scaled[,i], K=A,  theta=eigA$theta, Q=eigA$Q)
#     print(i)
#    }
#h2=(vcA[,1]/(vcA[,1]+vcA[,2]))
#    #print(median(h2))
#    #vcA.downsample[[pheno]]=h2



#vcA.mat= do.call('rbind', vcA.downsample)
#g.counts=apply(count.matrix[-invariant.inds,], 1, sum, na.rm=T)

#h2.med=sapply(vcA.downsample, median, na.rm=T)
#reads=c(3053776,1e6,5e5,1e5,5e4)

#png(file=paste0('/data/eQTL/plots/h2_vs_abundance_new.png'), width=1024, height=1024)
#smoothScatter(log10(g.counts), vcA.mat[1,], xlim=c(1,8), ylim=c(0,1), nbin=256, 
#              colramp=colorRampPalette(c('white', blues9, 'purple')),
#              nrpoints=7000, xlab='log10(counts per transcript)', ylab='h^2', main='h^2 vs abundance', pch=21, cex=.5, col='#00000066')
#dev.off()

#scatter.smooth(log10(g.counts), vcA.mat[1,],xlim=c(0,10), ylim=c(0,1), col='#00000066', cex=.7)


#png(file=paste0('/data/eQTL/plots/h2_vs_downsampling_new.png'), width=512, height=512)
#par(yaxs='i', xaxs='i')
#plot(reads/1e6, h2.med, ylim=c(0,.4), xlim=c(0,5),ylab='h^2', xlab='million reads per sample', main='h^2 vs downsampling', type='p', cex=2, pch=20)
#dev.off()
#
#
#read.bins=rle(as.character(sort(bins)))
#read.bins=paste(read.bins$values, read.bins$lengths, sep=' n=')
#png(file=paste0('/data/eQTL/plots/h2_vs_downsampling_by_abundance_class_new.png'), width=1920, height=1024)
#
#boxplot(vcA.mat[1,]~bins, border='black',at=c(((1:10)*5)-4), xlim=c(1,51), names=read.bins, xlab='observed log10(reads per transcript)', ylab='h^2')
#boxplot(vcA.mat[2,]~bins, border='blue',at=c(((1:10)*5)-3), xlim=c(1,50), add=T,names=rep('',10))
#boxplot(vcA.mat[3,]~bins, border='purple',at=c(((1:10)*5)-2), xlim=c(1,50), add=T,names=rep('',10))
#boxplot(vcA.mat[4,]~bins, border='green',at=c(((1:10)*5)-1), xlim=c(1,50), add=T,names=rep('',10))
#boxplot(vcA.mat[5,]~bins, border='red',at=c(((1:10)*5)), xlim=c(1,50), add=T,names=rep('',10))
#legend('topleft', legend=round(log10(reads),2),
#       title='downsampling (log10(reads per sample))',
#       col=c('black', 'blue', 'purple', 'green', 'red'),
#       fill=c('black', 'blue', 'purple', 'green', 'red'))
#dev.off()
#
#
#bins=cut(log10(as.vector(g.counts)),10)


# wl3=which(LODS>3, arr.ind=T)
##png(file=paste0('/data/eQTL/plots/LOD_subset_', ss, '.png'), width=768, height=768)
#par(xaxs='i', yaxs='i')
#    plot(wl3[,2],wl3[,1], xlab='marker position (marker genome index)', ylab='gene position (gene genome index)', pch=19, col='#00000022', cex=.5,
#         main=paste('n=', nrow(pheno.scaled), '   BYxRM eQTL', ss),)
#    abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='grey')
#    abline(h=cumsum(rle(gene.annot.df$chr)$lengths), col='grey')
##dev.off()
#LODSb=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled,betas=T)
#save(LODSb, file='/data/eQTL/RData/LODSb.tpm.bc.RData')
#chromosome vector for markers


# batch, A, batch X A, AA, batch X AA VC decomposition 
#gbe=list()
#for(pheno in rownames(tpm.matrix)[-c(1:1544)]){
#    print(pheno)
#    print( paste(round(match(pheno, rownames(tpm.matrix))/nrow(tpm.matrix)*100,2), '%'))
#    p=tpm.matrix[pheno,]
#    if(is.na(sd(p))) {next;}
#    #rr=regress(p~1, ~batch.mat+A+batch.by.g+AA+batch.by.gg, verbose=F)
#    rr=regress(p~1, ~batch.mat+A+batch.by.g, verbose=F)
#    print(summary(rr))
#    gbe[[pheno]]=rr$sigma
#}
#y2=sapply(gbe, function(x) x/sum(x))
#y=sapply(gbe, function(x) x/sum(x))
#vpm=((apply(y, 1, mean)[1:4])/sum(apply(y, 1, mean)[1:4]))
#vps=((apply(y, 1, sd)[1:4])/sum(apply(y, 1, mean)[1:4]))
#
#bp=barplot(vpm, las=1, ylim=c(0,.45), ylab='fraction of phenotypic variance',
#           names=c('batch', 'A', 'batch X A', 'E')
#           )
#segments(bp,      vpm-vps/sqrt(6623), bp,    vpm+vps/sqrt(6623), lwd=1.5, col='black')
#sum((apply(y, 1, mean)[2:4])/sum(apply(y, 1, mean)[2:4])) 
#save(y, file='/data/eQTL/RData/batch_mixed_effect_model.RData')


build.gcoord.key =function(filein) {
    ref.seq=read.fasta(filein)
    contig.lengths=sapply(ref.seq, length)
    names(contig.lengths)=paste0('chr', c(as.character(as.roman(1:16)), 'Mito'))
    gcoord.key=cumsum(c(0, as.vector(contig.lengths)))[-18]
    names(gcoord.key)=names(contig.lengths)
    return(gcoord.key)
}


get.contig.lengths =function(filein) {
    ref.seq=read.fasta(filein)
    contig.lengths=sapply(ref.seq, length)
    return(contig.lengths)
}


# t.tpm.matrix (normalized tpm measurements)
# pheno.scaled (residuals from tpm + batch + additional covariates)
# gdata (genotypes, unscaled)
# gdata.scaled (genotypes, scaled)
# additional.cov (are there additional covariates)
# covs (the additional covariates)




find.background.QTL = function (covariates, t.tpm.matrix, pheno.scaled, gdata, gdata.scaled ,bck.thresh=3.5,
                                chromosomes=NULL) {

    if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
    cvec=(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])

    # Iteration 1  (find background QTL effects)
    LODS=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
    # lists to hold LODs, markers, and scaled markers by chromosome----------------
    LODS.by.chr=list()
    gdata.s.by.chr=list()
    gdata.by.chr=list()
    for(cc in chromosomes) {   
            LODS.by.chr[[cc]]=LODS[,which(cvec %in% cc)]   
            gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
            gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
    }
    #-----------------------------------------------------------------------------
     
    # max LOD per chromosome
    max.lod.chr=lapply(LODS.by.chr, function(x) { apply(x,1,max) })
    # index of max LOD per chromosome
    max.lod.chr.ind=lapply(LODS.by.chr, function(x) { apply(x,1,which.max) })

    # pulls out the marker index (by chromosome) of the peak marker 
    pmarker.chr.ind =  mapply(function(x,y,z) { 
        h=y[x>bck.thresh]
        i=colnames(z)[h]
        names(i)=names(h)
        return(i)
        },  max.lod.chr, max.lod.chr.ind,gdata.by.chr)

    pmarker.chr.ind.df=data.frame(chr=rep(names(pmarker.chr.ind), sapply(pmarker.chr.ind, length)),
             gene=as.vector(do.call('c', lapply(pmarker.chr.ind, function(x) names(x)))),
             marker=as.vector(do.call('c', lapply(pmarker.chr.ind, function(x) (x)))), stringsAsFactors=F)
    pbg=split(pmarker.chr.ind.df$marker, pmarker.chr.ind.df$gene)

    # Iteration 2--------------------------------------------------------------------------------------------
    plist=lapply(rownames(tpm.matrix), function(i) { 
         if(!(i %in% names(pbg)) ) {  return(pheno.scaled[,i]) } 
         else{ residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,pbg[[i]]]) )  }
    })
    p.resid1=do.call('cbind', plist)
    colnames(p.resid1)=colnames(pheno.scaled)
    pheno.scaled=scale(p.resid1)
    LODS1=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
    LODS1.by.chr=list()
    for(cc in chromosomes) { LODS1.by.chr[[cc]]=LODS1[,which(cvec %in% cc)]   }
    max.lod1.chr = lapply(LODS1.by.chr, function(x) { apply(x,1,max) })
    max.lod1.chr.ind = lapply(LODS1.by.chr, function(x) { apply(x,1,which.max) })
    pmarker1.chr.ind = mapply(function(x,y,z) { 
        h=y[x>bck.thresh]
        i=colnames(z)[h]
        names(i)=names(h)
        return(i)
        },  max.lod1.chr, max.lod1.chr.ind,gdata.by.chr)
    pmarker1.chr.ind.df=data.frame(chr=rep(names(pmarker1.chr.ind), sapply(pmarker1.chr.ind, length)),
             gene=as.vector(do.call('c', lapply(pmarker1.chr.ind, function(x) names(x)))),
             marker=as.vector(do.call('c', lapply(pmarker1.chr.ind, function(x) (x)))), stringsAsFactors=F)
    pbg1=split(pmarker1.chr.ind.df$marker, pmarker1.chr.ind.df$gene)
    pbg.merge=pbg
    for(i in names(pbg1)){    pbg.merge[[i]]=c(pbg.merge[[i]], pbg1[[i]]) }

    # Iteration 3--------------------------------------------------------------------------------------------
    plist=lapply(rownames(tpm.matrix), function(i) { 
         if(!(i %in% names(pbg.merge)) ) {  return(pheno.scaled[,i])} 
         else{ residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,pbg.merge[[i]]]) )    }
    })
    p.resid2=do.call('cbind', plist)
    colnames(p.resid2)=colnames(pheno.scaled)
    pheno.scaled=scale(p.resid2)
    LODS2=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
    LODS2.by.chr=list()
    for(cc in chromosomes) {    LODS2.by.chr[[cc]]=LODS2[,which(cvec %in% cc)]   }
    max.lod2.chr=lapply(LODS2.by.chr, function(x) { apply(x,1,max) })
    max.lod2.chr.ind=lapply(LODS2.by.chr, function(x) { apply(x,1,which.max) })
    pmarker2.chr.ind =  mapply(function(x,y,z) { 
        h=y[x>bck.thresh]
        i=colnames(z)[h]
        names(i)=names(h)
        return(i)
        },  max.lod2.chr, max.lod2.chr.ind,gdata.by.chr)
    pmarker2.chr.ind.df=data.frame(chr=rep(names(pmarker2.chr.ind), sapply(pmarker2.chr.ind, length)),
             gene=as.vector(do.call('c', lapply(pmarker2.chr.ind, function(x) names(x)))),
             marker=as.vector(do.call('c', lapply(pmarker2.chr.ind, function(x) (x)))), stringsAsFactors=F)
    pbg2=split(pmarker2.chr.ind.df$marker, pmarker2.chr.ind.df$gene)

    pbg.final=pbg.merge
    for(i in names(pbg2)){    pbg.final[[i]]=c(pbg.final[[i]], pbg2[[i]]) }

    pm.chr.ind.df=rbind(pmarker.chr.ind.df, pmarker1.chr.ind.df, pmarker2.chr.ind.df)
    background.QTL=list()
    for(cc in chromosomes) {
        background.QTL[[cc]]=split(pm.chr.ind.df$marker[pm.chr.ind.df$chr!=cc], pm.chr.ind.df$gene[pm.chr.ind.df$chr!=cc]) 
    }
    return(background.QTL)
}





mapQTL=function(covariates, background.QTL, t.tpm.matrix, pheno.scaled, 
                gdata, gdata.scaled, 
                chromosomes=NULL,
                n.perm=1000, FDR.thresh=.05, Arandom=TRUE) {

    if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
    cvec=(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])

    gdata.s.by.chr=list()
    gdata.by.chr=list()
    for(cc in chromosomes) {   
            gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
            gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
    }

    peakList=list()
    #LODmatrix=list()
    for(cc in chromosomes) {  
         print(cc) 
         plist=lapply(colnames(t.tpm.matrix), function(i) { 
              if(!(i %in% names(background.QTL[[cc]])) ) { return(  residuals(lm(t.tpm.matrix[,i]~covariates)) ) }
              else{  residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,background.QTL[[cc]][[i]]]) ) }
         })
         preal=do.call('cbind', plist)
         colnames(preal)=colnames(t.tpm.matrix)
         
         presid=preal

         # remove residual additive effect
         if(Arandom) {
            Aloco=A.mat(do.call('cbind', (gdata.by.chr[-match(cc, names(gdata.by.chr))])))/2
            eigA=doEigenA_forMM(presid,Aloco)
            svdAloco=svd(Aloco)
            pb=txtProgressBar(min=1, max=ncol(presid), style=3)
            for(tp in 1:ncol(presid)){
                 setTxtProgressBar(pb,tp)
                rr=m.S(presid[,tp], K=Aloco,  theta=eigA$theta, Q=eigA$Q)
                ###W=solve(rr[1]*Aloco+rr[2]*diag(nrow(presid)))
                W=svdAloco$u %*% tcrossprod(diag(1/((svdAloco$d*rr[1])+(rr[2]))), svdAloco$v)
                if(rr[1]>0) {
                      blups=calc.BLUPS(rr[1]*Aloco,diag(nrow(presid)),W,presid[,tp],matrix(1,nrow(presid),1),0 )[,1]
                      presid[,tp]=as.vector(presid[,tp] - blups)
                }
            }
            rm(W)
            close(pb)
         }
         
         preal=presid
         presid=scale(presid)
         
         #pscale=pheno.scaled
         obsLOD = fasterLOD(nrow(preal),presid,gdata.s.by.chr[[cc]],betas=TRUE,pheno=preal)
         findingPeaks=TRUE
         jump=1

        while(findingPeaks==TRUE) {
            max.obsLOD=apply(obsLOD$LOD,1,max)   
            max.ind.obsLOD=apply(obsLOD$LOD,1,which.max)   

            permLOD=matrix(0, length(max.obsLOD), n.perm)
           
            print('doing permutations')
            pb =txtProgressBar(min = 1, max =n.perm, style = 3)
            for(i in 1:n.perm) {
                 setTxtProgressBar(pb, i)
                 permLOD[,i]=apply(fasterLOD(nrow(presid),presid[sample(nrow(presid)),], gdata.s.by.chr[[cc]]),1,max)
            }
           close(pb)
           
           obsPcnt = sapply(seq(1.5, 9, .05), function(thresh) { sum(max.obsLOD>thresh) }   )
           names(obsPcnt) = seq(1.5, 9, .05)
           if(sum(obsPcnt)<3) {break}
           # expected number of QTL peaks with LOD greater than threshold
           expPcnt = sapply(seq(1.5, 9, .05),  
                             function(thresh) { 
                                    #print(thresh); 
                                    mean(apply(permLOD, 2, function(ll) {sum(ll>thresh) }) )
                                } )
          names(expPcnt) = seq(1.5, 9, .05)
          pFDR = expPcnt/obsPcnt
          pFDR[is.na(pFDR)]=0
          pFDR[!is.finite(pFDR)]=0
          #to make sure this is monotonic
          pFDR = rev(cummax(rev(pFDR)))
          fdrFX=approxfun(pFDR, seq(1.5,9,.05))
          thresh=fdrFX(FDR.thresh)
          print(paste('FDR Thresh 5%', thresh)) 
          peak.gene.ind=max.ind.obsLOD[which(max.obsLOD>thresh)]

          if(length(peak.gene.ind)>0) { 
              findingPeaks=TRUE 
              pldf=data.frame(gene=names(peak.gene.ind),
                              pmarker=colnames(gdata.by.chr[[cc]])[peak.gene.ind], 
                              pcind=as.vector(peak.gene.ind),
                              stringsAsFactors=F)
              pldf$r=obsLOD$r[cbind(pldf$gene, pldf$pmarker)]
              pldf$LOD=obsLOD$LOD[cbind(pldf$gene, pldf$pmarker)]
              lodInt=matrix(NA,nrow(pldf),2)
              #Lmat=obsLOD$LOD[pldf$gene,]
              for(peak in 1:nrow(pldf)){
                 int.range=range(which(obsLOD$LOD[pldf[peak,]$gene,]>pldf[peak,]$LOD-1.5)) 
                 lodInt[peak,]=colnames(gdata.by.chr[[cc]])[int.range]
                 #if(is.null(nrow(Lmat))) {Lmat[!(c(1:length(Lmat)) %in% seq(int.range[1], int.range[2]))]=0 }
                 #else {Lmat[peak, !(c(1:ncol(Lmat)) %in% seq(int.range[1], int.range[2]))]=0 }
              }
              pldf$CI.l=lodInt[,1]
              pldf$CI.r=lodInt[,2]
              peakList[[cc]][[jump]]=pldf
              
              #LODmatrix[[cc]][[jump]]=Lmat
              presid.tmp=matrix(0, nrow(t.tpm.matrix), nrow(pldf))
              rownames(presid.tmp)=rownames(t.tpm.matrix)
              colnames(presid.tmp)=pldf$gene
             
              pldfm = do.call('rbind', peakList[[cc]])
              found.q.covs=split(pldfm$pmarker, pldfm$gene)
              
              for(i in pldf$gene){
                  # modified 041017
                  #if(!(i %in% names(background.QTL[[cc]])) ) {
                  #    presid.tmp[,i]=scale(residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,found.q.covs[[i]]]))) }
                  #else{
                  #   presid.tmp[,i]=scale(residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,background.QTL[[cc]][[i]]]+gdata[,found.q.covs[[i]]])) ) }
                  presid.tmp[,i]=scale(residuals(lm(presid[,i]~gdata[,found.q.covs[[i]]])))
              }
              obsLOD = fasterLOD(nrow(presid.tmp),presid.tmp,gdata.s.by.chr[[cc]],betas=TRUE, pheno=presid.tmp)
              presid=presid.tmp
              jump=jump+1 
          } else {findingPeaks=FALSE} 
    }
    }
    return(peakList)
}




buildPeakListDF=function(peakList.OD, gdata, gene.GR, marker.GR, peak.dist.cis=40000) {
    
    #use Frank's definition of padded gene (1000 upstream, 200 downstream) --------------------------------
    g2=gene.GR
    g2p=g2[strand(g2)=='+']
    g2m=g2[strand(g2)=='-']
    start(ranges(g2p))=start(ranges(g2p))-1000
    end(ranges(g2p))=end(ranges(g2p))+200
  
    start(ranges(g2m))=start(ranges(g2m))-200
    end(ranges(g2m))=end(ranges(g2m))+1000
    
    g2[strand(g2)=='+']=g2p
    g2[strand(g2)=='-']=g2m

    gene.GR=g2
    #---------------------------------------------------------------------------------------------------------

    all.peaks=do.call('rbind', lapply(peakList.OD, function(x) do.call('rbind', x)))
    all.peaks.g.C=    match(all.peaks$pmarker, colnames(gdata))
    all.peaks.g.L=    match(all.peaks$CI.l, colnames(gdata))
    all.peaks.g.R=    match(all.peaks$CI.r, colnames(gdata))
    intv=Intervals(cbind(all.peaks.g.L, all.peaks.g.R))
    intv=intv[order(all.peaks.g.C),]
    #plot(intv,  use_points=F, lwd=.3, col='#00000066', ylim=c(0,3500))
    poverlap.cnt=interval_overlap(1:ncol(gdata),intv)

    #png(file='/data/eQTL/plots/hotspot_interval_overlap.png', width=1080, height=512)
    plot(sapply(poverlap.cnt, length), type='h', col='grey', ylab='transcripts "regulated"')
    abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])$lengths), col='grey')
    #dev.off()

    pgi=match(all.peaks$gene, colnames(t.tpm.matrix))
    mgi=match(all.peaks$pmarker, colnames(gdata))
    mgiL=all.peaks.g.L
    mgiR=all.peaks.g.R

    peaks.gene.GR=gene.GR[match(all.peaks$gene, gene.GR$ORF)]
    peaks.marker.GR=marker.GR[match(all.peaks$pmarker, marker.GR$mname)]

    peakInt.GR=GRanges(seqnames= seqnames(marker.GR[mgiL]), ranges=IRanges(start=start(marker.GR[mgiL]), end=end(marker.GR[mgiR])))
    
    all.peaks$marker.gcoord=peaks.marker.GR$gcoord
    all.peaks$gene.gcoord=peaks.gene.GR$gcoord
           
    all.peaks$gene.to.CI.dist=distance(peaks.gene.GR, peakInt.GR)
    all.peaks$cis=(all.peaks$gene.to.CI.dist<1) & abs(all.peaks$gene.gcoord-all.peaks$marker.gcoord)<peak.dist.cis
    all.peaks$cis[is.na(all.peaks$cis)]=FALSE
    all.peaks$chr=sapply(strsplit(all.peaks$pmarker, ':'), function(x) x[1] )
    all.peaks$gcind=match(all.peaks$pmarker, colnames(gdata))
    return(all.peaks)
}


eQTL_bigPlot=function(all.peaks, gcoord.key, marker.GR, xlim.ind=NULL) {
    no.coord=F
    if(is.null(xlim.ind) ) {xlim.ind=1:length(marker.GR) 
    no.coord=T
    }
    with(all.peaks, { 
        par(xaxs='i', yaxs='i')
        if(no.coord) {
        plot(marker.gcoord,gene.gcoord, pch=20, cex=.5, type='n', xlab='marker position', ylab='transcript position'
             ,xaxt='n', yaxt='n', xlim=range( marker.GR$gcoord[xlim.ind])) 
        axis(1, at=gcoord.key, labels=names(gcoord.key))
        }
        else{
            
        plot(marker.gcoord,gene.gcoord, pch=20, cex=.5, type='n', xlab='marker position', ylab='transcript position'
             , yaxt='n',xaxt='n', xlim=range( marker.GR$gcoord[xlim.ind])) 
        # new for frank this was commented out
         axis(1, at=gcoord.key, labels=names(gcoord.key))
        }
        axis(2, at=gcoord.key, labels=names(gcoord.key))

        purp='#80008055'
        orng='#FFA50055'
               points(marker.gcoord[r>0],gene.gcoord[r>0], pch=20, cex=.5, col='#80008088')
        points(marker.gcoord[r<0],gene.gcoord[r<0], pch=20, cex=.5, col='#FFA50088')
        points(marker.gcoord[cis],gene.gcoord[cis], pch=20, cex=.5, col='red')
        #abline(v=match(as.vector(unlist(hotspots)), colnames(gdata)), lty=3, col='blue')
        abline(v=gcoord.key, col='grey' )
        abline(h=gcoord.key, col='grey' )
       
        #new for frank   
           # add intervals
        #   mgiL=marker.GR$gcoord[match(CI.l, marker.GR$mname)]
        #   mgiR=marker.GR$gcoord[match(CI.r, marker.GR$mname)]
        ##   segments(mgiL[r>0],gene.gcoord[r>0], mgiR[r>0],gene.gcoord[r>0], lwd=.5, col='#80008055')
        ##   segments(mgiL[r<0],gene.gcoord[r<0], mgiR[r<0],gene.gcoord[r<0], lwd=.5, col='#FFA50055')
          })
        #new for frank
        #gxl=(marker.GR$gcoord[xlim.ind])
        #cxl=start(marker.GR[xlim.ind])
        #gtoc=approxfun(cxl, gxl, rule=2)
        #pcxl=pretty(cxl)
        #gcxl=gtoc(pcxl)
        #axis(1, at =gcxl , labels=pcxl)

       #segments(mgiL[all.peaks$r>0],pgi[all.peaks$r>0], mgiR[all.peaks$r>0],pgi[all.peaks$r>0], lwd=.5, col='#80008055')
       #segments(mgiL[all.peaks$r<0],pgi[all.peaks$r<0], mgiR[all.peaks$r<0],pgi[all.peaks$r<0], lwd=.5, col='#FFA50055')
       #, xlim=c(3700,4500))#, col='grey') ##, xlim=c(3750,5500))
       # segments(mgiL,pgi, mgiR,pgi, lwd=.5, col='#00000055' )
       #, xlab='marker index', ylab='transcript index')#, col='grey') ##, xlim=c(3750,5500))
       #cumsum(rle(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])$lengths), col='grey')
       #abline(h=cumsum(rle(gene.annot.df$chr)$lengths), col='grey')
       #})
}

buildGeneticMap=function(gdata.by.chr) {
    genetic.map=list()
    for(i in 1:16){
        print(i)
        nmap=rep(0, ncol(gdata.by.chr[[i]]))
        names(nmap)=colnames(gdata.by.chr[[i]])
        nmap[1]=0
        for( j in 2:(ncol(gdata.by.chr[[i]]))) {
            nmap[j]=100*((1012-sum(gdata.by.chr[[i]][,j]==gdata.by.chr[[i]][,j-1]))/1012)
        }
        nmap=cumsum(nmap)
        genetic.map[[paste0('chr', as.roman(1:16))[i]]]=nmap
    }
 return(genetic.map) 
}




# code for hotspot analysis ------------------------------------------------------------------------------------------
model.QTL.effects = function( t.tpm.matrix, gdata,  keep.transcripts, peaks.per.gene, gbatch.fact, pmarker=NULL, cc, docoef=FALSE) {
       tmm=matrix(0, nrow(t.tpm.matrix), length(keep.transcripts))
       mid=match(keep.transcripts, colnames(t.tpm.matrix))
       colnames(tmm)=colnames(t.tpm.matrix)[mid]
       rownames(tmm)=rownames(t.tpm.matrix)
       #Correct transcripts for background QTL effects 
       q.cnt=rep(0, length(keep.transcripts))
       q.pos=list()
       names(q.cnt)=colnames(tmm)
       pb =txtProgressBar(min = 1, max =length(keep.transcripts), style = 3)
       for(i in 1:length(keep.transcripts)) { 
             setTxtProgressBar(pb, i)
           
             mid.i= mid[i]
             gene.i=colnames(tmm)[i]
             ppgi.a=peaks.per.gene[[gene.i]]
             q.cnt[i]=sum(ppgi.a$chr==cc)
             q.pos[[gene.i]]=ppgi.a$pcind[ppgi.a$chr==cc]
             ppgi=ppgi.a[ppgi.a$chr!=cc,]
             if(is.null(ppgi)) {
                 if(length(pmarker>1)){ 
                        l1= lm(t.tpm.matrix[,mid.i]~gbatch.fact+gdata[,pmarker])
                 } else {
                     l1=lm(t.tpm.matrix[,mid.i]~gbatch.fact)
                 }
             }else {
                 if(nrow(ppgi)==0) {
                     if(length(pmarker>1)){ 
                         l1=lm(t.tpm.matrix[,mid.i]~gbatch.fact+gdata[,pmarker])
                     } else {
                         l1=lm(t.tpm.matrix[,mid.i]~gbatch.fact)
                     }      
                 #next;
                 } else {
                     bQTL=gdata[,ppgi$gcind[ppgi$chr!=cc]]
                     if(length(pmarker>1)) {
                         #print(i)
                          l1=lm(t.tpm.matrix[,mid.i]~gbatch.fact+bQTL+gdata[,pmarker])

                     } else{
                         l1=lm(t.tpm.matrix[,mid.i]~gbatch.fact+bQTL)

                     }
                 }
             }
       tmm[,gene.i]=scale(residuals(l1))
       }
       close(pb)
       return(tmm)
}

# tmm is residualized phenotype
# max.obsLOD is observed peak LOD per transcript
do.hotspot.FDR = function( tmm, max.obsLOD, g.s, n.perm=100, fdr.thresh=.05) {
   permLOD=matrix(0, length(max.obsLOD), n.perm)
   print('doing permutations')
   pb =txtProgressBar(min = 1, max =n.perm, style = 3)
   for(i in 1:n.perm) {
             setTxtProgressBar(pb, i)
             permLOD[,i]=apply(fasterLOD(nrow(tmm),tmm[sample(1:nrow(tmm)),], g.s),1,max)
   }
   close(pb)
       
   obsPcnt = sapply(seq(1.5, 6, .01), function(thresh) { sum(max.obsLOD>thresh) }   )
   names(obsPcnt) = seq(1.5, 6, .01)
       
   if(sum(obsPcnt)<5) { return(999) } #break;}
               
   # expected number of QTL peaks with LOD greater than threshold
   expPcnt = sapply(seq(1.5, 6, .01),  
                         function(thresh) { 
                                #print(thresh); 
                                mean(apply(permLOD, 2, function(ll) {sum(ll>thresh) }) )
                            } )
   names(expPcnt) = seq(1.5, 6, .01)
   pFDR = expPcnt/obsPcnt
   pFDR = rev(cummax(rev(pFDR)))
   if(sum(is.na(pFDR))>2) {
       print('this')
       return(999) } #break;}

   fdrFX=approxfun(pFDR, seq(1.5,6,.01))
   thresh=fdrFX(fdr.thresh)
   print(paste(c('FDR Thresh', fdr.thresh, '%', thresh))) 
   return(thresh)
}

build.model.matrix=function(sub.t, g.s, l.boundary.index, r.boundary.index) {
    if(l.boundary.index == 0 ) {   mmatrix=model.matrix(lm(sub.t[,1]~ g.s[,r.boundary.index])) 
       } else if(r.boundary.index == ncol(g.s) ) { 
           mmatrix=model.matrix(lm(sub.t[,1]~ g.s[,l.boundary.index] ))  
       } else {
           mmatrix=model.matrix(lm(sub.t[,1]~g.s[, l.boundary.index] + g.s[,r.boundary.index]))
       }
    return(mmatrix)
}
   

calcMVstat=function(g.s, sub.t,  mmatrix, l.boundary.index, r.boundary.index) {
    sss=rep(0, ncol(g.s))
    sub.null=scale(residuals(lm(sub.t~mmatrix)))
    SST=crossprod(sub.null)
    t5=abs(eigen(SST)$values)
    pc.to.retain=which(cumsum(t5/sum(t5))<.999)
    pb =txtProgressBar(min = l.boundary.index, max =r.boundary.index, style = 3)
    nn=ncol(sub.t)
    for(i in l.boundary.index:r.boundary.index) { 
        setTxtProgressBar(pb, i)
        test=lm.fit(cbind(mmatrix,g.s[,i]), sub.t)
        t2=residuals(test)
        SSE=crossprod(t2)
        t4=abs(eigen(SSE)$values)
        X2=-(1012-1-((nn+2)/2))*(sum(log(t4[pc.to.retain]))-sum(log(t5[pc.to.retain])))
        sss[i]=X2
     } 
    close(pb) 
    names(sss)=colnames(g.s)
    return(sss)
}


findHotspots.iteration1=function(t.tpm.matrix, gene.annot.df, peaks.per.chr, genetic.map,
                                 gdata, gdata.s.by.chr, add.cov, do.save=T) { #gbatch.fact) {
    hotspots=list()
    total.hotspot.n=0
    unique.chrs=paste0('chr', as.roman(1:16))
    for(cc in unique.chrs) { 
       hs.n=1
       ppc=peaks.per.chr[[cc]]
       print(cc)
       keep.transcripts=colnames(t.tpm.matrix)
       # remove local QTL -------------------------------------------
       keep.transcripts=keep.transcripts[gene.annot.df$chr!=cc]  
       ppc2=ppc[ppc$gene %in% keep.transcripts,]
       keep.transcripts=ppc2$gene
       keep.transcripts=unique(keep.transcripts)
       #------------------------------------------------------------- 
       g.s=gdata.s.by.chr[[cc]]
       g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
       
       pmarker=c()
       repeat {
           tmm=model.QTL.effects(t.tpm.matrix, gdata, keep.transcripts, peaks.per.gene,add.cov, pmarker,cc)
           tmp= fasterLOD(nrow(tmm),tmm,g.s)
           max.obsLOD=apply(tmp,1,max)   
           #str(tmm)

           thresh=do.hotspot.FDR(tmm, max.obsLOD, g.s)
           
           if(thresh==999) {break; }
           #thresh
           keep.transcripts=names(which(max.obsLOD>thresh))
           print(length(keep.transcripts))
           
           if(length(keep.transcripts)<45) {break;}

           LOD.kt=tmp[keep.transcripts,]
           LOD.kt.wmax=apply(tmp[keep.transcripts,], 1, which.max)
           LOD.kt.max=apply(tmp[keep.transcripts,], 1, max)
          
           l15=list()
            for(peak in 1:length(LOD.kt.max)){
                 l15[[peak]]=range(which(LOD.kt[peak,]>(LOD.kt.max[peak]-1.5) ) )
            }
           l15=Intervals(do.call('rbind', l15))
           io=interval_overlap(1:ncol(g.s), l15)
           peak.index=which.max(sapply(io, length))  
           # make this binning proportional to chromosome length
           ##h=hist(LOD.kt.wmax, breaks=contig.lengths[cc]*.00013)
           ##peak.index=floor(h$mids[which.max(h$counts)])
           if(peak.index>ncol(g.s)) {peak.index=ncol(g.s) }
           r.boundary.index=findInterval(genetic.map[[cc]][peak.index]+5, genetic.map[[cc]], rightmost.closed = TRUE )
           l.boundary.index=findInterval(genetic.map[[cc]][peak.index]-5, genetic.map[[cc]] )

           peaks.in.peak.bin=interval_overlap(peak.index, l15)[[1]]
           #which((LOD.kt.wmax>=l.boundary.index) & (LOD.kt.wmax< r.boundary.index))
           to.scan=names(LOD.kt.max)[peaks.in.peak.bin]
           
           to.scan=LOD.kt[to.scan,peak.index]
           to.scan=sort(to.scan, decreasing=T)
           if(length(to.scan)>500) { to.scan=to.scan[1:500] }
           if(length(to.scan)<3)  {break;}
           
           nn=length(to.scan)
           #to.scan2=to.scan[to.scan<30]
           sub.t=tmm[,names(to.scan)]
            
           mmatrix=build.model.matrix(sub.t, g.s, l.boundary.index, r.boundary.index)
           sss =calcMVstat(g.s, sub.t, mmatrix, l.boundary.index, r.boundary.index)      

           print(colnames(g.s)[which.max(sss)])
           if(!is.finite(max(sss)))  {break;}
           if(min(sss)< (-1))        {break;}
           if((max(sss)/4.605)<5)    {break;}
           #mvlod[[cc]][[hs.n]]=sss
           #rl.w=rle(which(sss>(max(sss)-1.5 )))
           hs.int.range=range(which((sss/4.605)>(max(sss/4.605)-1.5 ))) 
           names(sss)[hs.int.range]
           
           #plot(g.s.pos, sss, col='red',ylim=range(sss[sss>0]))
           pmarker=c(pmarker, names(which.max(sss)))
           l=names(sss)[hs.int.range][1]
           r=names(sss)[hs.int.range][2]

           hotspots[[cc]][[hs.n]]=list(maxneglog10p=max(sss), 
                                             n.trans=length(keep.transcripts),
                                             peak=names(which.max(sss)), 
                                             cI.l=l,
                                             cI.r=r,
                                             to.scan=to.scan
                                             )
           total.hotspot.n=total.hotspot.n+1
           print(cc) 
           print(hs.n)
           print(hotspots[[cc]][[hs.n]]$peak)   
           hs.n=hs.n+1
            
           pmarker.new = names(which.max(sss))
           l.boundary=g.s.pos[match(pmarker.new, colnames(g.s))]-25000
           r.boundary=g.s.pos[match(pmarker.new, colnames(g.s))]+25000
           l.boundary.index  =  findInterval(l.boundary, g.s.pos)
           r.boundary.index  =  findInterval(r.boundary, g.s.pos, rightmost.closed=F)
           if(l.boundary.index<1) l.boundary.index=1
           if(r.boundary.index>(ncol(g.s))) r.boundary.index=ncol(g.s)
            #-----------------------------------------------------------------------------------------------

           print(l.boundary.index)
           print(r.boundary.index)
            attr(sub.t, 'chr')=cc
            attr(sub.t, 'hs.n')=total.hotspot.n
            attr(sub.t, 'l.boundary.index')=l.boundary.index
            attr(sub.t, 'r.boundary.index')=r.boundary.index
            if(do.save) {
                 #save(sub.t, file=paste0('/media/juno/bootstrap_hotspot/', total.hotspot.n, '.RData' ))
                 save(sub.t, file=paste0('/data/eQTL/RData/hotspots/', total.hotspot.n, '.RData' ))
            }
      }
    }
    return(hotspots)
}

findHotspots.iteration2=function(hotspots, t.tpm.matrix, gene.annot.df, peaks.per.chr, genetic.map,
                                 gdata, gdata.s.by.chr, gbatch.fact, do.save=T) {
    hotspots.refined=list()
    total.hotspot.n=0
    unique.chrs=paste0('chr', as.roman(1:16))
    for(cc in unique.chrs) { 
        #g=gdata.by.chr[[cc]]
        g.s=gdata.s.by.chr[[cc]]
        g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
           
        pmarker.all=as.character(sapply(hotspots[[cc]], function(x)x$peak))
        pmarker.new=pmarker.all
        #peaks=c()
        #do while
        for(pm in 1:length(pmarker.all)) {
            pmarker=pmarker.new[-pm]
            keep.transcripts=names(hotspots[[cc]][[pm]]$to.scan)
            sub.t=model.QTL.effects(t.tpm.matrix, gdata, keep.transcripts, peaks.per.gene, gbatch.fact, pmarker, cc)
                  
            l.boundary=g.s.pos[match(pmarker.all[pm], colnames(g.s))]-50000
            r.boundary=g.s.pos[match(pmarker.all[pm], colnames(g.s))]+50000
            l.boundary.index  =  findInterval(l.boundary, g.s.pos)
            r.boundary.index  =  findInterval(r.boundary, g.s.pos, rightmost.closed=F)
            if(l.boundary.index<1) l.boundary.index=1
            if(r.boundary.index>(ncol(g.s))) r.boundary.index=ncol(g.s)
              
            mmatrix=model.matrix(lm(sub.t[,1]~1))
            sss=calcMVstat(g.s, sub.t, mmatrix, l.boundary.index, r.boundary.index)
            print(names(which.max(sss))) 
            hotspots.refined[[cc]][[pm]]=data.frame(maxneglog10p=max(sss), 
                                                    n.trans=length(keep.transcripts),
                                                    peak=names(which.max(sss)))
            pmarker.new[pm]=names(which.max(sss))
            total.hotspot.n=total.hotspot.n+1
               
            #update boundaries ----------------------------------------------------------------------------
            l.boundary=g.s.pos[match(pmarker.new[pm], colnames(g.s))]-50000
            r.boundary=g.s.pos[match(pmarker.new[pm], colnames(g.s))]+50000
            l.boundary.index  =  findInterval(l.boundary, g.s.pos)
            r.boundary.index  =  findInterval(r.boundary, g.s.pos, rightmost.closed=F)
            if(l.boundary.index<1) l.boundary.index=1
            if(r.boundary.index>(ncol(g.s))) r.boundary.index=ncol(g.s)
            #-----------------------------------------------------------------------------------------------

            attr(sub.t, 'chr')=cc
            attr(sub.t, 'hs.n')=total.hotspot.n
            attr(sub.t, 'l.boundary.index')=l.boundary.index
            attr(sub.t, 'r.boundary.index')=r.boundary.index
            if(do.save) {
                 #save(sub.t, file=paste0('/media/juno/bootstrap_hotspot/', total.hotspot.n, '.RData' ))
                 save(sub.t, file=paste0('/data/eQTL/RData/hotspots/', total.hotspot.n, '.RData2' ))
            }
        }
    }
    return(hotspots.refined)
}


















####################do bootstraps locally or on Juno nodes
bootstrap.Hotspots=function(hotspot.number.in ,onjuno=FALSE) {
    # this is already in the workspace
    #load('~/bootstrap_hotspot/gdata.s.by.chr.RData')
    #load('/media/juno/bootstrap_hotspot/gdata.s.by.chr.RData')
    if(!onjuno) {
        dir.create(paste0('/data/eQTL/RData/hotspots/',hotspot.number.in))
        load(paste0('/data/eQTL/RData/hotspots/', hotspot.number.in,'.RData'))
    } else{
        dir.create(paste0('/home/jbloom/bootstrap_hotspot/output/',hotspot.number.in))
        load(paste0('/home/jbloom/bootstrap_hotspot/input/', hotspot.number.in,'.RData'))

    }
    cc=attr(sub.t, 'chr')
    hs.n=attr(sub.t, 'hs.n')
    l.boundary.index=attr(sub.t, 'l.boundary.index')
    r.boundary.index=attr(sub.t, 'r.boundary.index')

    n.boots=1000
    peaks=rep(NA, n.boots)

    pb =txtProgressBar(min = 1, max = n.boots, style = 3)
    for(bootstrap.iteration in 1:n.boots) {
        setTxtProgressBar(pb, bootstrap.iteration)
        #print(bootstrap.iteration)
        sss=rep(NA, ncol(gdata.s.by.chr[[cc]]))
        g=gdata.s.by.chr[[cc]]
        boot.ind=sample(1:nrow(g), replace=T)
        p=sub.t[boot.ind,]
        g=g[boot.ind,]
        nn=ncol(p)
        SST=crossprod(p)
        t5=abs(eigen(SST)$values)
        pc.to.retain=which(cumsum(t5/sum(t5))<.999)
        #print(length(pc.to.retain))

        for(i in l.boundary.index:r.boundary.index ) { #1:length(sss)){
               #print(i)   
               test=lm.fit(as.matrix(g[,i]), p)
               t2=residuals(test)
               SSE=crossprod(t2)
               t4=abs(eigen(SSE)$values)
               X2=-(1012-1-((nn+2)/2))*(sum(log(t4[pc.to.retain]))-sum(log(t5[pc.to.retain])))
               sss[i]=X2
        } 
    peak=colnames(gdata.s.by.chr[[cc]])[which.max(sss)]
    #if(!onjuno) {
    #    save(peak, file=paste0('/data/eQTL/RData/hotspots/',hotspot.number.in, '/', bootstrap.iteration))
    #}
    #else {
    #    save(peak, file=paste0('/home/jbloom/bootstrap_hotspot/output/',hotspot.number.in, '/', bootstrap.iteration))
    #}
    peaks[bootstrap.iteration]=peak
    }
    close(pb)
    save(peaks, file=paste0('/data/eQTL/RData/hotspots/',hotspot.number.in, '.bootstrap.peaks'))
    return(peaks)
}
######## end  code for hotspot analysis --------------------------------------------------------------------------------------

downsampleMarkers=function(gdata.by.chr, gdata.s.by.chr) {
    gdata.downsampled=list()
    for(cc in names(gdata.by.chr))  {
      print(cc)
      g1=gdata.by.chr[[cc]]
      g2=gdata.s.by.chr[[cc]]
      XXr=crossprod(g2)/(nrow(g2)-1)
      fCX=sort(findCorrelation(XXr, cutoff=.99, exact=F))
      gdata.downsampled[[cc]]=g1[,-fCX]  
    }
    return(gdata.downsampled)
}


removeAdditiveEffects =function(t.tpm.matrix, background.QTL, covariates, gdata) { 
         plist=lapply(colnames(t.tpm.matrix), function(i) { 
              if(!(i %in% names(background.QTL)) ) { return(  residuals(lm(t.tpm.matrix[,i]~covariates)) ) }
              else{  residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,background.QTL[[i]]$gcind]) ) }
         })
         preal=do.call('cbind', plist)
         colnames(preal)=colnames(t.tpm.matrix)

         presid=preal
         A=A.mat(gdata)/2
                 
         eigA=doEigenA_forMM(presid,A)
         svdA=svd(A)
         pb=txtProgressBar(min=1, max=ncol(presid), style=3)
         for(tp in 1:ncol(presid)){
             setTxtProgressBar(pb,tp)
             rr=m.S(presid[,tp], K=A,  theta=eigA$theta, Q=eigA$Q)
             #W=solve(rr[1]*A+rr[2]*diag(nrow(presid)))
             W=svdA$u %*% tcrossprod(diag(1/((svdA$d*rr[1])+(rr[2]))), svdA$v)
             if(rr[1]>0) {
                 blups=calc.BLUPS(rr[1]*A,diag(nrow(presid)),W,presid[,tp],matrix(1,nrow(presid),1),0 )[,1]
                 presid[,tp]=as.vector(presid[,tp] - blups)
             }
         }
         rm(W)
         close(pb)
         return(presid)
        # presid now has additive effects removed 
}



do2locusScan=function(pheno.scaled, gdata.downsampled, out.dir, marker.gap=10){
    n.pheno=nrow(pheno.scaled) 
    for(cc1 in 1:16) {
        for(cc2 in cc1:16) {
            print(cc1)
            print(cc2)
            
            g1=gdata.downsampled[[unique.chrs[cc1]]]
            g2=gdata.downsampled[[unique.chrs[cc2]]]
            
            # build vector of positions to test
            if(cc1==cc2) {
                bM=(as.matrix(bandSparse(ncol(g1),ncol(g1),c(0:marker.gap), symmetric=T)))
                bM[lower.tri(bM)]=1
                to.test=which(bM==0, arr.ind=T)
            }else {
                 bM=matrix(0, ncol(g1), ncol(g2))
                 to.test=which(bM==0, arr.ind=T)
            }
            gint=scale(g1[,to.test[,1]]*g2[,to.test[,2]])

            s2=fasterLOD(n.pheno, pheno.scaled, gint)
            a2=array(0,dim=c(ncol(pheno.scaled),ncol(g1), ncol(g2)))
            dimnames(a2)=list(colnames(pheno.scaled), colnames(g1), colnames(g2))
            attr(a2, '1')=unique.chrs[cc1]
            attr(a2, '2')=unique.chrs[cc2]
            for(k in 1:nrow(to.test)){   a2[,to.test[k,1], to.test[k,2]]=s2[,k]   }
            save(a2,file=paste0(out.dir, cc1, '_', cc2))
            #image.plot(a2[33,,])     
        }
    }
}


peakfinder2D = function(i2d, threshold,pdist.m=10) {
    #   i2d   = iLODsh[trait,,]
        ipeaks  = which(i2d>threshold, arr.ind=T)
        ipeaks = cbind(ipeaks, i2d[ipeaks])
        ipeaks = cbind(ipeaks, rep(NA, nrow(ipeaks)))
        ipeaks = cbind(ipeaks, ipeaks[,3])
        #colnames(ipeaks) = c('x', 'y', 'lod', 'group', 'glod')
        ipeaks.m1.name=dimnames(i2d)[[1]][ipeaks[,1]]
        ipeaks.m2.name=dimnames(i2d)[[2]][ipeaks[,2]]
        ipeaks=data.frame(ipeaks,  ipeaks.m1.name, ipeaks.m2.name, stringsAsFactors=F)
        colnames(ipeaks) = c('x', 'y', 'lod', 'group', 'glod', 'm1', 'm2')
        
        g=1
        while(sum(is.na(ipeaks[,'group']))>0) {
            peak = which.max(ipeaks[,'glod'])
            pdist = sqrt(abs(ipeaks[peak,'x']-ipeaks[,'x'])^2 + abs(ipeaks[peak,'y']-ipeaks[,'y'])^2)
            gpeaks = which(pdist< pdist.m & is.na(ipeaks[,'group']))
            if(length(gpeaks)>2) {
               ipeaks[gpeaks, 'group']=g
               ipeaks[gpeaks, 'glod']=0
               g=g+1
               #print(g)
            } else{
               ipeaks[gpeaks, 'group']=0
               ipeaks[gpeaks, 'glod']=0
            }
        }
        ipeaks=ipeaks[ipeaks[,'group']>0,]
        ips = split(data.frame(ipeaks), ipeaks[,'group'])
        ipsp=do.call('rbind', lapply(ips, function(x) {
                lmax= which.max(x$lod)
                x[lmax, c(1,2,3,4,6,7)]
         }))
        return(ipsp)
}

#aggregatePeaks= function(in.dir, out.file) {
#    peaks.2D.real=list()
#    for(cc1 in 1:16) {
#        print(cc1)
#        for(cc2 in cc1:16) {
#            print(cc2)
#            load(paste0(in.dir, cc1, '_', cc2))
#             peaks.2D.real[[paste0(cc1, '_', cc2) ]]=peakfinder3D(a2, 3, pdist.m=50)
#        }
#    }
#    saveRDS(peaks.2D.real, file=out.file) 
#}

# calculate average number of peaks from permutations and divide by observed number of peaks from un-permuted data

peakfinder3D = function(i2d, threshold,pdist.m=10) {
    #   i2d   = iLODsh[trait,,]
        ipeaks  = which(i2d>threshold, arr.ind=T)
        ipeaks = cbind(ipeaks, i2d[ipeaks])
        ipeaks = cbind(ipeaks, rep(NA, nrow(ipeaks)))
        ipeaks = cbind(ipeaks, ipeaks[,4])
        #colnames(ipeaks) = c('x', 'y', 'lod', 'group', 'glod')
        ipeaks.m1.name=dimnames(i2d)[[2]][ipeaks[,2]]
        ipeaks.m2.name=dimnames(i2d)[[3]][ipeaks[,3]]
        ipeaks=data.frame(ipeaks,  ipeaks.m1.name, ipeaks.m2.name, stringsAsFactors=F)
        colnames(ipeaks) = c('trait','x', 'y', 'lod', 'group', 'glod', 'm1', 'm2')
        ipeaks$trait=dimnames(i2d)[[1]][ipeaks$trait]

        peaks=list()
        ip=split(ipeaks, ipeaks$trait)
        for(tt in names(ip) ) {
            ipeaks=ip[[tt]] 
            g=1
             while(sum(is.na(ipeaks[,'group']))>0) {
                peak = which.max(ipeaks[,'glod'])
                pdist = sqrt(abs(ipeaks[peak,'x']-ipeaks[,'x'])^2 + abs(ipeaks[peak,'y']-ipeaks[,'y'])^2)
                gpeaks = which(pdist< pdist.m & is.na(ipeaks[,'group']))
                    if(length(gpeaks)>2) {
                       ipeaks[gpeaks, 'group']=g
                       ipeaks[gpeaks, 'glod']=0
                       g=g+1
                       #print(g)
                } else{
                       ipeaks[gpeaks, 'group']=0
                       ipeaks[gpeaks, 'glod']=0
                 }
              }
            ipeaks=ipeaks[ipeaks[,'group']>0,]
            ips = split(data.frame(ipeaks), ipeaks[,'group'])
            ipsp=do.call('rbind', lapply(ips, function(x) {
                lmax= which.max(x$lod)
                x[lmax, c(1:5, 7:8)]
            }))
            
            peaks[[tt]]=ipsp
        }
        return(peaks)
}

collapse2D=function(obs2D) {
    pdr=lapply(obs2D, function(x) do.call('rbind', x))
    pdr2=do.call('rbind',  pdr)
    return(pdr2)
}


#new function to aggregate 2D peaks without saving the raw matrix of all 2D lod scores (2D peak detection on the fly)
do2locusScan.reduced_output=function(pheno.scaled, gdata.downsampled, all.peaks.DS, marker.gap=20){
    n.pheno=nrow(pheno.scaled) 
    
    full.scan.peaks=list()
    marginal.scan.peaks=list()

    # could parallelize this
    #cc=which(lower.tri(matrix(1, 16, 16), diag=T), arr.ind=T)
    #cc.table=cbind(cc[,2], cc[,1])

    for(cc1 in 1:16) {
        for(cc2 in cc1:16) {
            print(cc1)
            print(cc2)
            
            g1=gdata.downsampled[[unique.chrs[cc1]]]
            g2=gdata.downsampled[[unique.chrs[cc2]]]
            
            # build vector of positions to test
            if(cc1==cc2) {
                bM=(as.matrix(bandSparse(ncol(g1),ncol(g1),c(0:marker.gap), symmetric=T)))
                bM[lower.tri(bM)]=1
                to.test=which(bM==0, arr.ind=T)
            }else {
                 bM=matrix(0, ncol(g1), ncol(g2))
                 to.test=which(bM==0, arr.ind=T)
            }
            gint=scale(g1[,to.test[,1]]*g2[,to.test[,2]])

            s2=fasterLOD(n.pheno, pheno.scaled, gint)
            a2=array(0,dim=c(ncol(pheno.scaled),ncol(g1), ncol(g2)))
            dimnames(a2)=list(colnames(pheno.scaled), colnames(g1), colnames(g2))
            attr(a2, '1')=unique.chrs[cc1]
            attr(a2, '2')=unique.chrs[cc2]
            for(k in 1:nrow(to.test)){   a2[,to.test[k,1], to.test[k,2]]=s2[,k]   }
            
            # the time to save and reload these arrays is pretty awful, switch to peak detection on the fly
            #save(a2,file=paste0(out.dir, cc1, '_', cc2))
            #image.plot(a2[33,,])     
        
        
            # full 2D matrix is never saved (or loaded)
            # get peak positions for full scan
            full.scan.peaks[[paste0(cc1, '_', cc2) ]] = peakfinder3D(a2, 3, pdist.m=50)

            # cc1 and cc2 can be found inside do2locusScan() 
            # for marginal scan need all.peaks.DS (full scan peaks snapped to downsampled markers) 
            peaks_on_this_chr_pair=all.peaks.DS[all.peaks.DS$chr %in% c(unique.chrs[cc1], unique.chrs[cc2]), ]
            pocp=split(peaks_on_this_chr_pair, peaks_on_this_chr_pair$gene)
            tind=dimnames(a2)[[1]][sort(match(names(pocp), dimnames(a2)[[1]]))]
            pocp=pocp[tind]
            a3=a2[tind,,]
            rm(a2)
            a6=array(0, dim(a3))
            dimnames(a6)=dimnames(a3)
            for(tt in tind) {
                   a4=a3[tt,,]
                   p.relevant=pocp[[tt]]
                   m1=match(p.relevant$pmarker.ds, dimnames(a4)[[1]])
                   m2=match(p.relevant$pmarker.ds, dimnames(a4)[[2]])
                   na.omit(m1)
                   na.omit(m2)
                   a5=matrix(0, dim(a4)[1],dim(a4)[2])
                   a5[na.omit(m1),]=a4[na.omit(m1),]
                   a5[,na.omit(m2)]=a4[,na.omit(m2)]
                   a6[tt,,]=a5
               }
               rm(a3); rm(a4); rm(a5)
               marginal.scan.peaks[[paste0(cc1, '_', cc2) ]] = peakfinder3D(a6, 3, pdist.m=50)
        
        }
    }
    return(list(full.scan.peaks=full.scan.peaks, marginal.scan.peaks=marginal.scan.peaks))
}






















#



aggregatePeaks= function(in.dir, out.file1, out.file2, all.peaks.DS) {
    unique.chrs=paste0('chr', as.roman(1:16))
    peaks.2D.real=list()
    peaks.2D.marginal=list()
    for(cc1 in 1:16) {
        print(cc1)
        for(cc2 in cc1:16) {
                print(cc2)
                # loads a2
                load(paste0(in.dir, cc1, '_', cc2))
                peaks.2D.real[[paste0(cc1, '_', cc2) ]]    =peakfinder3D(a2, 3, pdist.m=50)
                peaks_on_this_chr_pair=all.peaks.DS[all.peaks.DS$chr %in% c(unique.chrs[cc1], unique.chrs[cc2]), ]
                pocp=split(peaks_on_this_chr_pair, peaks_on_this_chr_pair$gene)
                tind=dimnames(a2)[[1]][sort(match(names(pocp), dimnames(a2)[[1]]))]
                pocp=pocp[tind]
                a3=a2[tind,,]
                rm(a2)
                a6=array(0, dim(a3))
                dimnames(a6)=dimnames(a3)
                for(tt in tind) {
                    a4=a3[tt,,]
                    p.relevant=pocp[[tt]]
                    m1=match(p.relevant$pmarker.ds, dimnames(a4)[[1]])
                    m2=match(p.relevant$pmarker.ds, dimnames(a4)[[2]])
                    na.omit(m1)
                    na.omit(m2)
                    a5=matrix(0, dim(a4)[1],dim(a4)[2])
                    a5[na.omit(m1),]=a4[na.omit(m1),]
                    a5[,na.omit(m2)]=a4[,na.omit(m2)]
                    a6[tt,,]=a5
                    #print(tt)
                }
                rm(a3)
                rm(a4)
                rm(a5)
                peaks.2D.marginal[[paste0(cc1, '_', cc2) ]]=peakfinder3D(a6, 3, pdist.m=50)
             gc()
        }
    }
    saveRDS(peaks.2D.real, file=out.file1) 
    saveRDS(peaks.2D.marginal, file=out.file2) 
}

## now given refined hotspots go back and 
#local.mediation=list()
#hotspot.peaks=do.call('c', htabler[,'peak']) #as.character(htabler$peak)
#local.mediation=foreach(h=1:length(hotspot.peaks)) %dopar%  {  #in hotspot.peaks ){
#   hp=hotspot.peaks[h]
#   print(hp)
#   mint=marker.GR[marker.GR$mname==hp]
#   # all genes within 10 kb of peak marker
#   genes.near.hotspot=subsetByOverlaps(gene.GR, mint+10000)
#
#   null= lm(t.tpm.matrix~covariates.OD)
#   full= lm(t.tpm.matrix~covariates.OD+gdata[,hp])
#   rn=residuals(null)
#   fn=residuals(full) 
#   rss1=colSums(rn^2)
#   rss2=colSums(fn^2)
#   Fstat=(rss1-rss2)/(rss2/(1012-14))
#   hotspot.model=pf(Fstat, 1,1012-14, lower.tail=F)
#   
#   # anova(lm(t.tpm.matrix[,i]~covariates) ,lm(t.tpm.matrix[,i]~covariates+gdata[,hp]) )
#   #  hotspot.model = rep(NA, 6290) 
#   # for(i in 1:6290) {
#   #     hotspot.model[i]=anova(lm(t.tpm.matrix[,i]~covariates) ,lm(t.tpm.matrix[,i]~covariates+gdata[,hp]) )$'Pr(>F)'[2]
#   # }
#   qhot=qvalue(hotspot.model)
#   sig.at.hotspot=which(qhot$qvalue<.05)
#    
#    print(length(genes.near.hotspot))
#    int.list=list()
#    for(cg in 1:length(genes.near.hotspot)) {
#        print(cg)
#        cisg=genes.near.hotspot[cg]$ORF
#        
#        sya=scale(t.tpm.matrix[,sig.at.hotspot])
#        rmcis.ind=match(cisg, colnames(sya))
#        if(!is.na(rmcis.ind)) {
#            sy=sya[,-rmcis.ind]
#        } else {sy=sya} 
#
#        rmcis.ind=match(cisg, colnames(t.tpm.matrix))
#        cis.e=scale(t.tpm.matrix[,rmcis.ind])
#        #testing if growth (OD)  mediates QTL effects on expression levels  -----------------------
#        #sy=scale(t.tpm.matrix[,gind])
#        #sOD=scale(OD.cov)
#        
#         M1= lm(sy  ~ covariates.OD +gdata[,hp]-1)
#          # check that M2 is actually significant (source of NAs?)  
#         M2= lm(cis.e ~ covariates.OD+gdata[,hp]-1)
#         M3= lm(sy  ~ covariates.OD+cis.e+gdata[,hp]-1)
#         
#         # iterate over all the trans genes or compare to best M4 model
#         # is reduction in M3 model greater than reduction in M4 model
#         M4=list()
#         for(k in colnames(sy)){
#         M4[[k]]= lm(cis.e  ~ covariates.OD+sy[,k]+gdata[,hp]-1)
#         }
#         sM4=lapply(lapply(M4, summary), coefficients)
#         # compare M4,M2 and M1 for reverse model
#
#         sM1= coefficients(summary(M1)) #$coefficients
#         sM2= coefficients(summary(M2)) #$coefficients
#         sM3= coefficients(summary(M3)) #$coefficients
#         
#         #https://en.wikipedia.org/wiki/Sobel_test-----------------------------------------------------------------------------------------------
#         # cis expression mediating trans expression
#         tau=sapply(sM1, function(x)    x[grepl('gdata', rownames(x)),'Estimate'])      #sM1[which(grepl('gdata', rownames(sM1))),'Estimate']
#         tau.se=sapply(sM1, function(x) x[grepl('gdata', rownames(x)),'Std. Error']) #sM1[which(grepl('gdata', rownames(sM1))),'Std. Error']
#
#         a=sM2[which(grepl('gdata', rownames(sM2))),'Estimate']
#         a.se=sM2[which(grepl('gdata', rownames(sM2))),'Std. Error']
#
#         b=sapply(sM3, function(x)    x[grepl('cis.e', rownames(x)),'Estimate'])   # sM3[which(grepl('cis.e', rownames(sM3))),'Estimate']
#         b.se=sapply(sM3, function(x)    x[grepl('cis.e', rownames(x)),'Std. Error']) 
#         #b.se=sM3[which(grepl('cis.e', rownames(sM3))),'Std. Error']
#         
#         tau.prime=sapply(sM3, function(x)    x[grepl('gdata', rownames(x)),'Estimate']) #sM3[which(grepl('gdata', rownames(sM3))),'Estimate']
#         tau.prime.se=sapply(sM3, function(x)    x[grepl('gdata', rownames(x)),'Std. Error'])
#         #sM3[which(grepl('gdata', rownames(sM3))),'Std. Error']
#        
#         pooled.se = a^2 * b.se^2 + b^2*a.se^2 + a.se^2+b.se^2
#         #Zscore=(tau-tau.prime)/pooled.se
#         Zscore=(a*b)/pooled.se
#         # double check logic here JB 6/22
#         Zscore[abs(tau.prime)>abs(tau)]=NA
#         p.med=pnorm(abs(Zscore), lower.tail=F)*2
#         dff=data.frame(tau, tau.se, tau.prime, tau.prime.se, a,a.se, b,b.se,Zscore, p.med)
#         #hist(dff$p.med)
#         rownames(dff)=colnames(sy)
#         #---------------------------------------------------------------------------------------------------------------------------------------
#
#
#         # compare M4,M2 and M1 for reverse model
#         # trans expression mediating cis expression 
#         r.tau=sM2[which(grepl('gdata', rownames(sM2))),'Estimate']
#         r.tau.se=sM2[which(grepl('gdata', rownames(sM2))),'Std. Error']
#
#         r.a=sapply(sM1, function(x)    x[grepl('gdata', rownames(x)),'Estimate'])      #sM1[which(grepl('gdata', rownames(sM1))),'Estimate']
#         r.a.se=sapply(sM1, function(x) x[grepl('gdata', rownames(x)),'Std. Error']) #sM1[which(grepl('gdata', rownames(sM1))),'Std. Error']
#        
#         r.b=sapply(sM4, function(x)    x[grepl('sy', rownames(x)),'Estimate'])
#         r.b.se=sapply(sM4, function(x)    x[grepl('sy', rownames(x)),'Std. Error']) 
#        
#         r.tau.prime=sapply(sM4, function(x)    x[grepl('gdata', rownames(x)),'Estimate']) #sM3[which(grepl('gdata', rownames(sM3))),'Estimate']
#         r.tau.prime.se=sapply(sM4, function(x)    x[grepl('gdata', rownames(x)),'Std. Error'])
#
#         r.pooled.se = r.a^2 * r.b.se^2 + r.b^2*a.se^2 + r.a.se^2+r.b.se^2
#         #Zscore=(tau-tau.prime)/pooled.se
#         r.Zscore=(r.a*r.b)/r.pooled.se
#         # double check logic here JB 6/22
#         r.Zscore[abs(r.tau.prime)>abs(r.tau)]=NA
#         r.p.med=pnorm(abs(r.Zscore), lower.tail=F)*2
#         r.dff=data.frame(r.tau, r.tau.se, r.tau.prime, r.tau.prime.se, r.a, r.a.se, r.b,r.b.se, r.Zscore, r.p.med)
#         #hist(dff$p.med)
#         rownames(r.dff)=colnames(sy)
#
#         d=cbind(dff,r.dff)
#         #----------------------------------------------------------------------------------------------
#         #local.mediation[[hp]][[cisg]]=d
#         int.list[[cisg]]=d
#    }
#   return(int.list)   
#}
#names(local.mediation)=hotspot.peaks
#save(local.mediation, file='~/Dropbox/Public/eQTL/local.mediation.RData')
#llt=lapply(local.mediation, function(x) 
#       {
#          sapply(x, function(y) sum(y[,'p.med']<1e-20 & ((y[,'r.p.med']> y[,'p.med']) | is.na(y[,'r.p.med']) ) , na.rm=T)/length(y$p.med))
#       }
#)


















#peakfinder 3d marginal
#cc1=4
#cc2=14
#
#   if(hs.n==1) {
#       g=gdata.s.by.chr[[cc]]
#       XXr=crossprod(g)/(ncol(tpm.matrix)-1)
#       #XYr=crossprod(g,tmm)/(ncol(tpm.matrix)-1)
#       #scanone.all=(-1012*log(1-XYr^2))/(2*log(10))
#       fCX=sort(findCorrelation(XXr, cutoff=.99, exact=F))
#       g.s=g[,-fCX]
#       #sXXr=crossprod(g.s)/(ncol(tpm.matrix)-1)
#       #sXYr=crossprod(g.s,tmm)/(ncol(tpm.matrix)-1)
#       #scanone.red=(-1012*log(1-sXYr^2))/(2*log(10))
#       #cor(apply(scanone.all,2,max), apply(scanone.red, 2, max))
#       #plot(apply(scanone.all,2,max), apply(scanone.red, 2, max), xlim=c(0,40), ylim=c(0,40) )
#       for( i in 1:ncol(tmm)) {
#          # print(i)
#           # with preset hyperparameters
#          # optimize this 
#           bl[[cc]][[colnames(tmm)[i]]]=EBglmnet(g.s, tmm[,i], hyperparameters=c(.1,.1),verbose=0)
#           # choose hyperparameters
#           #bl[[cc]][[colnames(tmm)[i]]]=cv.EBglmnet(g.s, tmm[,i],verbose=-1)
#       }
#
#       gind=match(names(bl[[cc]]), gene.annot.df[,1])
#       q.bl.c=lapply(bl[[cc]], function(x) nrow(x$fit))
#       q.bl.c=do.call('c', q.bl.c)
#       greps=rep(gind,as.vector(q.bl.c))
#       bl.loc=sapply(bl[[cc]], function(x) x$fit[,1] )
#       bl.loc=sapply(bl[[cc]], function(x) x$fit[,1] )
#       bl.beta=sapply(bl[[cc]], function(x) x$fit[,'beta'] )
#       bl.var=sapply(bl[[cc]], function(x) x$fit[,'posterior variance'] )
#       bl.t=sapply(bl[[cc]], function(x) x$fit[,'t-value'] )
#       bl.p=sapply(bl[[cc]], function(x) x$fit[,'p-value'] )
#       #plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.beta)))
#       #plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.var)))
#       #plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.t)))
#       #plot(unlist(as.vector(bl.loc)), -log10(unlist(as.vector(bl.p))))
# }











# bb=seq(0,ncol(tmp),30)
#   if(max(ncol(tmp))>max(bb)) { bb=c(bb, max(bb)+30) }
#   l.in.bins=cut(LOD.kt.wmax, bb)
#   r.l.in.bins=rle(sort(as.vector(l.in.bins)))
#   max.bin=r.l.in.bins$values[which.max(r.l.in.bins$lengths)]
#
# mb=gsub('\\(', '', max.bin)
#   mb=gsub('\\]', '', mb)
#   mb=as.numeric(strsplit(mb, ',')[[1]])
#
#   to.scan=names(LOD.kt.max)[which(l.in.bins %in% max.bin)]
#
#  #regress out effects of background QTL on chromosome
#   #for(ttt in colnames(sub.t)) {
#   #    #print(ttt)
#   #     bcl=bl[[cc]][[ttt]]$fit[,'locus1']
#   #     bq=bcl[g.s.pos[bcl]<l.boundary | g.s.pos[bcl]>r.boundary]
#   #     if(length(bq)>0) {
#   #         print(length(bq))
#   #         sub.t[,ttt]=scale(residuals(lm(sub.t[,ttt]~g.s[,bq])))[,1] 
#   #     }
#   #}
#
#  #max.to.scan=LOD.kt.max[to.scan]
#   #max.to.scan=(sort(max.to.scan, decreasing=T))
#   #max.to.scan=max.to.scan[max.to.scan>3.5]
#  #nn=min(800, length(keep.transcripts))
#   #if(cc=='chrXIV' & nn>799) {nn=500}
#   #if(cc=='chrXIV' & nn>799) {nn=500}
#
#   #sub.t=tmm[,names(sort(max.obsLOD,decreasing=T)[1:nn])]
#
#  if(mb[2]>ncol(g)) {mb[2]=ncol(g) }
#  # this can't be right 
#  if(mb[1]>1) {mb[1]=1 }
#   pc.to.retain=which(cumsum(t5$d/sum(t5$d))<.999)
#
#   print(length(pc.to.retain))
#   #  http://arxiv.org/abs/1305.5870
#   for(i in 1:length(sss)){
#       setTxtProgressBar(pb, i)
#       test=lm(sub.t~gdata.s.by.chr[[cc]][,i])
#       t2=residuals(test)
#       SSE=crossprod(t2)
#       t4=svd(SSE)
#       X2=-(1012-1-((nn+2)/2))*(sum(log(t4$d[pc.to.retain]))-sum(log(t5$d[pc.to.retain])))
#       sss[i]=X2
#   }
#   close(pb)
   #plot(-log10(pchisq(sss,nn,lower.tail=F)))
   #sss=-log10(pchisq(sss,nn,lower.tail=F))
   # plot(sss)










# Obsolete code ::

# load express output -------------------------------------------------------------------------------------------
#    transcript.fa.file=paste0(base.dir, 'reference/R64-1.1.cdna.ncrna.fa')
#    transcripts.fa=read.fasta(transcript.fa.file)
#    transcript.info=sapply(transcripts.fa, function(x) attr(x,'Annot'))

#dir.in=paste0(base.dir, 'star_out/')
#counts.list=list()
#for(sub.dir in list.files(dir.in)) {
#    tdf=read.delim(paste0(dir.in, sub.dir, '/F/results.xprs'), header=T, sep='\t')
#    print(paste(sub.dir, nrow(tdf)))
#    counts.list[[sub.dir]]=tdf[order(as.character(tdf$target_id)),]
#}
#raw.count.matrix = do.call('cbind', lapply(counts.list, function(x) (x$est_counts)))
#rownames(raw.count.matrix)=as.character(counts.list[[1]]$target_id)
#eff.len.matrix = do.call('cbind', lapply(counts.list, function(x) (x$eff_length)))
#rownames(eff.len.matrix)=as.character(counts.list[[1]]$target_id)

#plot(countToTpm(counts.list[[1]]$est_counts, counts.list[[1]]$eff_length), counts.list[[1]]$tpm)
#--------------------------------------------------------------------------------------------------------------







# Code to extract the worst 172 samples for repeat growth, extraction, and library prep

# we can extract the computed tpm values but easier to recompute them from counts later given we're going to shuffle genes around
# and subset the data in various ways
#tpm.matrix = do.call('cbind', lapply(counts.list, function(x) (x$tpm)))
#rownames(tpm.matrix)=as.character(counts.list[[1]]$target_id)
# ---------------------------------------------------------------------------------------------------------------
##repeat the worst 172 -----------------------------------------------------------------------------------------
#rep.172=sort(sample.sum)[1:172]
###name,ynb.growth,ynb.rank,old.plate,old.well,new.plate,new.well,vol
#rep.172.df=data.frame(do.call('rbind', strsplit(names(rep.172), '-'))[,1:5], stringsAsFactors=F)
#names(rep.172.df)=c('name', 'old.plate', 'old.well', 'int.plate', 'int.well')
#rep.172.df[order(rep.172.df$int.plate),]
#cols=matrix(rep(1:12,each=8),8,12)
#rows=t(matrix(rep(toupper(letters[1:8]),each=12), 12,8))
#wells=matrix(paste0(rows,cols), 8, 12)
#eQTL.12=rep.172.df[order(rep.172.df$int.plate),][1:86,]
#eQTL.13=rep.172.df[order(rep.172.df$int.plate),][87:172,]
#eQTL.12$new.plate='BYxRM_eQTL_12'
#eQTL.13$new.plate='BYxRM_eQTL_13'
#nw12=sample(wells,86)
#inc.12=sample(wells[(!(wells %in% nw12))])
#eQTL.12$new.well=nw12
#nw13=sample(wells,86)
#inc.13=sample(wells[(!(wells %in% nw13))])
#eQTL.13$new.well=nw13
#padd=data.frame(c(rep('1879', 5), rep('1950', 5)),       rep('BY_RM',10),       c(wells[1:5,1], wells[1:5,3]),      rep('-', 10),      rep('-', 10),      rep('BYxRM_eQTL_12', 10),      c(inc.12),stringsAsFactors=F)
#    names(padd)=c('name', 'old.plate', 'old.well', 'int.plate', 'int.well', 'new.plate', 'new.well')
#eQTL.12=rbind(eQTL.12, padd)
#padd=data.frame(c(rep('1879', 5), rep('1950', 5)), rep('BY_RM',10), c(wells[6:8,1], wells[1:2,2],  wells[6:8,3], wells[1:2,4] ),  rep('-', 10),      rep('-', 10),      rep('BYxRM_eQTL_13', 10),      c(inc.13),stringsAsFactors=F)
#    names(padd)=c('name', 'old.plate', 'old.well', 'int.plate', 'int.well', 'new.plate', 'new.well')
#eQTL.13=rbind(eQTL.13, padd)
#eQTL.12$vol=5
#eQTL.13$vol=5
#write.table(eQTL.12, file=paste0('/data/eQTL/rearrange/', 12, '.csv'),  row.names=F, sep=',',quote=F)
#write.table(eQTL.13, file=paste0('/data/eQTL/rearrange/', 13, '.csv'),  row.names=F, sep=',',quote=F)
#---------------------------------------------------------------------------------------------------------------------

#New, do 2-locus scan for every transcript with significant signal
#######################
#    g=gdata.s.by.chr[[cc]]
#    p=presid[, (names(peak.gene.ind))]
    
#    nmarkers.chr= (ncol(g)) 
#    stat.mat.nm=matrix(0, ncol(p2), nmarkers.chr)
#    stat.mat.am=array(0, c(ncol(p2), nmarkers.chr, nmarkers.chr))
#    stat.mat.fm=array(0, c(ncol(p2), nmarkers.chr, nmarkers.chr))
#    system.time({
#    for(j in 1:1) {# (nmarkers.chr-1) ) {
#        print(j)
#        nm=colSums(lm.fit(as.matrix(g[,j]), p2)$residuals^2)
#        stat.mat.nm[,j]=nm
#       for(k in (j+1):nmarkers.chr ){
#            #    print(k)
#        g2=g[,c(j, k)]
#        g3=cbind(g[,c(j, k)], g[,j]*g[,k])
#        am=colSums(lm.fit(g2, p2)$residuals^2)
#        fm=colSums(lm.fit(g3, p2)$residuals^2)
#        stat.mat.am[,j,k]=am
#        stat.mat.fm[,j,k]=fm
#        }
#        #1012*log(nm/am)/4.605
#    }
#    })
#######################

#.sig[,ilibrary(RcppArmadillo)
#library(Rcpp)
#library(inline)]
#
##http://math.stackexchange.com/questions/182309/block-inverse-of-symmetric-matrices
#
#te=XXr[c(1,100,600), c(1,100,600)]
#system.time( {replicate(10000, {invte=solve(te)}) } )
#
#tA=te[c(1,2), c(1,2)]
#invtA=solve(tA)
#
#system.time( 
#            {replicate(10000, {
#delta=te[c(1,2),3]
#idm= invtA %*% (delta)
#
#u=(te[3,3]-t(delta) %*% idm)[1,1]
#e1.3= -idm/u
#invte.b=invtA+((idm %*% t(delta) %*% invtA)/u)
#}
#
#)})
#Xinvs=list()
#Xn=list()
#for(i in 1:(794-20)) {
#    print(i)
#    for(j in (i+20):794) {
#    Xinvs[[paste(i,j,sep='_')]]=XXr[c(i,j),c(i,j)]
#    Xn[[paste(i,j,sep='_')]]=cbind(XYr[c(i,j)], c(0,0))
#    }
#}
#
#
#XB=bdiag(Xinvs)
#XN=bdiag(Xn)
#
#
#
#
#
#Xinvs=list()
#for(i in 1:(794-20-20-20)) {
#    print(i)
#    for(j in (i+20):(794-20-20)) {
#        for(k in (j+20):(794-20)) {
#              for(l in (k+20):794) {
#                    (t(XYr[c(i,j,k,l)])%*%solve(XXr[c(i,j,k,l),c(i,j,k,l)])%*% XYr[c(i,j,k,l)])[1,1]
#        }
#        }
#    }
#}
#
#
#
#
#
#             (t(XYr[c(i,j)])%*%solve(XXr[c(i,j),c(i,j)])%*% XYr[c(i,j)])[1,1]
#
#i=1
#j=50
#
#xn=cbind(XYr[c(i,j)],c(0,0))
#t(xn)%*%solve(XXr[c(i,j),c(i,j)])%*%xn
#
#i=20
#j=500
#   
#xn2=cbind(XYr[c(i,j)],c(0,0))
#t(xn2)%*%solve(XXr[c(i,j),c(i,j)])%*%xn2
#
#xn3=cbind(xn,xn2)
#bd=bdiag((XXr[c(1,50),c(1,50)]), (XXr[c(20,500),c(20,500)]) )
#xd=bdiag(xn,xn2)
#
#t(xd) %*%solve(bd) %*%(xd)
#
#xn3%*%bd%*%t(xn3)
#
#t(xn)%*%solve(XXr[c(i,j),c(i,j)])%*%xn
#
#xn2=
#
#
#
#
#
#allXXr2=
#
#
#save.image('~/Desktop/011116_eQTL.RData')
#library(EBglmnet)
#library(MASS)
#library(lassoshooting)
## cut 
##cut.all=split(1:5618, cut(1:5618,11))
##system.time({ test2.all=mcMap(function(i){  mr2.scantwo.additive(XYr[,i], XXr,n.marker.space) }, cut.all, mc.cores=11) })
##test2.all=do.call('rbind', test2.all)
#
#n.marker.space=20
#sy=scale(s)
#stpreds=scale(gforsim)
#
#image.plot(XXr)
#dXXr=as.dist(1-XXr^2)
#fit <- hclust(dXXr, method="ward.D")
##plot(fit) # display dendogram
#cg=cutree(fit, k=ceiling(ncol(stpreds)/20))
#abline(v=cumsum(rle(cg)$lengths)/794)
#
#
#
#A.chr.list=list()
# for(n in 1:max(cg)) {
#        print(n)
#        cluster.inds = which(cg==n)
#        A.chr.list[[as.character(n)]]=tcrossprod(stpreds[,cluster.inds])/length(cluster.inds)
#}
#
#fitall=regress(sy[,1]~1, ~A.chr.list[[1]]+A.chr.list[[2]]+A.chr.list[[3]]+A.chr.list[[4]]+A.chr.list[[5]]+
#A.chr.list[[6]]+A.chr.list[[7]]+A.chr.list[[8]]+A.chr.list[[9]]+A.chr.list[[10]]+
#A.chr.list[[11]]+A.chr.list[[12]]+A.chr.list[[13]]+A.chr.list[[14]]+A.chr.list[[15]]+
#A.chr.list[[16]]+A.chr.list[[17]]+A.chr.list[[18]]+A.chr.list[[19]]+A.chr.list[[20]], verbose=T, pos=rep(T,21))
#
#
#m3.test=mr3.scantwo.BI.additive.2(XYr[,1:2], XXr, 20)
#
#
#fm.in=paste0('~', paste(paste0('A.chr.list[[', 1:40, ']]'), collapse='+'))
#
#fitall=regress(sy[,2]~1, fm.in, verbose=T, pos=rep(T,41))
#[2,]  280  383  481 87.14
#
#barplot(fitall$sigma[-41], names.arg=1:40)
#
#
#
#
##library(pls)
#p1=pcr(sy[,1]~stpreds-1, ncomp=50)
#
#
#
#
##A.chr=tcrossprod(stpreds)/ncol(stpreds)
##Agroup=which(cg==7 | cg==11 | cg==12)
##Ao=tcrossprod(stpreds[,Agroup])/length(Agroup)
#
## estimate VC for simulated data
#library(leaps)
#library(caret)
#
#A=tcrossprod(stpreds)/ncol(stpreds)
#eigA=doEigenA_forMM(presid,A)
#vcEstMap=matrix(0, ncol(presid),2)
## calculate mixed model, one term for additive variance  -------------------------------------------
#n.cores=11
#for(i in 1:ncol(presid) ){
#    print(i)
#   vcEstMap[i,]= m.S(presid[,i], K=A,  theta=eigA$theta, Q=eigA$Q)
#}
#
#XXr=crossprod(stpreds)/(n-1)
#
#fCX=sort(findCorrelation(XXr, cutoff=.99, exact=F))
#stsub=stpreds[,-fCX]
#As=tcrossprod(stsub)/ncol(stsub)
#eigAs=doEigenA_forMM(presid,As)
#vcEstMapS=matrix(0, ncol(presid),2)
## calculate mixed model, one term for additive variance  -------------------------------------------
#n.cores=11
#for(i in 1:ncol(presid) ){
#    print(i)
#   vcEstMapS[i,]= m.S(presid[,i], K=As,  theta=eigAs$theta, Q=eigAs$Q)
#}
#
#
#XXrs=crossprod(stsub)/(1011)
#XYrs=crossprod(stsub,presid)/(1011)
#
#gind.cut=split(1:ncol(presid), cut(1:ncol(presid),11))
#s2.m=mcMap(function(i) { scantwo.additive.noInv(XYrs[,i], XXrs, mspace=5) }, gind.cut, mc.cores=11)
#s2.m=do.call('rbind', s2.m)
#
#
#plot(0,0, type='n', xlim=c(0, ncol(stsub)), ylim=c(0,6000))
#for(i in 1:nrow(s2)) {
#    points(s2[i,1],i)
#    points(s2[i,2],i)
#}
#
#gsind.cut=split(1:length(which(s2[,3]>.2)), cut(1:length(which(s2[,3]>.2)),11))
#
#s3=mcMap(function(i) { scanthree.additive.noInv(XYrs, XXrs, mspace=5) }, gind.cut, mc.cores=11)
#s3=do.call('rbind', s3)
#
#plot(0,0, type='n', xlim=c(0, ncol(stsub)), ylim=c(0,nrow(s3)))
#for(i in 1:nrow(s3)) {
#    points(s2[i,1],i)
#    points(s2[i,2],i)
#}
#
#
#
#
#
#
#
#}
#
#
#
#
#plot(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2]), vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2]), 
#     xlab='marker subset', ylab='all markers', main = 'r2>.998' )
#x11()
#plot(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2]), vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2]), 
#     xlab='marker subset', ylab='all markers', main = 'r2>.996' )
#x11()
#plot(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2]), vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2]), 
#     xlab='marker subset', ylab='all markers', main = 'r2>.994' )
#x11()
#plot(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2]), vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2]), 
#     xlab='marker subset', ylab='all markers', main = 'r2>.992' )
#x11()
#plot(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2]), vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2]), 
#     xlab='marker subset', ylab='all markers', main = 'r2>.99' )
#
#
#
#cor(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2]), vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2]))
#
#hist(-(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2])-vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2])), breaks=100)
#
#mean(-(vcEstMapS[,1]/(vcEstMapS[,1]+vcEstMapS[,2])-vcEstMap[,1]/(vcEstMap[,1]+vcEstMap[,2])))
#
#
#for(i in 1:6000) {
#    print(i)
#    zz=regsubsets(stsub, presid[,i], int=F, nbest=1, nvmax=4, really.big=T, all.best=F, method='exhaustive')
#}
#
#
#
#
#
#
#coef.inds1=matrix(0,20,1)
#coef.inds2=matrix(0,20,2)
#coef.inds3=matrix(0,20,3)
##coef.inds4=matrix(0,20,4)
##coef.inds5=matrix(0,20,5)
##coef.inds6=matrix(0,20,6)
#
#r2ind=matrix(0,20,6)
##bicind=matrix(0,20,4)
##cpind=matrix(0,20,4)
##adjr2ind=matrix(0,20,4)
#
#for(i in 1:20) {
#lp1=regsubsets(stpreds[,seq(i,794,2)], sy[,2], int=F, nbest=1, nvmax=2, really.big=T, all.best=F, method='exhaustive')
#lp1=regsubsets(stsub, sy[,2], int=F, nbest=1, nvmax=3, really.big=T, all.best=F, method='exhaustive')
#lp1=leaps(stsub, sy[,2], int=F, nbest=1, method='r2')
#
#n1=match(names(which(summary(lp1)$which[1,])), colnames(stpreds))
#n2=match(names(which(summary(lp1)$which[2,])), colnames(stpreds))
#n3=match(names(which(summary(lp1)$which[3,])), colnames(stpreds))
#n4=match(names(which(summary(lp1)$which[4,])), colnames(stpreds))
#n5=match(names(which(summary(lp1)$which[5,])), colnames(stpreds))
#n6=match(names(which(summary(lp1)$which[6,])), colnames(stpreds))
#
#coef.inds1[i,]=n1
#coef.inds2[i,]=n2
#coef.inds3[i,]=n3
#coef.inds4[i,]=n4
#coef.inds5[i,]=n5
#coef.inds6[i,]=n6
#
#r2ind[i,]=summary(lp1)$rsq[1:6]
#bicind[i,]=summary(lp1)$bic[1:4]
#cpind[i,]=summary(lp1)$cp[1:4]
#adjr2ind[i,]=summary(lp1)$adjr2[1:4]
#}
#
#data.frame(n1=coef.inds1,n2=coef.inds2, n3=coef.inds3,n4=coef.inds4,n5=coef.inds5,n6=coef.inds6, r2=r2ind)
#
#, bic=bicind, cp=cpind, adjr2=adjr2ind)
#
#2112  252  386  486
#
#lp1=regsubsets(stpreds[,seq(1,794,2)], sy[,8], int=F, nbest=1, nvmax=3, really.big=T, all.best=F, method='exhaustive')
#match(names(which(summary(lp1)$which[3,])), colnames(stpreds))
#
#lp1=regsubsets(stpreds[,seq(2,794,2)], sy[,8], int=F, nbest=1, nvmax=3, really.big=T, all.best=F, method='exhaustive')
#match(names(which(summary(lp1)$which[3,])), colnames(stpreds))
#
#
#r21=summary(lp1)$rsq[1:3]
#coef.inds[i,]=n1
#r2ind[i,]=r21
#print(n1)
#}
#
#eg=expand.grid(coef.inds[,1], coef.inds[,2])
#egvec=rep(0, 8000)
#for(i in 1:nrow(eg)){
#    print(i)
#    egvec[i]=( sum(lm.fit(stpreds[,as.vector(unlist(eg[i,]))], sy[,8])$residuals^2))
#}
#
#
#
#
#
#coef.inds=matrix(0, 20,3)
#r2ind=matrix(0,20,3)
#for(i in 1:20) {
#lp1=regsubsets(stpreds[,seq(i,794,20)], sy[,8], int=F, nbest=1, nvmax=4, really.big=T, all.best=F, method='exhaustive')
#n1=match(names(which(summary(lp1)$which[3,])), colnames(stpreds))
#r21=summary(lp1)$rsq[1:3]
#coef.inds[i,]=n1
#r2ind[i,]=r21
#print(n1)
#}
#
#lp2=regsubsets(stpreds[,sort(as.vector(coef.inds))], sy[,2], int=F, nbest=1, nvmax=4, really.big=T, all.best=F, method='seqrep')
#match(names(which(summary(lp2)$which[4,])), colnames(stpreds))
#
#
#
#
#eg1=paste0('X',colnames(stpreds)[sort(coef.inds[,1])])
#eg2=paste0('X',colnames(stpreds)[sort(coef.inds[,2])])
#eg3=paste0('X',colnames(stpreds)[sort(coef.inds[,3])])
#std=data.frame(stpreds)
#
#test=regsubsets(sy[,8]~eg1+eg2+eg3,data=std, nbest=1)
#
#eg=expand.grid(coef.inds[,1], coef.inds[,2], coef.inds[,3])
#egvec=rep(0, 8000)
#for(i in 1:nrow(eg)){
#    print(i)
#    egvec[i]=( sum(lm.fit(stpreds[,as.vector(unlist(eg[i,]))], sy[,8])$residuals^2))
#}
#eg[2112,]
#     Var1 Var2 Var3
#2112  252  386  486
#
#
#
#lp2=regsubsets(stpreds[,seq(2,794,5)], sy[,1], int=F, nbest=1, nvmax=3, really.big=T, all.best=F, method='exhaustive')
#n2=names(which(summary(lp2)$which[3,]))
#
#lp3=regsubsets(stpreds[,seq(3,794,5)], sy[,1], int=F, nbest=1, nvmax=3, really.big=T, all.best=F, method='exhaustive')
#n3=names(which(summary(lp3)$which[3,]))
#
#lp4=regsubsets(stpreds[,seq(4,794,5)], sy[,1], int=F, nbest=1, nvmax=3, really.big=T, all.best=F, method='exhaustive')
#n4=names(which(summary(lp4)$which[3,]))
#
#lp5=regsubsets(stpreds[,seq(5,794,5)], sy[,1], int=F, nbest=1, nvmax=3, really.big=T, all.best=F, method='exhaustive')
#n5=names(which(summary(lp5)$which[3,]))
#
#
#sort(match(c(n1,n2,n3,n4,n5), colnames(stpreds)))
#
#lpf=regsubsets(stpreds[,sort(match(c(n1,n2,n3,n4,n5), colnames(stpreds)))], sy[,1], int=F, nbest=1, nvmax=3, all.best=F, method='exhaustive')
#names(which(summary(lpf)$which[3,]))
#
#
#which(summary(lp)$which[2,])
#
#
#
#dup=FALSE
#stpreds.orig=stpreds
#ds=duplicated(stpreds, MARGIN=2)
#while(sum(ds>0) ) {
#    print('removing markers')
#    print(ncol(stpreds))
#    stpreds=stpreds[,-which(ds)]
#    ds=duplicated(stpreds, MARGIN=2)
#}
#stpreds2=stpreds[,-which(duplicated(stpreds))]
#
#library(leaps)
#lp=regsubsets(stpreds[,seq(1,794,5)], sy[,1], int=F, nbest=1, nvmax=2, really.big=T, all.best=F)
#system.time({lp=regsubsets(stpreds[,seq(1,794,2)], sy[,1], int=F, nbest=1, nvmax=10, really.big=T, all.best=F, method='forward')})
#
#which(summary(lp)$which[3,])
#
#
#
#
#
#lp=leaps(stpreds[,seq(1,794,3)], sy[,1], int=F, nbest=1, nvmax=3, really.big=T, all.best=F)
#
#
#ssy=svd(sy)
#
#ssyL=fasterLOD(1012,scale(ssy$u[,1:50]), stpreds)
#
#sy
#ssyr=svd(scale(t.tpm.matrix))
#ssyR=fasterLOD(1012,scale(ssyr$u[,1:50]), gdata.scaled)
#for(i in 1:50) {
#    plot(ssyR[i,])
#    readline()
#}
#
#
#
#low <- 1
#high <- 1
#
##or Ben Bolker's better alternative
#delta <- col(XXr) - row(XXr)
#XXr[delta==1] 
#
#
##runs of correlated markers
#c(1,which(XXr[delta==1]>.9979)+1)
#
#which(XXr[delta==1]<.9979)+1
#
#
##let's try reducing XXr
#
#
#
#
#
#lshoot=list()
#i=1
#lr=lm.ridge(sy[,i]~stpreds, lambda=vcEstMap[i,2]/(vcEstMap[i,1]/794))
#lsr=lassoshooting(X=stpreds, y=sy[,i] , lambda=sqrt(vcEstMap[i,2]/(vcEstMap[i,1]/794)))
#lbetas=fasterLOD(1012, sy, stpreds, betas=TRUE)
#abs(lbetas$r[1,])
#library(glmnet)
#library(CDLasso)
#i=2
#lsr2=l1.reg(t(stpreds), sy[,i], lambda=sqrt(vcEstMap[i,2]/(vcEstMap[i,1]/794)))
#lsr3=l1.reg(t(stpreds), sy[,i], lambda=36)
#lsr4=slim(stpreds, sy[,i], method='lq',lambda=.09)
#l
#
#for(i in 1:2000) {
#    print(i)
#    #e/(a/n) is lambda for RR
#    #sqrt for lasso?
#    lshoot[[as.character(i)]]=lassoshooting(X=stpreds, y=sy[,i] ,lambda=sqrt(vcEstMap[i,2]/(vcEstMap[i,1]/794)))
#}
#
#
#
#lasso.list=list()
#for(i in 1:2000) {
#    print(i)
#   lasso.list[[as.character(i)]]=cv.EBglmnet(stpreds, sy[,i] ,verbose=1, nfolds=5)
#  #  lasso.list[[as.character(i)]]=EBglmnet(stpreds, sy[,i] ,prior='lasso', hyperparameters=c(.01, verbose=1)
#
#}
#i=1
#
#lasso.list2=list()
#for(i in 1:6000) {
#    print(i)
#   lasso.list2[[as.character(i)]]=EBglmnet(stpreds, sy[,i] ,verbose=1, hyperparameters=c(-.1,1))
#}
#
#
#
#i=1
#x11()
#zz=sapply(lasso.list2, function(x) x$fit)
#plot(0,0, type='n', xlim=c(0, ncol(gforsim)), ylim=c(0,6000))
#for( i in 1:6000) {
#    points(zz[[i]][,1][zz[[i]][,6]<.0001], rep(i, sum(zz[[i]][,6]<.0001)))
#}
#
#
#for(i in 1:2000) {
#plot(simL[i,])
#points(lshoot[[i]]$coefficients*100, col='red')
#abline(v=add.qtl.ind)
#readline()
#}
#glmn=cv.glmnet(stpreds, sy[,1],nfolds=5)
#
#
#plot(lshoot[[2]]$coefficients)
#
#
#
#plot(lassoshooting(X=stpreds, y=sy[,3] ,lambda=sqrt(.91/(.09/794)))$coefficients)
#abline(v=add.qtl.ind)
#
#
#
#testc=EBglmnet(stpreds, sy[,i] ,verbose=1, hyperparameters=c(-.1,1))
#
#test=lassoshooting(X=stpreds, y=sy[,1] ,lambda=.91/.01)
##test=lassoshooting(X=stpreds, y=sy[,1] ,lambda=(.91/.09/794))
#test=lassoshooting(X=stpreds, y=sy[,1] ,lambda=sqrt(.91/(.09/794)))
#
#
#
#
#
#cvE=cv.EBglmnet(stpreds, sy[,4] ,verbose=1)
#
#
#
#ssub=seq(1,794,4)
#Xs=(stpreds[,ssub])
#XXsr=crossprod(Xs)/(n-1)
#A2=tcrossprod(stpreds[,ssub])/length(ssub)
#
#XYr=crossprod(Xs,sy)/(n-1)
#t3=mr3.scantwo.BI.additive.2(XYr[,1:5], XXsr, 5)
#gind.cut=split(gind, cut(gind,11))
#r2.3=mcMap(function(i){ mr3.scantwo.BI.additive.2(XYr[,i], XXr,n.marker.space) }, gind.cut, mc.cores=11)
#r23=do.call('rbind', r2.3)
#x11()
#
#plot(0,0, type='n', xlim=c(0, ncol(gforsim)), ylim=c(0,2000))
#    points(r23[,1], 1:1969)
#    points(r23[,2], 1:1969)
#    points(r23[,3], 1:1969)
#
#
#
#
#
#mind=apply(cbind(apply(simL,1 , which.max)-add.qtl.ind[1],
#apply(simL,1 , which.max)-add.qtl.ind[2],
#apply(simL,1 , which.max)-add.qtl.ind[3]),1, function(s) min(abs(s)))
#
#
#
#
#)
#
#
#
#
#r=regress(sy[,1]~1, ~A)
#r2=regress(sy[,1]~1, ~A2)
#
#
#sqrt(.09)/794)
#test=lars(x=stpreds, y=sy[,1],intercept=FALSE)
#testcv=cv.lars(x=stpreds, y=sy[,1], K=5)
#A=tcrossprod(stpreds)/ncol(stpreds)
#tB=BGLR(sy[,1], response_type='gaussian', ETA=list(list(X=stpreds, model='BayesC')),R2=.09)
#plot(tB)
#
#plot(tB$ETA[[1]]$b)
#
#
#zz=sapply(lasso.list, function(x) x$fit)
#plot(0,0, type='n', xlim=c(0, ncol(gforsim)), ylim=c(0,107))
#for( i in 1:204) {
#    points(zz[[i]][,1][zz[[i]][,6]<.01], rep(i, sum(zz[[i]][,6]<.01)))
#}
#cvR=rvm(stpreds, sy[,1], kpar=list())
#
#cvRe=sigest(sy[,1]~stpreds)
#
#
#
#system.time({     test21=mr2.scantwo.additive.noInv(XYr[,c(1:50)], XXr, n.marker.space)})
#
#system.time({     test3=mr3.scantwo.additive(XYr[c(1:300),c(1:2)], XXr[c(1:300),c(1:300)], n.marker.space)})
#system.time({     test3=mr3.scantwo.additive(XYr[c(1:300),c(1:2)], XXr[c(1:300),c(1:300)], n.marker.space)})
#system.time({     test3.B=mr3.scantwo.BI.additive(XYr[c(1:300),c(1:2)], XXr[c(1:300),c(1:300)], n.marker.space)})
#system.time({     test3.B.2=mr3.scantwo.BI.additive.2(XYr[c(1:300),c(1:2)], XXr[c(1:300),c(1:300)], n.marker.space)})
#
#stpreds0=stpreds
#
#calc.additive.pairwise=function(p, g, n.marker.space=20, n.perms=5) {
#    n.perms=5
#    n.marker.space=15
#    sy=scale(p)
#    g=stsub
#    stpreds=scale(g)
#    n=nrow(p)
#    XYr=crossprod(stpreds,sy)/(n-1)
#    sy.p=replicate(n.perms, sy[sample(1:nrow(stpreds)),])
#   
#    permXYr=list()
#    permXYr[[1]]=XYr
#    for(nn in 1:n.perms) {
#        permXYr[[nn+1]]=crossprod(stpreds, sy.p[,,nn])/(n-1)
#    }
#    #permXYr=replicate(n.perms, crossprod(stpreds[sample(1:nrow(stpreds)),],sy)/(n-1))
#
#    allXYr=abind(permXYr, along=3)
#    XXr=crossprod(stpreds)/(n-1)
#    
#    r2.a=mcMap(function(i){ scantwo.additive.noInv(allXYr[,,i], XXr,n.marker.space) }, 1:(n.perms+1), mc.cores=11)
#    r2.aa=do.call('abind', c( r2.a, along=3))
#    dimnames(r2.aa)[[1]]=colnames(p)
#    dimnames(r2.aa)[[2]]=c('i', 'j', 'LOD')
#    dimnames(r2.aa)[[3]]=c('obs', paste0('p', 1:(dim(r2.aa)[3]-1)))
#
#    r2.aa[,3,]=-(n*log(1-r2.aa[,3,]))/4.605
#    #nm=(-n*log(1-XYr^2))/(2*log(10))
#    nm=(-n*log(1-allXYr^2))/(2*log(10))
#    am.nm=matrix(0,ncol(p),dim(r2.aa)[3])
#    colnames(am.nm)=dimnames(r2.aa)[[3]]
#    for(i in 1:ncol(p)) {
#        for(j in 1:dim(r2.aa)[3]) { 
#        am.nm[i,j]=max( nm[r2.aa[i,1,j],i,j], nm[r2.aa[i,2,j],i,j]) 
#        }
#    }
#    a2.v.a=r2.aa[,3,]-am.nm
#
#    #which((am-am.nm)>4)
#    #w2d=colnames(p)[which((am-am.nm)>4)]
#    #gene.annot.df[match(w2d,gene.annot.df$name),]
#    #plot(r2.aa[which((am-am.nm)>4),1,'obs'], r2.aa[which((am-am.nm)>4),2,'obs'])
#
#    obsPcnt=sapply(seq(1.5,5,.01), function(thresh) {sum(a2.v.a[,1]>thresh)})
#    names(obsPcnt)=seq(1.5,5,.01)
#    expPcnt = sapply(seq(1.5, 5, .01),  
#                         function(thresh) { 
#                                #print(thresh); 
#                                mean(apply(a2.v.a[,2:ncol(a2.v.a)], 2, function(ll) {sum(ll>thresh) }) )
#                            } )
#    names(expPcnt) = seq(1.5, 5, .01)
#    pFDR = expPcnt/obsPcnt
#    fdrFX=approxfun(pFDR, seq(1.5,5,.01))
#    thresha2=fdrFX(.05)
#    #R> thresha2
#    #[1] 3.22
#    rsig=r2.aa[,,1][a2.v.a[,1] > thresha2,]
#
#
#    plot(rsig[,1],1:nrow(rsig), xlim=c(0,ncol(gforsim)), col='#00000055',cex=log(rsig[,3],base=4.5))
#    points(rsig[,2],1:(nrow(rsig)),col='#00000055',cex=log(rsig[,3],base=4.5))
#
#    
#
#
#
#   o.v.t= cbind(peak.gene.ind[which(a2.v.a[,1]>thresha2)], r2.aa[which(a2.v.a[,1]>thresha2),,'obs'])
#   ghosts=(o.v.t[,3]> o.v.t[,1] & o.v.t[,2]< o.v.t[,1])
#   gind=which(a2.v.a[,1]>thresha2)
#   is.on.chr=gene.annot.df[match(names(gind), gene.annot.df$name),'chr']=='chrXIV'
#
#   new.col=ghosts+is.on.chr
#
#   plot( peak.gene.ind[which(a2.v.a[,1]>thresha2)][is.on.chr], r2.aa[which(a2.v.a[,1]>thresha2),1,'obs'][is.on.chr], 
#        xlim=c(0,900), ylim=c(0,900), xlab='one locus peak', ylab='two locus peaks', col=ghosts+1, pch=19)
#   points( peak.gene.ind[which(a2.v.a[,1]>thresha2)], r2.aa[which(a2.v.a[,1]>thresha2),1,'obs'], 
#        xlim=c(0,900), ylim=c(0,900), xlab='one locus peak', ylab='two locus peaks', col=ghosts+1)
#
#   points( peak.gene.ind[which(a2.v.a[,1]>thresha2)][is.on.chr], r2.aa[which(a2.v.a[,1]>thresha2),2,'obs'][is.on.chr],col=ghosts+1, pch=19)
#   points( peak.gene.ind[which(a2.v.a[,1]>thresha2)], r2.aa[which(a2.v.a[,1]>thresha2),2,'obs'],col=ghosts+1)
#
#   segments( peak.gene.ind[which(a2.v.a[,1]>thresha2)], r2.aa[which(a2.v.a[,1]>thresha2),1,'obs'],
#                peak.gene.ind[which(a2.v.a[,1]>thresha2)], r2.aa[which(a2.v.a[,1]>thresha2),2,'obs'], lwd=.5)
#
#
#    #gind.cut=list(a=1:2, b=1:2, c=1:2,d=1:2, e=1:2, f=1:2, g=1:2, h=1:2, i=1:2, j=1:2, k=1:2)
#    #r2.3=mcMap(function(i){ mr3.scantwo.additive(XYr[1:500,i], XXr[1:500,1:500],n.marker.space) }, gind.cut, mc.cores=11)
#
#   
#    gind.cut=split(gind, cut(gind,11))
#    r2.3=mcMap(function(i){ mr3.scantwo.additive(XYr[,i], XXr,n.marker.space) }, gind.cut, mc.cores=11)
#    r2.3p1=mcMap(function(i){ mr3.scantwo.additive(allXYr[,i,2], XXr,n.marker.space) }, gind.cut, mc.cores=11)
#    r2.3ap1=do.call('rbind', r2.3p1)
#    r2.3p2=mcMap(function(i){ mr3.scantwo.additive(allXYr[,i,3], XXr,n.marker.space) }, gind.cut, mc.cores=11)
#    r2.3ap2=do.call('rbind', r2.3p2)
#    r2.3p3=mcMap(function(i){ mr3.scantwo.additive(allXYr[,i,4], XXr,n.marker.space) }, gind.cut, mc.cores=11)
#    r2.3ap3=do.call('rbind', r2.3p3)
#    r2.3ap=do.call('abind', c(list(r2.3ap1, r2.3ap2, r2.3ap3), along=3 ))
#
#    r2.3a=do.call('rbind', r2.3)
#    colnames(r2.3a)=c('i', 'j', 'k', 'LOD')
#    rownames(r2.3a)=names(gind)
#    max.r2.3an=rep(NA, nrow(r2.3a))
#    for(i in 1:nrow(r2.3a)) {     
#        print(i)
#        tw1=summary(lm(sy[,gind[i]]~stpreds[,r2.3a[i,1]]+stpreds[,r2.3a[i,2]]-1))$r.squared
#        tw2=summary(lm(sy[,gind[i]]~stpreds[,r2.3a[i,1]]+stpreds[,r2.3a[i,3]]-1))$r.squared
#        tw3=summary(lm(sy[,gind[i]]~stpreds[,r2.3a[i,2]]+stpreds[,r2.3a[i,3]]-1))$r.squared
#        max.r2.3an[i]=max(-n*(log(1-c(tw1,tw2,tw3)))/4.605)
#    }
#    r2.3a[,4]=((-n*log(1-r2.3a[,4]))/(4.605) )
#   
#    a3.v.a2=r2.3a[,4] - max.r2.3an
# 
#    max.r2.3ap=matrix(0,nrow(r2.3a),dim(r2.3ap)[3])
#    for(j in 1:dim(r2.3ap)[3] ) {
#        for(i in 1:nrow(r2.3a)) {     
#            print(i)
#            tw1=summary(lm(sy.p[,gind[i],j]~stpreds[,r2.3ap[i,1,j]]+stpreds[,r2.3ap[i,2,j]]-1))$r.squared
#            tw2=summary(lm(sy.p[,gind[i],j]~stpreds[,r2.3ap[i,1,j]]+stpreds[,r2.3ap[i,3,j]]-1))$r.squared
#            tw3=summary(lm(sy.p[,gind[i],j]~stpreds[,r2.3ap[i,2,j]]+stpreds[,r2.3ap[i,3,j]]-1))$r.squared
#            max.r2.3ap[i,j]=max(-n*(log(1-c(tw1,tw2,tw3)))/4.605)
#        }
#    }
#    perm3.max.lods=(-n*log(1-r2.3ap[,4,])/4.605)
#    r3.perm.lods= perm3.max.lods - max.r2.3ap
#   
#    obsPcnt=sapply(seq(1.5,5,.01), function(thresh) {sum(a3.v.a2>thresh)})
#    names(obsPcnt)=seq(1.5,5,.01)
#    expPcnt = sapply(seq(1.5, 5, .01),  
#                         function(thresh) { 
#                                #print(thresh); 
#                                mean(apply(r3.perm.lods, 2, function(ll) {sum(ll>thresh) }) )
#                            } )
#    names(expPcnt) = seq(1.5, 5, .01)
#    pFDR = expPcnt/obsPcnt
#    fdrFX=approxfun(pFDR, seq(1.5,5,.01))
#    thresha3=fdrFX(.05)
#
#
#   o.v.t= cbind(peak.gene.ind[which(a2.v.a[,1]>thresha2)], r2.aa[which(a2.v.a[,1]>thresha2),,'obs'])
#   
#    #r2.3p4=mcMap(function(i){ mr3.scantwo.additive(allXYr[,i,5], XXr,n.marker.space) }, gind.cut, mc.cores=11)
#    #r2.3p5=mcMap(function(i){ mr3.scantwo.additive(allXYr[,i,6], XXr,n.marker.space) }, gind.cut, mc.cores=11)
#
#
#    obsPcnt=sapply(seq(1.5,5,.01), function(thresh) {sum(a2.v.a[,1]>thresh)})
#    
#    test.scan3=mr3.scantwo.additive(XYr[,gind.sub],XXr,n.marker.space)
#
#
#   
#   YPR036W   470 386 482 124.895
#    boxplot(p[,'YPR036W']~g[,470], names=c('470_R', '470_B'))
#    boxplot(p[,'YPR036W']~g[,386]+g[,482], names=c('386_R','386_B', '482_R','482_B'))
#
#
#    #positions vs uni
#
#    test=sy[,3426]
#
#   summary( lm(sy[,'YAL035W']~stpreds[,430]-1))
#   summary(    lm(sy[,'YAL035W']~stpreds[,274]+stpreds[,406]-1))
#   #summary(    lm(sy[,'YAL035W']~stpreds[,274]+stpreds[,406]+stpreds[,430]-1))
#
#   summary(    lm(sy[,'YAL035W']~stpreds[,274]+stpreds[,375]+stpreds[,480]-1))
#   summary(    lm(sy[,'YAL035W']~stpreds[,407]+stpreds[,431]+stpreds[,481]-1))
#
#
#   ltest=lars(stpreds, sy[,'YAL035W'], intercept=FALSE)
#   cvltest=cv.lars(stpreds, sy[,'YAL035W'], intercept=FALSE,trace=TRUE, type='lasso')
#   ltest3=lassoshooting(stpreds, sy[,'YAL035W'], lambda=75)
#   plot(ltest3$coefficients)
#   abline(v=c(274,375,480))
#   
#   require(doMC)
#   registerDoMC(cores=6)
#   ltest4=cv.glmnet(stpreds, sy[,'YAL035W'], parallel=TRUE, intercept=FALSE)
#   ltest4$lambda.min
##ltest4$lambda.min
#
#
#    summary(lm(sy[,'YAL035W']~stpreds[,c(274,375,406,430,480,794)]-1))
#
#R> which(ltest3$coefficients>0)
#[1] 375 406 430 480 794
#R> which(ltest3$coefficients<0)
#[1] 274
#
#
#
#
#}
#
#image.plot(XXr)
#
#dXXr=as.dist(1-XXr^2)
#fit <- hclust(dXXr, method="ward.D")
##plot(fit) # display dendogram
#cg=cutree(fit, k=ceiling(ncol(stpreds)/3))
#abline(v=cumsum(rle(cg)$lengths)/794)
#
#A.chr=tcrossprod(stpreds)/ncol(stpreds)
#
#Agroup=which(cg==7 | cg==11 | cg==12)
#Ao=tcrossprod(stpreds[,Agroup])/length(Agroup)
#
#regress(sy[,1]~1, ~A,verbose=T)
#regress(sy[,1]~1, ~Ao, verbose=T)
#regress(sy[,1]~1, ~A.chr.list[[7]]+A.chr.list[[11]]+A.chr.list[[12]],verbose=T)
#
#A2=tcrossprod(stpreds[,add.qtl.ind])/3
#regress(sy[,1]~1, ~A2,verbose=T)
#
#
#
#
#
#
#
#
#calcWindowVCs=function(cg, stpreds, sy, n.cores=5){
#    vcA.chr=list()
#    for(n in 1:max(cg)) {
#        print(n)
#        cluster.inds = which(cg==n)
#        A.chr=tcrossprod(stpreds[,cluster.inds])/length(cluster.inds)
#        eigA=doEigenA_forMM(sy,A.chr)
#        
#        # calculate mixed model, one term for additive variance  -------------------------------------------
#        mcvc=mcMap(function(i) m.S(sy[,i], K=A.chr,  theta=eigA$theta, Q=eigA$Q), i=colnames(sy), mc.cores=n.cores)
#        mcvc=do.call('rbind', mcvc)
#        mcvc1=mcvc[,1]/(mcvc[,1]+mcvc[,2])
#        vcA.chr[[as.character(n)]]=mcvc1
#    }
#    allVC.windows=do.call('cbind', vcA.chr)
#    return(allVC.windows)
#}
#obsVCwindows=calcWindowVCs(cg, stpreds, sy)
#expVCwindows=replicate(5, calcWindowVCs(cg, stpreds, sy[sample(1:nrow(sy)),]))
#
#
#A.chr.list=list()
# for(n in 1:max(cg)) {
#        print(n)
#        cluster.inds = which(cg==n)
#        A.chr.list[[as.character(n)]]=tcrossprod(stpreds[,cluster.inds])/length(cluster.inds)
#}
#
#fitall=regress(sy[,1]~1, ~A.chr.list[[1]]+A.chr.list[[2]]+A.chr.list[[3]]+A.chr.list[[4]]+A.chr.list[[5]]+
#A.chr.list[[6]]+A.chr.list[[7]]+A.chr.list[[8]]+A.chr.list[[9]]+A.chr.list[[10]]+
#A.chr.list[[11]]+A.chr.list[[12]]+A.chr.list[[13]]+A.chr.list[[14]]+A.chr.list[[15]]+
#A.chr.list[[16]]+A.chr.list[[17]]+A.chr.list[[18]]+A.chr.list[[19]]+A.chr.list[[20]], verbose=T, pos=rep(T,21))
#-471.7
#
#
##target----
#-470.2
##----------
#
#fitall=regress(sy[,1]~1, ~A.chr.list[[1]]+A.chr.list[[2]]+A.chr.list[[3]]+A.chr.list[[4]]+A.chr.list[[5]]+
#A.chr.list[[6]]+A.chr.list[[7]]+A.chr.list[[8]]+A.chr.list[[9]]+A.chr.list[[10]]+
#A.chr.list[[11]]+A.chr.list[[12]]+A.chr.list[[13]]+A.chr.list[[14]]+A.chr.list[[15]]+
#A.chr.list[[16]]+A.chr.list[[17]]+A.chr.list[[18]]+A.chr.list[[19]]+A.chr.list[[20]], verbose=T, pos=rep(T,21))
#
##threshold at sigma >.001
#fitsome=regress(sy[,1]~1, ~A.chr.list[[7]]+A.chr.list[[9]]+A.chr.list[[11]]+A.chr.list[[12]]+A.chr.list[[16]]+A.chr.list[[19]], verbose=T, pos=rep(T,7))
#fitsome=regress(sy[,1]~stpreds[,add.qtl.ind[3]], ~A.chr.list[[7]]+A.chr.list[[9]]+A.chr.list[[11]]+A.chr.list[[12]]+A.chr.list[[16]]+A.chr.list[[19]], verbose=T, pos=rep(T,7))
#
#
#
#471.8
#
#fitsome=regress(sy[,1]~1, ~A.chr.list[[7]]+A.chr.list[[9]]+A.chr.list[[11]]+A.chr.list[[12]], verbose=T, pos=rep(T,5))
#
#Agroup2=which(cg==7 |  cg==11 | cg==12 )
#Ao2=tcrossprod(stpreds[,Agroup2])/length(Agroup2)
#
#fitg1=regress(sy[,1]~1, ~Ao2, verbose=T)
#-472.4             
#               
#               +A.chr.list[[4]]+A.chr.list[[5]]+
#A.chr.list[[6]]+A.chr.list[[7]]+A.chr.list[[8]]+A.chr.list[[9]]+A.chr.list[[10]]+
#A.chr.list[[11]]+A.chr.list[[12]]+A.chr.list[[13]]+A.chr.list[[14]]+A.chr.list[[14]]+
#A.chr.list[[16]]+A.chr.list[[17]]+A.chr.list[[18]]+A.chr.list[[19]]+A.chr.list[[20]], verbose=T, pos=rep(T,21))
#
#
#
#    for(i in 1:ncol(sy)){
#        if(is.na(sd(sy[,i]))) {
#        next;
#        }
#        vcA.chr[[as.character(n)]][i,]=m.S(sy[,i], K=A,  theta=eigA$theta, Q=eigA$Q)
#    }


#    h2=(vcA[,1]/(vcA[,1]+vcA[,2]))
#    print(median(h2))
#    vcA.downsample[[pheno]]=h2
#}
 #k2=Mclust(t(tmm), G=1:9)
           #for( k in 2:12) {
           #     k2=kmeans(t(tmm), 3)
           #}
           #gchr=gdata.s.by.chr[[cc]]
           #Achr=tcrossprod(gchr)/ncol(gchr)
           #vcAchr=calcA(t(k2$centers), Achr)
           #h2chr=vcAchr[,1]/(vcAchr[,1]+vcAchr[,2])
#findHotspots.iteration2=function(hotspots, t.tpm.matrix, gene.annot.df, peaks.per.chr, genetic.map,
#                                 gdata, gdata.s.by.chr, gbatch.fact, do.save=T) {
#    hotspots.refined=list()
#    total.hotspot.n=0
#    unique.chrs=paste0('chr', as.roman(1:16))
#    for(cc in unique.chrs) { 
#        # attempt to refind hotspots given detected hotspots
#        ppc=peaks.per.chr[[cc]]
#        print(cc)
#        keep.transcripts=colnames(t.tpm.matrix)
#        # remove local QTL
#        keep.transcripts=keep.transcripts[gene.annot.df$chr!=cc]  
#        ppc2=ppc[ppc$gene %in% keep.transcripts,]
#        keep.transcripts=ppc2$gene
#        keep.transcripts=unique(keep.transcripts)
#            
#        #g=gdata.by.chr[[cc]]
#        g.s=gdata.s.by.chr[[cc]]
#        g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
#           
#        pmarker.all=as.character(sapply(hotspots[[cc]], function(x)x$peak))
#        pmarker.new=pmarker.all
#        #peaks=c()
#        #do while
#        for(pm in 1:length(pmarker.all)) {
#            pmarker=pmarker.new[-pm]
#            tmm=model.QTL.effects(t.tpm.matrix, gdata, keep.transcripts, peaks.per.gene, gbatch.fact, pmarker, cc)
#                  
#               tmp= fasterLOD(nrow(tmm),tmm,g.s)
#               max.obsLOD=apply(tmp,1,max)   
#              
#               LOD.kt=tmp #[skeep.transcripts,]
#               LOD.kt.wmax=apply(tmp, 1, which.max)
#               LOD.kt.max=apply(tmp, 1, max)
#               l15=list()
#               for(peak in 1:length(LOD.kt.max)){
#                         l15[[peak]]=range(which(LOD.kt[peak,]>(LOD.kt.max[peak]-1.5) ) )
#               }
#               l15=Intervals(do.call('rbind', l15))
#               io=interval_overlap(1:ncol(g.s), l15)
#               peak.index=which.max(sapply(io, length))  
#            
#               l.boundary=g.s.pos[match(pmarker.all[pm], colnames(g.s))]-50000
#               r.boundary=g.s.pos[match(pmarker.all[pm], colnames(g.s))]+50000
#               l.boundary.index  =  findInterval(l.boundary, g.s.pos)
#               r.boundary.index  =  findInterval(r.boundary, g.s.pos, rightmost.closed=F)
#               if(l.boundary.index<1) l.boundary.index=1
#               if(r.boundary.index>(ncol(g.s))) r.boundary.index=ncol(g.s)
#               #peaks.in.peak.bin=which((LOD.kt.wmax>=l.boundary.index) & (LOD.kt.wmax< r.boundary))
#               #to.scan=names(LOD.kt.max)[peaks.in.peak.bin]
#              
#               peaks.in.peak.bin=interval_overlap(peak.index, l15)[[1]]
#               #which((LOD.kt.wmax>=l.boundary.index) & (LOD.kt.wmax< r.boundary.index))
#               to.scan=names(LOD.kt.max)[peaks.in.peak.bin]
#           
#               to.scan=LOD.kt[to.scan,peak.index]
#               to.scan=sort(to.scan, decreasing=T)
#               if(length(to.scan)>500) { to.scan=to.scan[1:500] }
#               if(length(to.scan)<3)  {break;}
#
#               #to.scan=LOD.kt.max
#               #to.scan=sort(to.scan, decreasing=T)
#               #if(length(to.scan)>500) { to.scan=to.scan[1:500] }
#               #if(length(to.scan)<3)  {break;}
#               sub.t=tmm[,names(to.scan)]
#               #to.scan=skeep.transcripts
#               #gdata.s.by.chr[[cc]][, l.boundary.index] + gdata.s.by.chr[[cc]][,r.boundary.index])))
#               mmatrix=model.matrix(lm(sub.t[,1]~1))
#               sss=calcMVstat(g.s, sub.t, mmatrix, l.boundary.index, r.boundary.index)
#               print(names(which.max(sss))) 
#               hotspots.refined[[cc]][[pm]]=data.frame(maxneglog10p=max(sss), n.trans=length(to.scan), peak=names(which.max(sss)))
#               pmarker.new[pm]=names(which.max(sss))
#               total.hotspot.n=total.hotspot.n+1
#               
#               #update boundaries ----------------------------------------------------------------------------
#               l.boundary=g.s.pos[match(pmarker.new[pm], colnames(g.s))]-50000
#               r.boundary=g.s.pos[match(pmarker.new[pm], colnames(g.s))]+50000
#               l.boundary.index  =  findInterval(l.boundary, g.s.pos)
#               r.boundary.index  =  findInterval(r.boundary, g.s.pos, rightmost.closed=F)
#               if(l.boundary.index<1) l.boundary.index=1
#               if(r.boundary.index>(ncol(g.s))) r.boundary.index=ncol(g.s)
#               #-----------------------------------------------------------------------------------------------
#
#               attr(sub.t, 'chr')=cc
#               attr(sub.t, 'hs.n')=total.hotspot.n
#               attr(sub.t, 'l.boundary.index')=l.boundary.index
#               attr(sub.t, 'r.boundary.index')=r.boundary.index
#               if(do.save) {
#                    #save(sub.t, file=paste0('/media/juno/bootstrap_hotspot/', total.hotspot.n, '.RData' ))
#                    save(sub.t, file=paste0('/data/eQTL/RData/hotspots/', total.hotspot.n, '.RData2' ))
#               }
#        }
#    }
#    return(hotspots.refined)
#}



#hotspots=findHotspots.iteration1(t.tpm.matrix, gene.annot.df, 
#                                 peaks.per.chr, genetic.map, gdata, gdata.s.by.chr, gbatch.fact) 
##save(hotspots, file='/data/eQTL/RData/hotspots_042616.RData')
#load('/data/eQTL/RData/hotspots_042616.RData')
#
##Iteration 2 
#hotspots.refined=findHotspots.iteration2(hotspots, t.tpm.matrix, gene.annot.df, 
#                                         peaks.per.chr, genetic.map, gdata, gdata.s.by.chr, gbatch.fact ,do.save=T) 
##save(hotspots.refined, file='/data/eQTL/RData/hotspots_refined_060116.RData')
#load('/data/eQTL/RData/hotspots_refined_060116.RData')
##save(hotspots.refined, file='/data/eQTL/RData/hotspots_refined_042616.RData')
##load('/data/eQTL/RData/hotspots_refined_042616.RData')
#
#hotspot.boots=list()
#for( i in 1:85) {    
#    print(i)
#    hotspot.boots[[i]]=bootstrap.Hotspots(i)  
#}
## on Juno
#hotspot.boots=list()
#for( i in 43:1) {    #85:1 # 64:1 43:1
#    print(i)
#    hotspot.boots[[i]]=bootstrap.Hotspots(i, onjuno=TRUE)  
#}
#
##sapply(sapply(split(marker.GR, seqnames(marker.GR)), function(x) (ranges(x[length(x),]))), function(x) x[[1]][1])
##   chrI   chrII  chrIII   chrIV    chrV   chrVI  chrVII chrVIII   chrIX    chrX   chrXI  chrXII chrXIII  chrXIV   chrXV  chrXVI 
## 202825  798782  302979 1521369  563829  266145 1068261  519219  424866  724256  643662 1058607  914575  765330 1067520  928107
#
#hotspot.boot.peaks=list()
#hotspot.boot.intervals=matrix(NA, 85,5)
#for(nn in c(1:85)[-44]) {
#    setwd(paste0('/data/eQTL/RData/hotspots/', nn))
#    f=list.files('.')
#    r=c()
#    for(n in f) {load(n); r=c(r, peak) }
#    str(r)
#    r=gsub(':', '_', r)
#    pos=as.numeric(sapply(strsplit(r, '_'), function(x)x[2]))
#    hotspot.boot.peaks[[nn]]=pos
#    hotspot.boot.intervals[nn,]=c(quantile(pos, .025), quantile(pos, .05),  quantile(pos, .5),  quantile(pos, .95),quantile(pos, .975))
#}
#
#write.table(
#            data.frame(cbind(as.character(do.call('rbind', 
#                               lapply(hotspots.refined, function(x) do.call('rbind',x)))[,3]),hotspot.boot.intervals)),
#            file='/data/eQTL/RData/hotspots_062016.txt', sep='\t' ,quote=F, row.names=FALSE)
#phname=cbind(as.character(do.call('rbind', 
#                               lapply(hotspots.refined, function(x) do.call('rbind',x)))[,3]))
#pdf(file='/data/eQTL/RData/hotspots_ints.pdf', width=10, height=10)
#    for(nn in c(1:85)[-44]) {    hist(hotspot.boot.peaks[[nn]], breaks=100, main=phname[nn], sub=hotspot.boot.intervals[nn,4]-hotspot.boot.intervals[nn,2])
#      abline(v=hotspot.boot.intervals[nn,2], col='blue', lwd=2, lty=2)
#      abline(v=hotspot.boot.intervals[nn,4], col='blue', lwd=2, lty=2)
#    }
#dev.off()
#
#
#hist(peaks.per.chr[[11]]$pcind, breaks=100)
#abline(v=match(as.character(do.call('rbind', hotspots.refined[[11]])$peak), colnames(gdata.by.chr[[11]])), col='blue')
#
##----------------------------------------------------------------------------------------------------------
#big.plot1='/home/jbloom/Dropbox/Public/eQTL/eQTL_map_1024x1024.png'
#png(file=big.plot1, width=1024, height=1024)
#eQTL_bigPlot(all.peaks, gcoord.key, marker.GR)
#dev.off()
#
#big.plot2='/home/jbloom/Dropbox/Public/eQTL/eQTL_map_1920x1920.png'
#png(file=big.plot2, width=1920, height=1920)
#eQTL_bigPlot(all.peaks, gcoord.key, marker.GR)
#dev.off()
#
#big.plot3='/home/jbloom/Dropbox/Public/eQTL/eQTL_map_10x10.pdf'
#pdf(file=big.plot3, width=10, height=10)
#eQTL_bigPlot(all.peaks, gcoord.key, marker.GR)
#dev.off()
#
#for(chr in unique.chrs[-17]) {
#big.plot.by.chr=paste0('/home/jbloom/Dropbox/Public/eQTL/eQTL_map_', chr, '.pdf')
#pdf(file=big.plot.by.chr, width=10, height=10)
#eQTL_bigPlot(all.peaks, gcoord.key, marker.GR, xlim.ind= which(as.character(seqnames(marker.GR)) %in% chr)   )
#dev.off()
#
#}
#eQTL_bigPlot(all.peaks, gcoord.key, marker.GR, xlim.ind= which(as.character(seqnames(marker.GR)) %in% chr)   )
#xlim.ind=which(as.character(seqnames(marker.GR)) %in% 'chrVIII'
#htable=do.call('rbind', lapply(hotspots, function(x) do.call('rbind', x)))
#
#htabler=do.call('rbind', lapply(hotspots.refined, function(x) do.call('rbind', x)))
#png(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/030816/eQTL_hotspots_crunch.png', width=1024, height=256)
#eQTL_bigPlot(all.peaks, gcoord.key)
#abline(v=marker.GR$gcoord[match(htabler$peak, marker.GR$mname)])
#dev.off()
#
#big.plot.by.chr=paste0('/home/jbloom/Dropbox/Public/eQTL/eQTL_map_chrVIII_fixed_axis.pdf')
#pdf(file=big.plot.by.chr, width=10, height=10)
#
#dev.off()
# GO Enrichment
doGO=function(all.genes, set.of.interest, ontology='BP') {
#for( thisOntology in c("BP", "MF", "CC") ) {
    geneList=set.of.interest
    testGenes= factor(0+(all.genes %in% geneList))
    names(testGenes)=all.genes
    #colnames(t.tpm.matrix)
    GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
    GOresult = runTest(GOData, algorithm="classic", statistic="fisher")
    gt=GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult)))
    gt$result1=as.numeric(gsub('< ', '', gt$result1))
    genesinterms=genesInTerm(GOData, gt[,1])
    genes.enriched.list=lapply(genesinterms, function(x) x[x%in%names(testGenes[testGenes==1])])
    genes.enriched.list.simple=lapply(genes.enriched.list, function(x) as.character(SYS2ORF.key[x]))
    gt$Genes=as.vector(sapply(genes.enriched.list.simple, paste, collapse=','))
    gt$GenesSystematic=   as.vector(sapply(genes.enriched.list, paste, collapse=','))
    print(head(gt,20))
    # pdf(file=paste0('/data/eQTL/RData/GO/', 'epistatic_transcripts', '_', thisOntology, '.pdf'), width=25, height=25)
    #     plotGOToTree(GOData, GOresult, sigThres = 0.00005)
    #  dev.off()
    #write.table(gt, file=paste0('/data/eQTL/RData/GO/', 'epistatic_transcripts', '_', thisOntology, '.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
    return(gt)
}

#multivariate scanone
# returns logDetRSSfull
# 01 12 17 new hotspot code 
mvn.scanone=function(g.s, Y,add.cov=NULL, roi1=NULL) {
         
           ldetRSSf.1=rep(NA, ncol(g.s))
           if(!is.null(roi1)) { ss=roi1 } else { ss=1:ncol(g.s) } 
           for(i in min(ss):max(ss)) {
               #ldetRSSf.1[i]=det_AtA(residuals(.lm.fit(as.matrix(g.s[,i]), Y)))
               if(is.null(add.cov) ) {
                RSSf.1=crossprod(residuals(.lm.fit(as.matrix(g.s[,i]), Y)))
               } else {
                RSSf.1=crossprod(residuals(.lm.fit(cbind(add.cov, g.s[,i]), Y)))
               }
               ldetRSSf.1[i]=determinant(RSSf.1, logarithm=T)$modulus
           }
           return(ldetRSSf.1)
       }

# function to do scantwo
# add code to add additional covariates
    mvn.scantwo=function(g.s, Y, roi1 ,add.cov=NULL) {
        ldetRSSf.2=matrix(NA, length(roi1), length(roi1))
        #XX=as.matrix(add.cov)
        for(i in min(roi1):(max(roi1)-10) ) {
           for(j in (i+10):max(roi1) ){
            #print(i)
            # could add in additional covariates here
            if(is.null(add.cov) ) {
               RSSf.2=crossprod(residuals(.lm.fit(cbind(g.s[,i],g.s[,j]),Y)))
            } else{
                RSSf.2=crossprod(residuals(.lm.fit(cbind(add.cov,g.s[,i],g.s[,j]),Y)))
            }
            ldetRSSf.2[(i-(min(roi1)-1)),(j-(min(roi1)-1))]=determinant(RSSf.2, logarithm=T)$modulus
            #mlvec.1[i]=1012*(ldetRSSn.1-ldetRSSf.1)/(2*log(10))
            }
        }
           return(ldetRSSf.2)
     }

mvn.scanthree=function(g.s, Y, roi1, mv2Dpeak, add.cov=NULL) {
    ldetRSSf.3=rep(NA, ncol(g.s)) #length(roi))
     for(i in roi1) {
            if( (abs(i-mv2Dpeak[1])>4) & (abs(i-mv2Dpeak[2])>4) ) {
             if(is.null(add.cov) ) {
               RSSf.3=crossprod(residuals(.lm.fit(cbind(g.s[,mv2Dpeak[1]], g.s[,mv2Dpeak[2]], g.s[,i]),Y)))
            } else{
                RSSf.3=crossprod(
                                residuals(.lm.fit(cbind(add.cov, g.s[,mv2Dpeak[1]], g.s[,mv2Dpeak[2]], g.s[,i]),Y) )
                                 )
            }
            #   mtest=.lm.fit(cbind(g.s[,mv2Dpeak[1]], g.s[,mv2Dpeak[2]], g.s[,i]), Y) #~g.s[,i])    
            #RSSf.3=crossprod(residuals(mtest))
            ldetRSSf.3[i]=determinant(RSSf.3, logarithm=T)$modulus
            }
       }
    return(ldetRSSf.3)
}



getROI=function(mvLOD.1, g.s, n=40) {
         #roi1=range(which(mvLOD.1> max(mvLOD.1)-50))
         wm=which.max(mvLOD.1)
         m=max(mvLOD.1)
         if(m>1 & m<ncol(g.s)) {
            if(mvLOD.1[wm-1]==0 | mvLOD.1[wm+1]==0 ) { 
            mvLOD.1[wm]=0
            wm=which.max(mvLOD.1)
            }
         }
               
         roi.l1=ifelse( (wm-n)<1,1,wm-n)
         roi.r1=ifelse( (wm+n)>ncol(g.s), ncol(g.s), wm+n  )
         roi1=roi.l1:roi.r1 #[2] #[1]:roi[2]
         return(roi1 )#[1]:roi1[2])
     }


permKmeans=function(Y, n.kmeans=15, n.perm=10) {
    #perm.maxdiff=rep(NA, n.perm)
    perm.maxdiff=foreach(n=1:n.perm, .combine=c ) %dopar% {
    #Y=fC$r
        print(n)
        permr=apply(Y,1, function(x) sample(x) )
        kmean.listN=list()
        for(i in 1:n.kmeans) {
               kmean.listN[[as.character(i)]]=kmeans(permr, centers=i)
        }
        #perm.maxdiff=
        return(max(diff(sapply(kmean.listN, function(x) x$betweenss)/sapply(kmean.listN, function(x) x$totss))))
    }
    return(perm.maxdiff)
}

# refactor this ... calculating this when X is large needs to be slimmed down
calc.composite.MV.LOD=function(detected1D, Ysub, g.s)  {
    print('calculating composite LOD')
    compositeMVN.1D=list()
    for(ii in 2:ncol(detected1D) ) { 
        #colnames(detected1D)[-1] ) {
           i=colnames(detected1D)[ii]
           pm=match(i, colnames(detected1D))
           pnull=determinant(crossprod(residuals(.lm.fit(detected1D[,-pm],Ysub))), logarithm=T)$modulus
           pscan=mvn.scanone(g.s, Ysub, add.cov=detected1D[,-pm])
           compositeMVN.1D[[i]]=( 1012*(pnull-pscan)/(2*log(10)))

           print(paste(i, colnames(g.s)[which.max( compositeMVN.1D[[i]])]))
        }
       compositeMVN.1D=do.call('rbind',  compositeMVN.1D)
       rownames( compositeMVN.1D) = colnames(detected1D)[-1]
       return(compositeMVN.1D)
}

           #compositeMVN.1D=foreach(i=colnames(detected1D)[-1]  ) %dopar% {

plotMV.LOD=function(compositeMVN.1D, g.s, g.s.pos, steps=FALSE)  {
    p.detected=match(rownames(compositeMVN.1D), colnames(g.s))
    plot(g.s.pos, compositeMVN.1D[1,], ylim=c(0,max(compositeMVN.1D)), ylab='mvLOD', xlab='chromosome position')
     for(j in 2:nrow(compositeMVN.1D)){ 
         if(steps) {         readline();    }
         points(g.s.pos, compositeMVN.1D[j,], col=j)      }
     abline(v=g.s.pos[p.detected], col=1:nrow(compositeMVN.1D))
}

relocate.peaks=function(compositeMV.LOD, g.s, g.s.pos ,detected1D, max.move=75000) {
       # use this to relocalize signal
       
       peak.relocate=colnames(g.s)[apply(compositeMV.LOD,1, which.max)]
       
       peak.o.pos = g.s.pos[colnames(detected1D)[-1]]
       peak.n.pos = g.s.pos[apply(compositeMV.LOD,1, which.max)]
        
       detected1D.2=cbind(detected1D[,1], g.s[,match(peak.relocate, colnames(g.s))] )
       colnames(detected1D.2)=c('Intercept', peak.relocate)
       
       pdiff = abs(peak.o.pos-peak.n.pos)
       pdiff.ind= as.numeric(which(pdiff<max.move))+1
       return( detected1D.2[,c(1,pdiff.ind)] )

}


calcIndividualEffects_afterMV=function(tmms, detected1D, dfin){

      fullmodel=mclapply(colnames(tmms), {function(i) lm(tmms[,i]~.-1, data=dfin)}, mc.cores=60)
      names(fullmodel)=colnames(tmms)
      dropterms=mclapply(colnames(tmms), {function(i) dropterm(fullmodel[[i]], test='F') }, mc.cores=60)
      allp=lapply(dropterms, function(x)x$"Pr(F)"[-1])
      allpv=unlist(allp)
      
      bet=list()
      for(i in 1:ncol(tmms)) {
            #print(i)
            tokeep=which(allp[[i]]<.05)
            dfr=dfin[,tokeep]
            if(is.null(dim(dfr))) { dfr=data.frame(dfr); names(dfr)=colnames(dfin)[tokeep]; }
            if(ncol(dfr)>0 ) {
                bet[[colnames(tmms)[i]]]=(lm(tmms[,i]~.-1., data=dfr))
            }
      }
      return(bet)

}

doMV.bootstraps=function(detected1D, compositeMV.LOD, Ysub, g.s) {
    boots=list()
      for(peaks in colnames(detected1D)[-1] ) {
          print(peaks)
          pm = match(peaks, colnames(detected1D))
          roi1=getROI(compositeMV.LOD[(pm-1),],g.s)
          bootpos=rep(NA,1000)
          pb =txtProgressBar(min = 1, max = 1000, style = 3)

          for(bb in 1:1000) {
             smpme=sample(1:1012, replace=T)
             bnull=determinant(crossprod(residuals(lm(Ysub[smpme,]~1))), logarithm=T)$modulus
             pM=mvn.scanone(g.s[smpme,], Ysub[smpme,],add.cov=detected1D[smpme,-pm] , roi1 )
             bmvLOD.1=1012*(bnull-pM)/(2*log(10))
             setTxtProgressBar(pb, bb)
             bootpos[bb]=which.max(bmvLOD.1)
        }
        close(pb)
        boots[[peaks]]=bootpos
        print(colnames(g.s)[quantile(bootpos, c(.025, .975))])
      }
    return(boots)
}

#stat.mat is dim n x number of permutations
doFDR = function (stat.range, stat.vec, stat.mat, threshold) {
     stest=stat.range
     maxr2= stat.vec
     maxr2.perm= stat.mat

     obsPcnt=sapply(stest, function(thresh) sum(maxr2>thresh, na.rm=T))
     names(obsPcnt)=stest   
     expPcnt = sapply(stest,  
                         function(thresh) { 
                            mean(apply(maxr2.perm, 2, function(ll) {sum(ll>thresh, na.rm=T) } ))
                         })
     names(expPcnt) = stest
     
     pFDR = expPcnt/obsPcnt
     pFDR[is.na(pFDR)]=0
     pFDR = rev(cummax(rev(pFDR)))
    # if(sum(is.na(pFDR))>5) {
    #   print('this')
    #   return(999) } #break;}

   fdrFX=approxfun(pFDR, stest)
   #seq(1.5,6,.01))
   thresh=fdrFX(threshold)
   return(thresh)

}




#      dboot=names(dfin)[6]
#      tt=drop.terms(terms(bet2[[1]]), dropx=4, keep.response=T)
#      t1=update(bet2[[1]], tt )
#      a1[i]=add1(t1, ~. + g.s[,i], test='Chisq')[2,5] #

# do CV ... create 13 random splits ... keep track of splits, find peaks in 12 of 13 sets and then estimate variance explained in the 13th set
# refit full model (keep AOV info)
# hard crash?? start here
# recreate this
#load('/data/eQTL/RData/062016.RData')
# do cross validation -------------------------------------------------------------------------------------------------------------------(run once)
doCV = function(gbatch.fact, covariates.OD, t.tpm.matrix, pheno.scaled.OD, gdata) {
    n.groups=length(levels(gbatch.fact))
    for(icv in 1:n.groups) {
        rmgroup1=which(gbatch.fact %in% levels(gbatch.fact)[icv])
        covariates.OD.c1=covariates.OD[-rmgroup1,]
        t.tpm.matrix.c1=t.tpm.matrix[-rmgroup1,]
        pheno.scaled.OD.c1=pheno.scaled.OD[-rmgroup1,]
        gdata.c1=gdata[-rmgroup1,]
        gdata.scaled.c1=scale(gdata.c1)
        background.QTL.OD.c1 = find.background.QTL(covariates.OD.c1,t.tpm.matrix.c1, pheno.scaled.OD.c1, gdata.c1, gdata.scaled.c1)
        peakList.OD.c1=mapQTL(covariates.OD.c1, background.QTL.OD.c1,
                          t.tpm.matrix.c1, pheno.scaled.OD.c1, gdata.c1, gdata.scaled.c1,
                          n.perm=100, FDR.thresh=.05)
        saveRDS(peakList.OD.c1, file=paste0('/data/eQTL/RData/cross_validation_peaks',icv))
        print(sum(sapply(sapply(peakList.OD.c1, function(x) { sapply(x, function(y) nrow(y) ) } ), sum)))
    }
}
#sum(sapply(sapply(peakList.OD.c1, function(x) { sapply(x, function(y) nrow(y) ) } ), sum))



# load cross-validation ----------------------------------------------------------------------------------------------------------------------------
loadCV = function(gbatch.fact, covariates.OD, t.tpm.matrix, gdata) {
    n.groups=length(levels(gbatch.fact))
    cvVE=list()
    for(icv in 1:n.groups) {
        print(icv)
        peakList.OD.c1=readRDS(file=paste0('/data/eQTL/RData/cross_validation_peaks',icv))
        all.peaks.OD.c1=buildPeakListDF(peakList.OD.c1,gdata,gene.GR, marker.GR)
        peaks.per.gene.c1=split(all.peaks.OD.c1, all.peaks.OD.c1$gene)
        rmgroup1=which(gbatch.fact %in% levels(gbatch.fact)[icv])

        for(g in names(peaks.per.gene.c1)) {
        #g=names(peaks.per.gene.c1)[100]
            ppg=peaks.per.gene.c1[[g]][!duplicated(peaks.per.gene.c1[[g]]$pmarker),]
            #ppg=peaks.per.gene[[g]][!duplicated(peaks.per.gene[[g]]$pmarker),]
            
            apeaks = match(ppg$pmarker, colnames(gdata))
            #print(apeaks)
            #match(peaks.per.gene[[g]]$pmarker, colnames(gdata))

            X2=gdata[,apeaks]
            #X=data.frame(covariates.OD, gdata[,apeaks])
            yr=residuals(lm(t.tpm.matrix[,g]~covariates.OD))
            yr.train=yr
            yr.test=yr
            yr.train[rmgroup1]=NA
            fitme=lm(yr.train~.-1, data=data.frame(X2))
            if(is.null(dim(X2))){ 
                predicted=X2[rmgroup1]*coef(fitme)
            }else {
                predicted=X2[rmgroup1,]%*%coef(fitme)
            }
            cvVE[[g]]=c(cvVE[[g]], cor(yr.test[rmgroup1], predicted)^2)
        }
    }
    return(cvVE)
}
#-------------------------------------------------------------------------------------------------------------------
