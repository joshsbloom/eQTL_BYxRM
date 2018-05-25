load('/data/eQTL/RData/mv.bootstrap.pos.RData')
mmm=t(apply(mv.boots,1, function(x) quantile(x,c(.025, .975))))
#boot.GRanges=GRanges(seqnames= seqnames(marker.GR[mgiL]), ranges=IRanges(start=start(marker.GR[mgiL]), end=end(marker.GR[mgiR])))
#sapply(strsplit(rownames(mmm), ':'), function(x)x[1])
boot.GRanges=GRanges(seqnames=sapply(strsplit(rownames(mmm), ':'), function(x)x[1]), ranges=IRanges(start=mmm[,1], end=mmm[,2]))

acis=all.peaks.OD[all.peaks.OD$cis,]
atrans=all.peaks.OD[!all.peaks.OD$cis,]

#Extend gene.GR such that each gene interval starts and ends at the midpoint between genes
gsplit=split(gene.GR, as.vector(seqnames(gene.GR)))
gsplit=lapply(gsplit, function(x) {
    x=x[order(x$gcoord),]
    e=end(x)
    s=start(x)
    esd=((start(x[-1])-end(x))/2)
    esd[esd<0]=0
    e=e+(floor(esd))-1
    s[-1]=s[-1]-ceiling(esd[-length(esd)])
    start(x)=s
    end(x)=e
    return(x)
})
geneExtend.GR=stack(GRangesList(gsplit))


# two GRanges objects, e.g. return for each cis bin the number of trans overlapping 
buildOverlapVector=function(cis, trans, namevec) {
    fo=findOverlaps(cis, trans)
    sfo=split(subjectHits(fo), queryHits(fo))
    blength=sapply(sfo, length)
    bvec=rep(0, length(cis))
    bvec[as.numeric(names(blength))]=blength
    names(bvec)=namevec
    return(bvec)
}
#genome cis definition
#cisGenes=acis$gene                             

#pointwise cis definition
cisGenes=colnames(t.tpm.matrix)[qvalue(cisModel)$qvalue<.05]

cis.GRanges=geneExtend.GR[match(cisGenes, geneExtend.GR$ORF)]
notcis.GRanges=geneExtend.GR[-match(cisGenes, geneExtend.GR$ORF)]

# define trans by confidence intervals
trans.GRanges=GRanges(seqnames=atrans$chr, ranges=IRanges(start= as.numeric(sapply(strsplit(atrans$CI.l, ':|_'), function(x)x[2])),
                                                          end  = as.numeric(sapply(strsplit(atrans$CI.r, ':|_'), function(x)x[2]))), gene=atrans$gene)
# all trans that overlap a hotspot
hot.trans=(unique(queryHits(findOverlaps(trans.GRanges, boot.GRanges))))

# redefine trans by peak only                 
trans.GRanges=GRanges(seqnames=atrans$chr, ranges=IRanges(start= as.numeric(sapply(strsplit(atrans$pmarker, ':|_'), function(x)x[2])),
                                                          end  = as.numeric(sapply(strsplit(atrans$pmarker, ':|_'), function(x)x[2]))), gene=atrans$gene)

# all trans that do not fall in a hotspot 
notHot=trans.GRanges[-hot.trans]

#cis.GRanges2kb=GRanges(seqnames=acis$chr, ranges=IRanges(start= as.numeric(sapply(strsplit(acis$pmarker, ':|_'), function(x)x[2]))-100,
#                                                         end  = as.numeric(sapply(strsplit(acis$pmarker, ':|_'), function(x)x[2]))+100), gene=acis$gene, gcoord=acis$marker.gcoord)
# Cleanup analysis of trans not in hotspots ... falling on cis and falling on genome bins (with permutations)

#genome bins
gbins=sapply(contig.lengths, function(x) {
           bin.size=10000
           # x=round(x, digits=-4)
           yy=seq(bin.size,x,bin.size)
           xx=seq(1,x, bin.size)
           if(length(xx)>length(yy)) {xx=xx[-length(xx)] } else {yy=yy[-length(yy)] }
           cbind(xx,yy)
})
G.Granges=GRanges(seqnames=rep(names(gbins), sapply(gbins, nrow)), ranges=IRanges(start=unlist(sapply(gbins, function(x) x[,1])),
                                                                                    end=unlist(sapply(gbins, function(x) x[,2]))))

# Redfine intervals to not overlap bootstrap 
G.Granges=G.Granges[-unique(queryHits(findOverlaps(G.Granges, boot.GRanges)))]

cis.GRanges=cis.GRanges[-unique(queryHits(findOverlaps(cis.GRanges, boot.GRanges)))]
notcis.GRanges=notcis.GRanges[-unique(queryHits(findOverlaps(notcis.GRanges, boot.GRanges)))]
nhotmarker.GR=marker.GR[-unique(queryHits(findOverlaps(marker.GR, boot.GRanges)))]


### Excluding bins that overlap hotspots!!!
# cis as pointwise 

# overlap of "not hotstpot" trans with genome bins
notHot_genome=buildOverlapVector(G.Granges, notHot, '') 
plot(gcoord.key[as.vector(seqnames(G.Granges))]+start(G.Granges), notHot_genome, type='h')
abline(v=marker.GR[match(rownames(mv.boots), marker.GR$mname)]$gcoord, col='red', lty=3)
abline(v=gcoord.key, lty=2, col='lightblue', lwd=3)

#not Hot G.GRanges
#not.hot.GRanges=G.Granges
#not.hot.GRanges$not_hot_count=as.vector(notHot_genome)
#fnho=findOverlaps(G.Granges, notHot)
#nhgr=notHot[subjectHits(fnho)]
#fnho=DataFrame(fnho)
#fnho$trans_gene_name=nhgr$gene
#fnho$bins=G.Granges[fnho$queryHits]
#fnhos=sapply(split(fnho, fnho$queryHits), nrow)
#
#top40=(names(sort(fnhos, decreasing=T)))[1:40]
#fnhos[top40]

#xx=do.call('rbind', split(fnho, fnho$queryHits)[(names(sort(fnhos, decreasing=T)))[1:40]])
#write.table(do.call('rbind', split(fnho, fnho$queryHits)[(names(sort(fnhos, decreasing=T)))[1:40]]), file='~/not_hot_hotspots.txt', sep='\t', quote=F)

#yy=fnho[match(top40, fnho$queryHits),]
#write.table(yy, file='~/not_hot_hotspots_bins_sorted.txt', sep='\t', quote=F)
#queryHits(fnho)

phmax=rep(NA,1000)
for(i in 1:1000) {
    pnotHot=marker.GR[sample(1:10337, length(notHot), replace=T)]
    pnotHot_genome=buildOverlapVector(G.Granges, pnotHot, '') 
    phmax[i]=max(pnotHot_genome)

}
par(mfrow=c(2,1))
hist(notHot_genome, breaks=1000, xlim=c(0,175), main='(not hotspot) trans counts per 10kb genome bin', ylim=c(0,100))
abline(v=max(phmax), col='red')
hist(pnotHot_genome, breaks=1000, xlim=c(0,175), main='randomized positions (not hotspot) trans counts per 10kb genome bin')


# overlap of not hotspot trans with cis genes -------------------------
par(mfrow=c(2,2))
cisTransOver=buildOverlapVector(cis.GRanges,  notHot, cis.GRanges$ORF) 
#rle(sort(cisTransOver))   1900
pzero=rep(NA,1000)
pzm=rep(NA,1000)
pz5=rep(NA,1000)
tctos=list()
for(i in 1:1000) {
    print(i)
    pcis=marker.GR[sample(1:10337, length(notHot), replace=T)]
    pcisTransOver=buildOverlapVector(cis.GRanges,  pcis, cis.GRanges$ORF) 
    pzero[i]=rle(sort(pcisTransOver))$lengths[1]
    pzm[i]=max(pcisTransOver)
    pz5[i]=sum(pcisTransOver>5)

    tcto=table(pcisTransOver)
    tcto[7]=sum(tcto[7:length(tcto)])
    tcto=tcto[1:7]
    tctos[[as.character(i)]]=tcto
}
#x11()
#
# figure histogram overlap -------------------
pdf(file='/data/eQTL/FiguresForPublication/histogram_figure.pdf', width=6, height=6,useDingbats=F )

atcto=data.frame(do.call('rbind', tctos))
stripchart(atcto, vertical=T, pch=19, cex=.5, ylim=c(0,3000), method='jitter', 
           jitter=.25, at=c(0:6), group.names=c(0:6), xlim=c(-.25,6.25), col='#00000022',
           xlab='(trans eQTL peaks not in hotspots) per local eQTL bin', ylab='number of local eQTL bins'
           )
tcto=table(cisTransOver)
tcto[7]=sum(tcto[7:length(tcto)])
tcto=tcto[1:7]
points(c(0:6), tcto, pch='-', cex=3, col='red')
dev.off()
#--------------------------------------------

expected.tcto=colMeans(atcto)
expected.tcto.sd=apply(atcto,2, sd)
etcto=data.frame(tcto=tcto-expected.tcto, tcto.sd=expected.tcto.sd)
ggplot(etcto, aes(x=tcto.cisTransOver, y=tcto.Freq))+geom_bar(stat="identity", fill='grey')+
    geom_errorbar(width=0, aes(ymin=tcto.Freq-2*tcto.sd, ymax=tcto.Freq+2*tcto.sd))+theme_bw()+
    xlab('Non-hotspot trans eQTL per local eQTL bin')+ylab('(observed - expected) number of local eQTL bins')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(file='/data/eQTL/FiguresForPublication/histogram_figure6_updated.pdf')

#barplot(tcto, xlab='(not hotspot) trans eQTL peaks per local eQTL bin', add=T)
#
#
#par(mfrow=c(2,1))
#hist(cisTransOver, breaks=10000, xlim=c(0,10), main='(not hotspot) trans counts per cis bin', ylim=c(0,3000))
#hist(pcisTransOver, breaks=10000, xlim=c(0,10),main='randomized positions (not hotspot) trans counts per cis bin', ylim=c(0,3000))
#
#
#
## overlap of not hotspot trans with not-cis genes ------------------------
#notcisTransOver=buildOverlapVector(notcis.GRanges,  notHot, notcis.GRanges$ORF) 
#rle(sort(notcisTransOver))
##rle(sort(notcisTransOver))   2141
#pnzero=rep(NA,1000)
#for(i in 1:1000) {
#pnotcis=marker.GR[sample(1:10337, length(notHot), replace=T)]
#pnotcisTransOver=buildOverlapVector(notcis.GRanges,  pnotcis, notcis.GRanges$ORF) 
#pnzero[i]=rle(sort(pnotcisTransOver))$lengths[1]
#}
##x11()
##par(mfrow=c(2,1))
#hist(notcisTransOver, breaks=10000,xlim=c(0,10),main='(not hotspot) trans counts per not cis bin', ylim=c(0,3000))
#hist(pnotcisTransOver, breaks=10000, xlim=c(0,10),main='randomized positions (not hotspot) trans counts per not cis bin', ylim=c(0,3000))
##----------------------------------------------------------------------------
#
#
#
#cisAll=buildOverlapVector(cis.GRanges, trans.GRanges, '') 
#pzeroCis=rep(NA,1000)
#for(i in 1:1000) {
#print(i)
#pcisALL=marker.GR[sample(1:11530, length(trans.GRanges), replace=T)]
#pcisALL=buildOverlapVector(cis.GRanges,  pcisALL, '') 
#pzeroCis[i]=rle(sort(pcisALL))$lengths[1]
#}
#
#
#
#
#
#
#
##G.Granges=G.Granges[1:1200]
#
#
#
#cisTransOver=buildOverlapVector(cis.GRanges2kb, trans.GRanges, cis.GRanges2kb$ORF) 
#cisNotTransOver=buildOverlapVector(cis.GRanges2kb, notHot, cis.GRanges2kb$ORF) 
#
#plot(cis.GRanges2kb$gcoord, cisTransOver, type='h', ylim=c(0,250))
#
##notCISgenes=colnames(t.tpm.matrix)[!(colnames(t.tpm.matrix) %in% acis$gene)]
##notcis.GRanges2kb=gene.GR[match(notCISgenes, gene.GR$ORF)]+100
##notcisTransOver=buildOverlapVector(notcis.GRanges2kb, trans.GRanges, notcis.GRanges2kb$ORF) 
#
##plot(notcis.GRanges2kb$gcoord, notcisTransOver, type='h', ylim=c(0,250))
#
#par(mfrow=c(2,1))
#xx=hist(cisTransOver,breaks=100000, xlim=c(0,50), ylim=c(0,1200), main='trans peak counts at CIS genes')
#hist(notcisTransOver,breaks=100000, xlim=c(0,50), ylim=c(0,1200), main='trans peak counts at not CIS genes')
#
#
#ptrans=marker.GR[sample(1:11530, length(trans.GRanges), replace=T)]
#permFullOver=buildOverlapVector(ptrans,  trans.GRanges, '')
#
#x11()
#
#ptrans=marker.GR[sample(1:11530,2950, replace=T)]
#    #ptrans=marker.GR[sample(1:11530, sum(cisTransOver<40), replace=T)]
#permSubOverNH=buildOverlapVector(cis.GRanges2kb, ptrans, '')
#
#odist=cisTransOver[cisTransOver<20]
#rodist=rle(sort(odist))
#rodist.n=rodist$lengths/(sum(rodist$lengths))
#plot(rodist.n, type='b')
#       points(pdist.n, type='b', col='blue')
#ddist=rep(0,200)
#for(t in 2:200 ) {
#    ptrans=marker.GR[sample(1:11530,sum(cisTransOver[cisTransOver<t]), replace=T)]
#    #ptrans=marker.GR[sample(1:11530, sum(cisTransOver<40), replace=T)]
#    permSubOver=buildOverlapVector(cis.GRanges2kb, ptrans, '')
#    pdist=rle(sort(permSubOver))
#    pdist.n=pdist$lengths/(sum(pdist$lengths))
#
#    ddist[t]=sum((pdist.n-rodist.n)^2)
#
#    #df=approxfun(density(permSubOver, bw='nrd0'))
#    #ddist[t]=sum(na.omit(-log(df(cisTransOver[cisTransOver<20]))))
#    #print()
#}
#X=0:20
#points(X+1, .28*exp(-X/2.5), type='b', col='green', lty=2)
#
#r1=rle(sort(cisTransOver))
#plot(r1$values, cumsum(r1$lengths*r1$values))
#
#plot(log10(r1$values), (cumsum(r1$lengths*r1$values)))
#
#
#
#yy=hist(permSubOver,breaks=100000, xlim=c(0,50)) #, ylim=c(0,1200))
#
#
#fb=findOverlaps(cis.GRanges2kb, boot.GRanges)
#notBootCIS.GRanges=cis.GRanges2kb[-unique(queryHits(fb)),]
##ptrans=marker.GR[sample(1:11530, sum(cisTransOver<40), replace=T)]
#notBootOver=buildOverlapVector(notBootCIS.GRanges, trans.GRanges, '')
#hist(notBootOver,breaks=100000)
#
#plot(density(cisTransOver, n=512, bw='nrd0', adjust=1), xlim=c(0,50), ylim=c(0,.2))
#points(density(permSubOver, bw='nrd0'), col='blue', type='l')
#
#, xlim=c(0,50)) #, ylim=c(0,1200))
#
#
#notBootCisGenome=buildOverlapVector(G.Granges, notBootCIS.GRanges, '')
#
#
## gene +/- 1kb
#
#notCISgenes=colnames(t.tpm.matrix)[!(colnames(t.tpm.matrix) %in% acis$gene)]
#cis.GRanges2kb=gene.GR[match(notCISgenes, gene.GR$ORF)]+1000
#
#fb=findOverlaps(cis.GRanges2kb, boot.GRanges)
#
#fo=findOverlaps(cis.GRanges2kb, trans.GRanges)
#fo=findOverlaps(notBootCIS.GRanges, trans.GRanges)
#
#
#plot(cis.GRanges2kb$gcoord, bvec, type='h', ylim=c(0,250))
#abline(v=gcoord.key, lty=2, col='blue')
#hist(bvec, breaks=1000,  main='trans peaks per cis gene +/- 1kb interval' )
#hist(bvec, breaks=100000, xlim=c(0,50),  main='trans peaks per cis gene +/- 1kb interval' )
#
#
#abline(v= marker.GR[match(rownames(B), marker.GR$mname)]$gcoord, col='red')
#
##permutations
##ptrans=marker.GR[sample(1:11530, 2607, replace=T)]+1000
#ptrans=marker.GR[sample(1:11530, length(trans.GRanges), replace=T)]
#ptrans=marker.GR[sample(1:11530, sum(bvec[bvec<40]), replace=T)]
#ptrans=marker.GR[sample(1:11530, sum(bvec<40), replace=T)]
#ptrans=marker.GR[sample(1:11530, 2950, replace=T)]
#
#fp=findOverlaps(G.Granges, ptrans)
#sfp=split(subjectHits(fp), queryHits(fp))
#fblength=sapply(sfp, length)
#fbvec=rep(0, length(G.Granges))
#fbvec[as.numeric(names(fblength))]=fblength
#x11()
#hist(fbvec, breaks=1000, xlim=c(0,15))
#quantile(fbvec,.99)
#
#
#
#plot(cis.GRanges2kb$gcoord, bvec, type='h', ylim=c(0,500))
#abline(v=gcoord.key, lty=2, col='blue')
#abline(v= marker.GR[match(rownames(B), marker.GR$mname)]$gcoord, col='red')
#
#plot(cis.GRanges2kb$gcoord, fbvec, type='h')
#abline(v=gcoord.key, lty=2, col='blue')
#
#
#
#cis.GRanges=GRanges(seqnames=acis$chr, ranges=IRanges(start= as.numeric(sapply(strsplit(acis$CI.l, ':|_'), function(x)x[2])),
#                                                       end = as.numeric(sapply(strsplit(acis$CI.r, ':|_'), function(x)x[2]))), gene=acis$gene)
#
#length(unique(queryHits(findOverlaps(cis.GRanges, boot.GRanges))))
#library(topGO)
#myGene2GO.full.table=read.delim('/data/eQTL/reference/go_annotations_sgd.txt',  sep='\t', header=F, stringsAsFactors=F)
#
## data frame of gene then GO
#myGene2GO=cbind( sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,5])
#myGene2GO=na.omit(myGene2GO)
#
#SYS2ORF=unique(cbind(sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,3]))
#SYS2ORF.key=list()
#SYS2ORF.key[SYS2ORF[,1]]=SYS2ORF[,2]
#gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
#names(gene2GOList) = unique(myGene2GO[,1])
#
#
#
#geneList=cis.GRanges2kb$ORF[bvec==0]
#testGenes= factor(0+(colnames(t.tpm.matrix) %in% geneList)) # colnames(t.tpm.matrix)[(which(vcA[,1]>.7))]) #totest )
#names(testGenes)=colnames(t.tpm.matrix)
#testGeneSet=names(testGenes[testGenes==1])
#
# for( thisOntology in c("BP", "MF", "CC") ) {
#        GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
#        GOresult = runTest(GOData, algorithm="classic", statistic="fisher")
#        pdf(file=paste0('/data/eQTL/RData/GO/', filename.clean(colnames(ppa)[i]),'_', thisOntology, '.pdf'), width=25, height=25)
#        plotGOToTree(GOData, GOresult, sigThres = 0.00005)
#        dev.off()
#        gt=GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult)))
#        
#        genesinterms=genesInTerm(GOData, gt[,1])
#        genes.enriched.list=lapply(genesinterms, function(x) x[x%in%names(testGenes[testGenes==1])])
#        genes.enriched.list.simple=lapply(genes.enriched.list, function(x) as.character(SYS2ORF.key[x]))
#        gt$Genes=as.vector(sapply(genes.enriched.list.simple, paste, collapse=','))
#        gt$GenesSystematic=   as.vector(sapply(genes.enriched.list, paste, collapse=','))
#        write.table(gt, file=paste0('/data/eQTL/RData/GO/', filename.clean(colnames(ppa)[i]),'_', thisOntology, '.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
#    }
#
