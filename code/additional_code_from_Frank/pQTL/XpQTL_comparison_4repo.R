library("GenomicRanges")

# read gene locations
geneInfo = read.table("ensemblGenes_ensembl83_160307_MOD.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)
rownames(geneInfo) <- geneInfo[,"geneID"]
allNames <- geneInfo[, "geneName"]
names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]

allNamesInv <- names(allNames)
names(allNamesInv) <- allNames



getGcoords = function ( chr , pos, spacing=0, sgd.table="sacCer3ChromLengths.txt" ) {
    offind = as.vector(cumsum(read.table(sgd.table, header=FALSE, sep="\t")[,2] + spacing))
    offind=    offind[-length(offind)]
    offind= c(0, offind)
    names(offind) = as.character(read.table(sgd.table, header=FALSE, sep="\t")[,1])
    chr.off=as.numeric(sapply(chr, function(x) {offind[[x]]}))
    return(chr.off+pos)
}
sepBetweenChr <- 1e5


################
# compare to X-pQTL

# load the table of distant peaks from our paper
# previously processed pre BoG
# local XQTL were excluded then
load("R_XpQTLDistantsacCer3.RData")

# load the current peak table, with padding, redefined cis definitions
load("R_allPeaksODPad_161213.RData")

# and the correct lm.coeff
load("peakModel.RData")

# these aren't strictly needed here because we won't work with the fold changes per se
# but still better to have them
# attach the proper coefficients
# note that the code has to be different than in the ASE analyses because
naturalLMCoeff <- sapply(rownames(allPeaksODPad), function(x){
    thisPM <- peakModel[[allPeaksODPad[x, "gene"]]]
    thisPM[x, "lm.coeff"]
})

allPeaksODPad <- cbind(allPeaksODPad, naturalLMCoeff)


# load the tpm matrix to see which genes were tested
load("log2_t.tpm.matrix.RData")


# split into distant and local peaks (analyze those separately)
# turn each table into a GRanges object

eQTLGRangesPad <- GRanges(seqnames = paste0(allPeaksODPad[,"chr"], "_", allPeaksODPad$gene), ranges = IRanges(start=sapply(allPeaksODPad[,"CI.l"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])}), end=sapply(allPeaksODPad[,"CI.r"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})))
eQTLGRangesPad$gene <- allPeaksODPad$gene
eQTLGRangesPad$cis <- allPeaksODPad$cis
eQTLGRangesPad$r <- allPeaksODPad$r
eQTLGRangesPad$LOD <- allPeaksODPad$LOD
eQTLGRangesPad$VE <- allPeaksODPad$"var.exp"
eQTLGRangesPad$LM <- allPeaksODPad$"naturalLMCoeff"
eQTLGRangesPad$peakPos <- sapply(allPeaksODPad[,"pmarker"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})

distantGRanges = eQTLGRangesPad[!eQTLGRangesPad$cis]


XpQTL = GRanges(seqnames = paste0(XpQTLDistantsacCer3[,"chromosome"], "_", XpQTLDistantsacCer3[,"gene"]), ranges = IRanges(start=as.numeric(XpQTLDistantsacCer3[,"X2LODIntervalLeft"]), end=as.numeric(XpQTLDistantsacCer3[,"X2LODIntervalRight"])), peakPos = as.numeric(XpQTLDistantsacCer3[,"peakPosition"]), r =  -as.numeric(XpQTLDistantsacCer3[,"alleleFrequencyDifference"]), gene =  XpQTLDistantsacCer3[,"gene"], LOD = XpQTLDistantsacCer3[,"LOD"])

length(unique(XpQTL$gene))
# 155
length(unique(XpQTL$gene)[unique(XpQTL$gene) %in% distantGRanges$gene])
# 154 genes
# all but one genes tested in XQTL are present (i.e. have a distant peak!) in the eQTL data
# which one is missing?
unique(XpQTL$gene)[!unique(XpQTL$gene) %in% distantGRanges$gene]
# does this guy have any eQTL?
"YFR056C" %in% allPeaksODPad$gene
# FALSE
# YFR056C (a dubious ORF) has no eQTL at all - it had one in pre-stranded data, that eQTL was probably due to sense strand reads from an overlapping gene

"YFR056C" %in% colnames(t.tpm.matrix)
# FALSE => mapping wasn't even attempted for this gene

XpQTL[XpQTL$gene == "YFR056C"]
# it did have one X-pQTL with LOD=8.8


# restrict both sets to the genes that are present in both
distantGRanges = distantGRanges[distantGRanges$gene %in% XpQTL$gene]
# now 1059 eQTL

XpQTL = XpQTL[XpQTL$gene %in% distantGRanges$gene]
# 1024/1025 XpQTL




# how many pQTL overlap an eQTL?
summary(overlapsAny(XpQTL, distantGRanges))
# 314 out of 1024 = 30.7%
summary(overlapsAny(distantGRanges, XpQTL))
# 321 / 1059 = 30.3%

# the numbers of matching QTL do not agree.
# probably because an QTL can overlap multiple in the other set?
# do the math:
table(countOverlaps(XpQTL, distantGRanges))
#0   1   2
#710 297  17

table(countOverlaps(distantGRanges, XpQTL))
#0   1   2   3
#738 312   8   1

distantGRanges[which(countOverlaps(distantGRanges, XpQTL) == 3)]
# chrIV_YGR209C [162498, 1080944], notice that this spans >900kb, LOD is 3.2

XpQTL[subjectHits(findOverlaps(distantGRanges[which(countOverlaps(distantGRanges, XpQTL) == 3)], XpQTL))]
# yep, three XpQTL there


# how many eQTL are there with effects such that expected detection power is beyond:
# 99% power
length(distantGRanges[abs(distantGRanges$VE) >= 0.035])
# 236 out of 1059
summary(overlapsAny(distantGRanges[abs(distantGRanges$VE) >= 0.035], XpQTL))
# 111 / 236 => 47%


# converse is difficult because have no power calcs for X-pQTL. Use high LOD instead
summary(overlapsAny(XpQTL[as.numeric(XpQTL$LOD) >=50], distantGRanges))
# 28/38 => 74%

# can lower this threshold to make more comparable to the number of "strong" eQTL (236)
summary(overlapsAny(XpQTL[as.numeric(XpQTL$LOD) >=15], distantGRanges))
# LOD 20
# 82/144 => 57%
# LOD 15
# 108/218 => 50%


# find the actual overlaps
pOverlapE = findOverlaps(XpQTL, distantGRanges)
eOverlapP = findOverlaps(distantGRanges, XpQTL)

# determine how often the effect goes the same direction
summary(sign(XpQTL[queryHits(pOverlapE)]$r) == sign(distantGRanges[subjectHits(pOverlapE)]$r))
# 254 / 331 (redundancy due to multiple-eQTL-overlap!). 77%


# plot the effects for the QTL that DO overlap

outlierThres = 0.15

pdf("pQTL_distantLODComparisonAtPeaks.pdf")

outliersEtoX <- sign(distantGRanges[queryHits(eOverlapP)]$r) != sign(XpQTL[subjectHits(eOverlapP)]$r) & abs(distantGRanges[queryHits(eOverlapP)]$r) > outlierThres & abs(XpQTL[subjectHits(eOverlapP)]$r) > outlierThres
plot(distantGRanges[queryHits(eOverlapP)]$r, XpQTL[subjectHits(eOverlapP)]$r, ylab="X-pQTL delta AF", xlab="eQTL effect", main="effects at shared QTL", col="#00000022", pch=19)
points(distantGRanges[queryHits(eOverlapP)]$r[outliersEtoX], XpQTL[subjectHits(eOverlapP)]$r[outliersEtoX], col="blue")

text(distantGRanges[queryHits(eOverlapP)]$r[outliersEtoX], XpQTL[subjectHits(eOverlapP)]$r[outliersEtoX], pos=3, cex=0.5, labels=paste(sapply(as.character(seqnames(distantGRanges[queryHits(eOverlapP)]))[outliersEtoX], function(x){strsplit(x, "_")[[1]][1]}), distantGRanges[queryHits(eOverlapP)][outliersEtoX]$peakPos, distantGRanges[queryHits(eOverlapP)][outliersEtoX]$gene, sep="_"))

abline(h=0, col="grey", lty=2)
abline(v=0, col="grey", lty=2)

dev.off()

as.data.frame(distantGRanges[queryHits(eOverlapP)][outliersEtoX])
# 19 total
# 5 HAP1, 11 MKT1
# 3 others
# OK, so disagreements mostly on HAP1 & MKT1

as.data.frame(XpQTL[subjectHits(eOverlapP)][outliersEtoX])

cor.test(XpQTL[queryHits(pOverlapE)]$r, distantGRanges[subjectHits(pOverlapE)]$r, method="s")
# r: rho=0.50



########################
# moving away from peak overlap
# for every peak position, what is its r/LOD/AF in other epxeriment?

# XpQTL results for all genes and all genome positions
load("R_diffsSmoothedAll.RData")

# make the names just the ORF so they can be accessed
names(diffsSmoothedAll) <- sapply(names(diffsSmoothedAll), function(x){ strsplit(x, "_")[[1]][1]})


# load the lifted-over & fixed file of coordinates linking sacCer2 (Albert 2014) to sacCer3 (now):
sacCer3SNPsForXQTL <- as.character(read.table("allXQTLsacCer2coordsAfterLiftover.txt", head=FALSE, sep="\t", stringsAsFactors=FALSE)[,1])
# turn them into an GRanges object so that we can easily find the closest to any given position
sacCer3SNPsForXQTLDF <- data.frame(t(sapply(sacCer3SNPsForXQTL, function(x){
    c(strsplit(x, ":")[[1]][1], strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])
})), stringsAsFactors=FALSE)

sacCer3SNPsForXQTLGR <- GRanges(seqnames = sacCer3SNPsForXQTLDF[,1], ranges = IRanges(start=as.numeric(sacCer3SNPsForXQTLDF[,2]), end=as.numeric(sacCer3SNPsForXQTLDF[,2])))

# restrict to eQTL on different chromosomes from gene location
distantGRanges_localRemoved <- distantGRanges[sapply(seqnames(distantGRanges), function(x){strsplit(x, "_")[[1]][1]}) != geneInfo[sapply(seqnames(distantGRanges), function(x){strsplit(x, "_")[[1]][2]}), "chr"]]
# down to 992 eQTL, from 1059

# now we should be able to pull the closest marker from the XQTL result for any position
XatE <- sapply(distantGRanges_localRemoved, function(x){
    thisChr <- strsplit(as.character(seqnames(x)), "_")[[1]][1]
    thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
    thisRange <- GRanges(seqnames = thisChr, ranges=IRanges(start=x$peakPos, end=x$peakPos))
    diffsSmoothedAll[[thisGene]][nearest(thisRange, sacCer3SNPsForXQTLGR)]
})

# OK, now plot this:
# NOTE THE INVERTED SIGN IN THE XQTL EFFECTS
# needed here because the diffsSmoothedAll were not yet inverted
# NOT needed below because "XpQTL" was inverted at construction

# make tables
XAtETable <- as.data.frame(distantGRanges_localRemoved)
# NOTE THE MINUS SIGN
XAtETable <- cbind(XAtETable, -XatE)
XAtETable <- cbind(XAtETable, overlapsAny(distantGRanges_localRemoved, XpQTL), overlapsAny(distantGRanges_localRemoved+10000, XpQTL+10000))
names(XAtETable)[(ncol(XAtETable)-2):ncol(XAtETable)] <- c("XpQTL_allele_difference", "XQTLOverlap", "XQTLOverlap_10kPad")
XAtETable <- XAtETable[order(XAtETable$r, decreasing=TRUE),]

write.table(XAtETable, file="pQTL_XAtETable_170331.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


plotCols = rep("#00000022", nrow(XAtETable))
plotColsRings <- rep("#FFFFFF00", nrow(XAtETable))
plotCols[XAtETable$XQTLOverlap_10kPad] <- "#FF000044"
plotColsRings[XAtETable$XQTLOverlap] <- "red"

outliers <- XAtETable[!XAtETable$XQTLOverlap_10kPad & abs(XAtETable$r) > 0.35,]

pdf("pQTL_XatE_Distant.pdf")
plot(XAtETable$r, XAtETable$XpQTL_allele_difference, xlab="eQTL effect", ylab="X-pQTL allele frequency difference", pch=19, col=plotCols, main="from mRNA QTL to protein effect", cex=log10(XAtETable$LOD))
points(XAtETable$r, XAtETable$XpQTL_allele_difference, col=plotColsRings, cex=log10(XAtETable$LOD))

text(outliers$r, outliers$XpQTL_allele_difference, pos=3, cex=0.5, labels=paste(sapply(as.character(outliers$seqnames), function(x){strsplit(x, "_")[[1]][1]}), outliers$peakPos, allNames[outliers$gene], sep="_"))
points(outliers$r, outliers$XpQTL_allele_difference, col="blue", cex=log10(outliers$LOD))

abline(h=0, col="grey", lty=2, lwd=2)
abline(v=0, col="grey", lty=2, lwd=2)
dev.off()

cor.test(XAtETable$r, XAtETable$XpQTL_allele_difference, method="s")
# rho=0.40, p < 2.2e-16 with all locals removed

# the outliers are heavily concentrated at MKT1, and somewhat on HAP1



#######################
# and the same, the other direction

load("scanoneLODS_OD_stranded.RData")


eQTLSNPDF <- data.frame(t(sapply(colnames(scanoneLODS.OD[[1]]), function(x){
    c(strsplit(x, ":")[[1]][1], strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])
})), stringsAsFactors=FALSE)

eQTLSNPGR <- GRanges(seqnames = eQTLSNPDF[,1], ranges = IRanges(start=as.numeric(eQTLSNPDF[,2]), end=as.numeric(eQTLSNPDF[,2])))


XpQTL_localRemoved <- XpQTL[sapply(seqnames(XpQTL), function(x){strsplit(x, "_")[[1]][1]}) != geneInfo[sapply(seqnames(XpQTL), function(x){strsplit(x, "_")[[1]][2]}), "chr"]]
# 946 out of 1024

# now we should be able to pull the closest marker from the XQTL result for any position
EatX <- sapply(XpQTL_localRemoved, function(x){
    thisChr <- strsplit(as.character(seqnames(x)), "_")[[1]][1]
    thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
    #thisRange <- GRanges(seqnames = thisChr, ranges=ranges(x))
    thisRange <- GRanges(seqnames = thisChr, ranges=IRanges(start=x$peakPos, end=x$peakPos))
    scanoneLODS.OD[[1]][thisGene,][nearest(thisRange, eQTLSNPGR)]
})


# MAKE TABLES
EAtXTable <- as.data.frame(XpQTL_localRemoved)
# NOTE THAT NO MINUS NEEDED BY CONSTRUCTION OF XpQTL OBJECT
EAtXTable <- cbind(EAtXTable, EatX)
EAtXTable <- cbind(EAtXTable, overlapsAny(XpQTL_localRemoved, distantGRanges), overlapsAny(XpQTL_localRemoved+10000, distantGRanges+10000))
names(EAtXTable)[(ncol(EAtXTable)-2):ncol(EAtXTable)] <- c("eQTL_r", "eQTLOverlap", "eQTLOverlap_10kPad")
names(EAtXTable)[names(EAtXTable) == "r"] <- c("XpQTL_allele_difference")
EAtXTable <- EAtXTable[order(EAtXTable$XpQTL_allele_difference, decreasing=TRUE),]

#write.table(EAtXTable, file="pQTL_EAtXTable_170331.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



# want to color red if there was peak overlap
plotCols = rep("#00000022", nrow(EAtXTable))
plotColsRings <- rep("#FFFFFF00", nrow(EAtXTable))
plotCols[EAtXTable$eQTLOverlap_10kPad] <- "#FF000044"
plotColsRings[EAtXTable$eQTLOverlap] <- "red"

outliers <- EAtXTable[!EAtXTable$eQTLOverlap_10kPad & abs(EAtXTable$XpQTL_allele_difference) > 0.3,]

pdf("pQTL_EatX_Distant.pdf")
plot(EAtXTable$eQTL_r, EAtXTable$XpQTL_allele_difference, xlab="eQTL effect", ylab="X-pQTL allele frequency difference", pch=19, col=plotCols, main="from protein QTL to mRNA effect", cex=log10(EAtXTable$LOD))
points(EAtXTable$eQTL_r, EAtXTable$XpQTL_allele_difference, col=plotColsRings, cex=log10(EAtXTable$LOD))
points(outliers$eQTL_r, outliers$XpQTL_allele_difference, col="blue", cex=log10(outliers$LOD))
text(outliers$eQTL_r, outliers$XpQTL_allele_difference, pos=3, cex=0.5, labels=paste(sapply(as.character(outliers$seqnames), function(x){strsplit(x, "_")[[1]][1]}), outliers$peakPos, allNames[outliers$gene], sep="_"))

abline(h=0, col="grey", lty=2, lwd=2)
abline(v=0, col="grey", lty=2, lwd=2)
dev.off()

cor.test(EAtXTable$eQTL_r, EAtXTable$XpQTL_allele_difference, method="s")
# rho=0.37, p < 2.2e-16


################################
# hotspot plots for eQTL and X-pQTL

chromosomeDividers <- c(0, 230218, 1043402, 1360022, 2891955, 3468829, 3738990, 4829930, 5392573, 5832461, 6578212, 7245028, 8323205, 9247636, 10031969, 11123260, 12071326)



pdf("hotspotLocationsHistogram.pdf", width=9, height=5)
hist(allPeaksODPad$marker.gcoord[!allPeaksODPad$cis], breaks=600, main="", xaxt="n", xlab="genome position (chromosome name)", col="black", ylab="number of genes influenced")
for(i in 2:(length(chromosomeDividers)-1)){
    abline(v=chromosomeDividers[i], col="grey", lty=2)
}
#axis(1, at=chromosomeDividers[-17], labels=paste0("chr", as.roman(1:16)))
axis(1, at=sapply(1:16, function(i){chromosomeDividers[i] + ((chromosomeDividers[i+1] - chromosomeDividers[i])/2)}), labels=as.roman(1:16), tick=FALSE)
dev.off()

doubleHistPlot <- function(xVals, yVals1, yVals2, thisYLim = NA, thisMain=""){
    if(is.na(thisYLim)){thisYLim <- c(-max(yVals2), max(yVals1))}
    
    plot(xVals, yVals1, type="h", main=thisMain, xaxt="n", xlab="genome position (chromosome name)", col="black", ylab="fraction of genes influenced", yaxt="n",
    ylim=thisYLim)
    points(xVals, -yVals2, type="h")
    for(i in 2:(length(chromosomeDividers)-1)){
        abline(v=chromosomeDividers[i], col="grey", lty=2)
    }
    axis(1, at=sapply(1:16, function(i){chromosomeDividers[i] + ((chromosomeDividers[i+1] - chromosomeDividers[i])/2)}), labels=as.roman(1:16), tick=FALSE)
    axis(2, at = seq(from=-1, to=1, by=0.2), labels=abs(seq(from=-1, to=1, by=0.2)))
    legend("topleft", legend="eQTL", box.lty=0)
    legend("bottomleft", legend="pQTL", box.lty=0)

}

pdf("pQTL_hotspotLocationsHistograms.pdf", width=9, height=5)

eQTLHist <- hist(allPeaksODPad$marker.gcoord[!allPeaksODPad$cis], breaks=600, plot=FALSE)
pQTLHist <- hist(getGcoords(XpQTLDistantsacCer3$chromosome, as.numeric(XpQTLDistantsacCer3$peakPosition)), breaks=eQTLHist$breaks, plot=FALSE)

doubleHistPlot(eQTLHist$mids, eQTLHist$counts / length(unique(allPeaksODPad$gene)), pQTLHist$counts / length(unique(XpQTLDistantsacCer3$gene)), thisMain = "all hotspots, all genes")

eQTLHist <- hist(allPeaksODPad$marker.gcoord[(!allPeaksODPad$cis) & allPeaksODPad$chr != "chrIII"], breaks=600, plot=FALSE)
pQTLHist <- hist(getGcoords(XpQTLDistantsacCer3$chromosome, as.numeric(XpQTLDistantsacCer3$peakPosition))[XpQTLDistantsacCer3$chromosome != "chrIII"], breaks=eQTLHist$breaks, plot=FALSE)

doubleHistPlot(eQTLHist$mids, eQTLHist$counts / length(unique(allPeaksODPad$gene)), pQTLHist$counts / length(unique(XpQTLDistantsacCer3$gene)), thisMain = "without chrIII, all genes")
doubleHistPlot(eQTLHist$mids, eQTLHist$counts / length(unique(allPeaksODPad$gene)), pQTLHist$counts / length(unique(XpQTLDistantsacCer3$gene)), thisMain = "without chrIII, all genes, zoomed", thisYLim = c(-0.2, 0.2))

eQTLHist_pQTLOnly <- hist(allPeaksODPad$marker.gcoord[(!allPeaksODPad$cis) & allPeaksODPad$chr != "chrIII" & allPeaksODPad$gene %in% XpQTLDistantsacCer3$gene], breaks=eQTLHist$breaks, plot=FALSE)
pQTLHist_pQTLOnly <- hist(getGcoords(XpQTLDistantsacCer3$chromosome, as.numeric(XpQTLDistantsacCer3$peakPosition))[XpQTLDistantsacCer3$chromosome != "chrIII"], breaks=eQTLHist$breaks, plot=FALSE)

doubleHistPlot(eQTLHist_pQTLOnly$mids, eQTLHist_pQTLOnly$counts / length(which(unique(allPeaksODPad$gene) %in% XpQTLDistantsacCer3$gene)), pQTLHist_pQTLOnly$counts / length(unique(XpQTLDistantsacCer3$gene)), thisMain = "without chrIII, X-pQTL genes only")

dev.off()


