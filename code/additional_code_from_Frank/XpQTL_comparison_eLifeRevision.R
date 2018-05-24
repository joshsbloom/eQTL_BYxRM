library("GenomicRanges")
library(qvalue)
library(DescTools)

# read gene locations
geneInfo = read.table("ensemblGenes_ensembl83_160307_MOD.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)
rownames(geneInfo) <- geneInfo[,"geneID"]
allNames <- geneInfo[, "geneName"]
names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]

allNamesInv <- names(allNames)
names(allNamesInv) <- allNames



getGcoords = function ( chr , pos, spacing=0, sgd.table="sacCer3ChromLenghts.txt" ) {
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

load("R_XpQTLDistantsacCer3")
load("R_allPeaksODPad_161213.RData")
load("peakModel.RData")

# attach the proper coefficients
naturalLMCoeff <- sapply(rownames(allPeaksODPad), function(x){
    thisPM <- peakModel[[allPeaksODPad[x, "gene"]]]
    thisPM[x, "lm.coeff"]
})
allPeaksODPad <- cbind(allPeaksODPad, naturalLMCoeff)


# load the tpm matrix to see which genes were tested
serverPrefix = ""
load(paste0(serverPrefix, "log2_t.tpm.matrix.RData", sep=""))


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

# restrict both sets to the genes that are present in both
distantGRanges = distantGRanges[distantGRanges$gene %in% XpQTL$gene]
# now 1059 eQTL

XpQTL = XpQTL[XpQTL$gene %in% distantGRanges$gene]
# 1024/1025 XpQTL




######################
# local eQTL

localPeaks = allPeaksODPad[allPeaksODPad$cis,]

localXQTL = read.table("localXpQTL.txt", stringsAsFactors=FALSE, sep="\t", head=TRUE)
nrow(localXQTL)
# 41 genes tested
length(which(localXQTL$localLOD > 3))
# 21

length(which(unique(localPeaks$gene) %in% localXQTL$Gene))
# 27 genes have a local eQTL

localPeaks[localPeaks$gene %in% localXQTL$Gene[localXQTL$Gene %in% localPeaks$gene],1:7]
# YBR067C had a second cis eQTL, now it's only one with stranded data

# now restrict local genes to the ones that can be tested
localPeaks = localPeaks[localPeaks$gene %in% localXQTL$Gene[localXQTL$Gene %in% localPeaks$gene],]
rownames(localPeaks) = localPeaks$gene


nrow(localXQTL[localXQTL$Gene %in% localPeaks$gene & localXQTL$localLOD > 3,])
# 16 out of 21, i.e. 76%

nrow(localXQTL[localXQTL$localLOD > 3 & !localXQTL$Gene %in% localPeaks$gene,])
# 5

fisher.test(cbind(
    c(nrow(localXQTL[localXQTL$Gene %in% localPeaks$gene & localXQTL$localLOD >= 3,]), nrow(localXQTL[!(localXQTL$Gene %in% localPeaks$gene) & localXQTL$localLOD >= 3,])),
    c(nrow(localXQTL[localXQTL$Gene %in% localPeaks$gene & localXQTL$localLOD < 3,]), nrow(localXQTL[!(localXQTL$Gene %in% localPeaks$gene) & localXQTL$localLOD < 3,]))
))
# 16    11
# 5     9
# p=0.2, odds = 2.6

# some association here, but not perfect
# in fact, not even > chance

pdf("pQTL_localLODcomparisonAtPeaks.pdf", width=8, height=4)
par(mfrow=c(1,2))
plot(localXQTL$localLOD[localXQTL$Gene %in% localPeaks$gene], localPeaks[localXQTL$Gene[localXQTL$Gene %in% localPeaks$gene], "LOD"], xlab="X-pQTL LOD", ylab="eQTL LOD")
plot(-localXQTL$alleleFrequencyDifference[localXQTL$Gene %in% localPeaks$gene], localPeaks[localXQTL$Gene[localXQTL$Gene %in% localPeaks$gene], "r"], xlab="X-pQTL delta AF", ylab="eQTL r")
abline(h=0, col="grey")
abline(v=0, col="grey")
dev.off()

cor.test(-localXQTL$alleleFrequencyDifference[localXQTL$Gene %in% localPeaks$gene], localPeaks[localXQTL$Gene[localXQTL$Gene %in% localPeaks$gene], "r"], method="s")
# rho=0.53, p = 0.005

# see at end for analysis NOT at respective peaks, but at the gene start (very similar to this, but slightly different markers)


########################
# moving away from peak overlap
# for every peak position, what is its r/LOD/AF in other experiment?

# these are the exact files used in the XpQTL paper
load("R_data4AcrossAnalyses.RData")

# make the names just the ORF so they can be accessed
names(countData) <- sapply(names(countData), function(x){ strsplit(x, "_")[[1]][1]})

# load the lifted-over & fixed file (s. previous scripts where this came from):
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
XatE <- lapply(distantGRanges_localRemoved, function(x){
    thisChr <- strsplit(as.character(seqnames(x)), "_")[[1]][1]
    thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
    #thisRange <- GRanges(seqnames = thisChr, ranges=ranges(x))
    thisRange <- GRanges(seqnames = thisChr, ranges=IRanges(start=x$peakPos, end=x$peakPos))
    thisNearest <- nearest(thisRange, sacCer3SNPsForXQTLGR)
    list(
        cbind(countData[[thisGene]][[1]][thisNearest,], countData[[thisGene]][[2]][thisNearest,]), # marker itself
        cbind(colSums(countData[[thisGene]][[1]][(thisNearest - 10):(thisNearest + 10),]), colSums(countData[[thisGene]][[2]][(thisNearest - 10):(thisNearest + 10),])), # sum of closest 21 markers
        cbind(apply(countData[[thisGene]][[1]][(thisNearest - 10):(thisNearest + 10),], 2, median), apply(countData[[thisGene]][[2]][(thisNearest - 10):(thisNearest + 10),], 2, median)) # median of closest 21 markers
    )
})

# do a G-test on each marker:

pValsInXAtEPositions <- t(sapply(XatE, function(x){
    c(GTest(x[[1]])$p.value, GTest(x[[2]])$p.value, GTest(x[[3]])$p.value)
}))

pdf("pQTL_pValsInXAtEPositions.pdf")
hist(pValsInXAtEPositions[,1], breaks=100)
hist(pValsInXAtEPositions[,2], breaks=100)
hist(pValsInXAtEPositions[,3], breaks=100)
dev.off()

1 - qvalue(pValsInXAtEPositions[,1])$pi0
# 0.1366334
# hmm, that is low
# this is of course extremely conservative because each SNP is so noisy

1 - qvalue(pValsInXAtEPositions[,2])$pi0
# 0.65545
# this SEEMS more believable, but is suboptimal because of the arbitrary SNP window
# let's look at a few random sets:

sapply(1:5, function(i){
    print(i)
    theseCounts <- lapply(distantGRanges_localRemoved, function(x){
        thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
        # pick a random marker:
        thisRange <- sample(sacCer3SNPsForXQTLGR, 1)
        thisNearest <- nearest(thisRange, sacCer3SNPsForXQTLGR)
        thisDat <- cbind(c(10, 10), c(10, 10))
        # sometimes fails when we hit a random marker close a chr end, therefore the 'try'
        try({thisDat <- cbind(colSums(countData[[thisGene]][[1]][(thisNearest - 10):(thisNearest + 10),]), colSums(countData[[thisGene]][[2]][(thisNearest - 10):(thisNearest + 10),]))}) # sum of closest 21 markers
        print(thisDat)
        thisDat
    })
    1 - qvalue(sapply(theseCounts, function(x){
        GTest(x)$p.value
    }))$pi0
})
#  0.4615135 0.3664628 0.5167904 0.5228282 0.4398202
# so there is inflation by summing the counts
# but the real data still stands out
# in fact, it stands above the randoms by ~0.2, similar to below
# not all of this is inflation – more that the random background is so high because so many QTL (optimistic interpretation)


# restrict to "strong" eQTL (VE >= 0.035)
1 - qvalue(pValsInXAtEPositions[,2][distantGRanges_localRemoved$VE >= 0.035])$pi0
# only 204 due to requirement to not be on gene chromosome
# 0.9039334 !!!

# make gene set exactly the same as for count-based analyses in paper:
XatE_2 <- lapply(distantGRanges, function(x){
    thisChr <- strsplit(as.character(seqnames(x)), "_")[[1]][1]
    thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
    #thisRange <- GRanges(seqnames = thisChr, ranges=ranges(x))
    thisRange <- GRanges(seqnames = thisChr, ranges=IRanges(start=x$peakPos, end=x$peakPos))
    thisNearest <- nearest(thisRange, sacCer3SNPsForXQTLGR)
    list(
    cbind(countData[[thisGene]][[1]][thisNearest,], countData[[thisGene]][[2]][thisNearest,]), # marker itself
    cbind(colSums(countData[[thisGene]][[1]][(thisNearest - 10):(thisNearest + 10),]), colSums(countData[[thisGene]][[2]][(thisNearest - 10):(thisNearest + 10),])), # sum of closest 21 markers
    cbind(apply(countData[[thisGene]][[1]][(thisNearest - 10):(thisNearest + 10),], 2, median), apply(countData[[thisGene]][[2]][(thisNearest - 10):(thisNearest + 10),], 2, median)) # median of closest 21 markers
    )
})

# do a G-test on each marker:

pValsInXAtEPositions_2 <- t(sapply(XatE_2, function(x){
    c(GTest(x[[1]])$p.value, GTest(x[[2]])$p.value, GTest(x[[3]])$p.value)
}))

1 - qvalue(pValsInXAtEPositions_2[,2][distantGRanges$VE >= 0.035])$pi0
# 236 genes
# 0.9198847
# wow!

# random?
sapply(1:5, function(i){
    print(i)
    theseCounts <- lapply(distantGRanges[distantGRanges$VE >= 0.035], function(x){
        thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
        # pick a random marker:
        thisRange <- sample(sacCer3SNPsForXQTLGR, 1)
        thisNearest <- nearest(thisRange, sacCer3SNPsForXQTLGR)
        thisDat <- cbind(c(10, 10), c(10, 10))
        # sometimes fails when we hit a random marker close a chr end, therefore the 'try'
        try({thisDat <- cbind(colSums(countData[[thisGene]][[1]][(thisNearest - 10):(thisNearest + 10),]), colSums(countData[[thisGene]][[2]][(thisNearest - 10):(thisNearest + 10),]))}) # sum of closest 21 markers
        #print(thisDat)
        thisDat
    })
    1 - qvalue(sapply(theseCounts, function(x){
        GTest(x)$p.value
    }))$pi0
})
# 0.6805868 0.3266363 0.5294731 0.4921313 0.2767041
# yes, real data is more than expected


#######################
# and the same, the other direction
# at all distant pQTL positions, pull tpm & genotypes, and compute pi1

load("gbatch_fact.RData")
load("ODcov.RData")
load("marker.LD.RData")
load("gdata_42k.RData")

phenoBatchOD = apply(t.tpm.matrix, 2, function(y){scale(lm(y ~ gbatch.fact + OD.cov)$res)})


# if want to use this, need lookup function:
getLDBlock <- function(markerDat, thisMarker, LDCut = 1){
    thisChr <- strsplit(thisMarker, ":")[[1]][1]
    allLDmarkers <- colnames(markerDat[[thisChr]])[markerDat[[thisChr]][thisMarker,] >= LDCut]
    return(c(allLDmarkers[1], allLDmarkers[length(allLDmarkers)]))
}

getClosestMarker <- function(markerDat, thisChr, thisPosition){
    thesePositions <- sapply(colnames(markerDat[[thisChr]]), function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})
    #print(head(thesePositions))
    #print(which.min(abs(thesePositions - thisPosition)))
    colnames(markerDat[[thisChr]])[which.min(abs(thesePositions - thisPosition))][1]
}

getPositionOfLinkedMarker <- function(markerDat, thisChr, thisPosition, LDCut = 1, leftOrRight = "left"){
    theseLDMarkers <- getLDBlock(markerDat, getClosestMarker(markerDat, thisChr, thisPosition), LDCut)
    if (leftOrRight == "left"){return(as.integer(strsplit(strsplit(theseLDMarkers[1], ":")[[1]][2], "_")[[1]][1]))}
    if (leftOrRight == "right"){return(as.integer(strsplit(strsplit(theseLDMarkers[length(theseLDMarkers)], ":")[[1]][2], "_")[[1]][1]))}
}


# remove any chance of working with a local QTL by accident - remove any on same chr
XpQTL_localRemoved <- XpQTL[sapply(seqnames(XpQTL), function(x){strsplit(x, "_")[[1]][1]}) != geneInfo[sapply(seqnames(XpQTL), function(x){strsplit(x, "_")[[1]][2]}), "chr"]]
# 946 out of 1024

pValsInEAtXPositions <- sapply(XpQTL_localRemoved, function(x){
    thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
    thisChr <- strsplit(as.character(seqnames(x)), "_")[[1]][1]
    thisMarker <- getClosestMarker(marker.LD, thisChr, x$peakPos)
    t.test(phenoBatchOD[,thisGene][gdata[,thisMarker] == 1], phenoBatchOD[,thisGene][gdata[,thisMarker] == -1])$p.value
})

pdf("pQTL_pValsInEAtXPositions.pdf")
hist(pValsInEAtXPositions, breaks=100)
dev.off()

1 - qvalue(pValsInEAtXPositions)$pi0
# pi1 = 0.46777
# at least 47% of pQTL overlap an eQTL

# better than at random marker?
sapply(1:5, function(i){
    print(i)
    randomPositionPs <- sapply(XpQTL_localRemoved, function(x){
        thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
        #thisChr <- strsplit(as.character(seqnames(x)), "_")[[1]][1]
        # pick a random marker:
        thisMarker <- sample(colnames(gdata), 1)
        t.test(phenoBatchOD[,thisGene][gdata[,thisMarker] == 1], phenoBatchOD[,thisGene][gdata[,thisMarker] == -1])$p.value
    })
    1 - qvalue(randomPositionPs)$pi0
})
# 0.2856617 0.2119434 0.3321604 0.3632650 0.2654992
# these are worse than observed



# That was for all genes, what about for "strong" pQTL?
# defined as LOD > 15 in paper

1 - qvalue(pValsInEAtXPositions[XpQTL_localRemoved$LOD >= 15])$pi0
# 0.6151075
# there's a slight discrepancy here:
# this is 207 pQTL total
# for counting overlaps as in paper, we used "XpQTL" from the distant experiment without explicitly chucking QTL on the gene chromosome, as done here for XpQTL_localRemoved

# to be perfectly consistent:
pValsInEAtXPositions_2 <- sapply(XpQTL, function(x){
    thisGene <- strsplit(as.character(seqnames(x)), "_")[[1]][2]
    thisChr <- strsplit(as.character(seqnames(x)), "_")[[1]][1]
    thisMarker <- getClosestMarker(marker.LD, thisChr, x$peakPos)
    t.test(phenoBatchOD[,thisGene][gdata[,thisMarker] == 1], phenoBatchOD[,thisGene][gdata[,thisMarker] == -1])$p.value
})

1 - qvalue(pValsInEAtXPositions_2[XpQTL$LOD >= 15])$pi0
# 0.6346427
# OK, even got us 2%



#################
# now that we have the full eQTL result, can redo the local analyses:

# let's use the gene start as the marker to search
# note gene strand!
geneMarker <- geneInfo[,"start"]
geneMarker[geneInfo[,"strand"] == "-"] <- geneInfo[,"end"][geneInfo[,"strand"] == "-"]
names(geneMarker) <- rownames(geneInfo)

localE <- apply(localXQTL, 1, function(x){
    thisChr <- geneInfo[x[1], "chr"]
    thisGene <- x[1]
    thisRange <- GRanges(seqnames = thisChr, ranges=IRanges(start=geneMarker[thisGene], end=geneMarker[thisGene]))
    scanoneLODS.OD[[1]][thisGene,][nearest(thisRange, eQTLSNPGR)]
})

pdf("pQTL_localLODcomparisonAll.pdf")
plot(localE, -localXQTL[,"alleleFrequencyDifference"], ylab="X-pQTL allele frequency difference", xlab="eQTL r", main="effects at gene position", pch=19, col="#00000088")#, cex=log10(localXQTL[,"localLOD"]))
abline(h=0, col="grey", lty=2, lwd=2)
abline(v=0, col="grey", lty=2, lwd=2)
dev.off()

cor.test(-localXQTL[,"alleleFrequencyDifference"], localE, method="s")
# rho = 0.57, p = 0.00011

# the two worst off-diagonals are YAL060W (top left) and YBR067C (bottom right)
# their local eQTL look legit when plotting as above - although YAL060W may be close to the chrI hotspot



