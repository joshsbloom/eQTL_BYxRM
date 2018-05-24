library(igraph) # network plotting
library(GenomicRanges)
library(topGO)
library(Rgraphviz)
library(Gviz)
library(qvalue)
library(ggplot2)
library(reshape2)


##########################################
# load objects and annotation

load("gbatch_fact.RData")
load("ODcov.RData")
load("log2_t.tpm.matrix.RData")
# 1012 x 5720

# gdata_42k has 42052 markers, many of which are very highly correlated
load("gdata_42k.RData")

hotspot.boot.intervals <- read.table("MV.12417_mod.txt", sep="\t", head=TRUE, stringsAsFactors=FALSE)
rownames(hotspot.boot.intervals) <- hotspot.boot.intervals[,1]
hotspot.boot.intervals <- hotspot.boot.intervals[,2:6]

load("mv.bootstrap.pos.RData")
hotspot.list <- list(rownames(hotspot.boot.intervals), mv.boots, hotspot.boot.intervals)
names(hotspot.list[[2]]) <- hotspot.list[[1]]
rownames(hotspot.list[[3]]) <- hotspot.list[[1]]

# hotspot effect coefficients:
load("B.Forward.RData")
load("B.cis.Forward.RData")
# B is 102 x 5629


geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]

allNamesInv <- names(allNames)
names(allNamesInv) <- allNames


# big object, careful on laptop
# list with one netry per chromosome
# contains pairwise correlations between markers on each chromosome
load("marker.LD.RData")

# if want to use marker.LD, need lookup function:
# gives boundaries of markers that within some LD threshold of each other (used here for markers in perfect LD)
getLDBlock <- function(markerDat, thisMarker, LDCut = 1){
    thisChr <- strsplit(thisMarker, ":")[[1]][1]
    allLDmarkers <- colnames(markerDat[[thisChr]])[markerDat[[thisChr]][thisMarker,] >= LDCut]
    return(c(allLDmarkers[1], allLDmarkers[length(allLDmarkers)]))
}
# example:
#getLDBlock(marker.LD, "chrI:48890_C/T", 1)

# also need a lookup function to get the closest marker to some position
# because the boot intervals are defined just by position
getClosestMarker <- function(markerDat, thisChr, thisPosition){
    thesePositions <- sapply(colnames(markerDat[[thisChr]]), function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})
    colnames(markerDat[[thisChr]])[which.min(abs(thesePositions - thisPosition))][1]
}
# example
#getClosestMarker(marker.LD, "chrI", 153969)

# and also a convenience function weaving the above two together
getPositionOfLinkedMarker <- function(markerDat, thisChr, thisPosition, LDCut = 1, leftOrRight = "left"){
    theseLDMarkers <- getLDBlock(markerDat, getClosestMarker(markerDat, thisChr, thisPosition), LDCut)
    if (leftOrRight == "left"){return(as.integer(strsplit(strsplit(theseLDMarkers[1], ":")[[1]][2], "_")[[1]][1]))}
    if (leftOrRight == "right"){return(as.integer(strsplit(strsplit(theseLDMarkers[length(theseLDMarkers)], ":")[[1]][2], "_")[[1]][1]))}
}
#getPositionOfLinkedMarker(marker.LD, "chrI", 153971)


# pad the hotspots with fully linked markers
leftHotspotFlank <- sapply(rownames(hotspot.boot.intervals), function(tH){getPositionOfLinkedMarker(marker.LD, strsplit(tH, ":")[[1]][1], hotspot.boot.intervals[tH,1], leftOrRight="left")})
rightHotspotFlank <- sapply(rownames(hotspot.boot.intervals), function(tH){getPositionOfLinkedMarker(marker.LD, strsplit(tH, ":")[[1]][1], hotspot.boot.intervals[tH,5], leftOrRight="right")})
hotspotChrs <- sapply(rownames(hotspot.boot.intervals), function(x){strsplit(x, ":")[[1]][1]})

# turn the hotspots into genomic ranges:
hotspotGRanges <- GRanges(seqnames = hotspotChrs, ranges = IRanges(start=leftHotspotFlank, end=rightHotspotFlank))
names(hotspotGRanges) <- hotspot.list[[1]]

# same for the genes, with padding to catch promoters and 3'UTRs
geneGRanges <- GRanges(seqnames = geneAnnotation[,"chr"], ranges = IRanges(start=geneAnnotation[,"start"], end=geneAnnotation[,"end"]), strand = geneAnnotation[,"strand"], gene= geneAnnotation[,"geneID"])
geneGRangesExtended = resize(geneGRanges, width(geneGRanges)+200)
geneGRangesExtended = resize(geneGRangesExtended, fix="end", width(geneGRangesExtended)+1000)
names(geneGRangesExtended) <- geneGRangesExtended$gene

# replace original object, only work with expanded definitions throughout
geneGRanges <- geneGRangesExtended



# polymorphism data
# note that the first line in the file is commented out with a '#' and therefore does not get read by default => do not need header=TRUE
# quote="" is important, otherwise it breaks halfway down
# this is because of the gene "IMP2'" (note that ' is part of the gene name) - removed those tick marks from the file to be safe
allVariants <- read.table("gdata_42k_VEP_wExtraMarkersAttached.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
allVariants <- cbind(allVariants, sapply(allVariants[,1], function(x){strsplit(x, ":")[[1]][1]}))
colnames(allVariants)[ncol(allVariants)] <- "chr"
allVariants <- cbind(allVariants, as.integer(sapply(allVariants[,1], function(x){strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1]})))
colnames(allVariants)[ncol(allVariants)] <- "pos"

# cannot use markers as row names since a few markers are there several times:
nrow(allVariants)
# 48254
length(unique(allVariants[,1]))
# 47754
# fewer than all, because some genes overlap so that the variant gets included in the table twice


# this object has eQTl flanking markers padded to include all markers in perfect LD
# and bases the local definition on overlapping with padded genes (incl 1kb promoter, 200 bp downstream)
load("R_allPeaksODPad_161213.RData")
eQTLGRangesPad <- GRanges(seqnames = allPeaksODPad[,"chr"], ranges = IRanges(start=sapply(allPeaksODPad[,"CI.l"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])}), end=sapply(allPeaksODPad[,"CI.r"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})))
eQTLGRangesPad$gene <- allPeaksODPad$gene
eQTLGRangesPad$cis <- allPeaksODPad$cis


# with the padded definition, multiple peaks can locally overlap a gene:
length(unique(allPeaksODPad$gene[allPeaksODPad$cis]))
# 2884
length(allPeaksODPad$gene[allPeaksODPad$cis])
# 2969
# note that there are multiple "cis" per some genes in the ASE analyses
# no need to fix here, is fixed in ASE analyses




############################################################
# various comparisons of effect sizes of distant vs local eQTL


summary(allPeaksODPad$"var.exp"[allPeaksODPad$cis])
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0004321 0.0209743 0.0497970 0.1110371 0.1335431 0.9231440
summary(allPeaksODPad$"var.exp"[!allPeaksODPad$cis])
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000001 0.0122704 0.0178880 0.0337849 0.0312528 0.9001627

t.test(allPeaksODPad$"var.exp"[allPeaksODPad$cis], allPeaksODPad$"var.exp"[!allPeaksODPad$cis])
# p < 2.2e-16; same for wilcox

#ratio of medians:
median(allPeaksODPad$"var.exp"[allPeaksODPad$cis])/median(allPeaksODPad$"var.exp"[!allPeaksODPad$cis])
# 2.78

# plots of these distributions - probably use Josh's versions instead in paper

# effects summed per gene
# use the 5643 genes that have an eQTL at all
varExpSummed <- t(sapply(unique(allPeaksODPad$gene), function(x){
    c(sum(allPeaksODPad$"var.exp"[allPeaksODPad$cis & allPeaksODPad$gene == x]),
    sum(allPeaksODPad$"var.exp"[(!allPeaksODPad$cis) & allPeaksODPad$gene == x]))
}))
colnames(varExpSummed) <- c("local", "distant")
summary(varExpSummed)
#Min.   :0.000000   Min.   :0.00000
#1st Qu.:0.000000   1st Qu.:0.09136
#Median :0.009331   Median :0.16238
#Mean   :0.058421   Mean   :0.20074
#3rd Qu.:0.054344   3rd Qu.:0.26641
#Max.   :0.923144   Max.   :0.90422

# paired ratio trans/cis:
summary(varExpSummed[,2] / varExpSummed[,1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.000   2.415  39.840     Inf     Inf     Inf

# plot this:
pdf("effectSizesLocalVDistant_scatterplot.pdf", height=5, width=5)
smoothScatter(varExpSummed, pch=19, col="#00000033", cex=0.5, xlab="Variance from local eQTL per gene", ylab="Variance from summed distant eQTL per gene", nrpoints=0, xaxs="i", yaxs="i", xlim=c(0,1), ylim=c(0,1))
abline(0, 1, col="grey", lty=2)
dev.off()

# with violins:
pdf("effectSizesLocalVDistant_violins.pdf", height=4, width=4)
plotDat <- melt(varExpSummed)
names(plotDat) <- c("gene", "cisOrTrans", "varExplained")
p <- ggplot(plotDat, aes(x=cisOrTrans, y=varExplained)) +
geom_violin(trim=TRUE, fill='#A4A4A4', color="black") +
labs(title="all genes", x="eQTL type", y = "Variance from summed eQTL") +
theme_minimal()
p

dev.off()


# the paired comparison is perhaps the fairer one, but what if we just compare the distributions:
t.test(varExpSummed[,1], varExpSummed[,2])
# p < 2.2e-16 of course
median(varExpSummed[,2]) / median(varExpSummed[,1])
# 17 times more trans when just comparing the medians
# but this comparison is not good since half of the genes have no local eQTL
# also applies to the paired ratio above


# what if we restrict to genes that HAVE a local eQTL
varExpSummedWithLocal <- t(sapply(unique(allPeaksODPad$gene[allPeaksODPad$cis]), function(x){
    c(sum(allPeaksODPad$"var.exp"[allPeaksODPad$cis & allPeaksODPad$gene == x]),
    sum(allPeaksODPad$"var.exp"[(!allPeaksODPad$cis) & allPeaksODPad$gene == x]))
}))
summary(varExpSummedWithLocal)
# Min.   :0.0004321   Min.   :0.00000
#1st Qu.:0.0227748   1st Qu.:0.08162
#Median :0.0529158   Median :0.14764
#Mean   :0.1143097   Mean   :0.17916
#3rd Qu.:0.1376453   3rd Qu.:0.24335
#Max.   :0.9231440   Max.   :0.88931

summary(varExpSummedWithLocal[,2] / varExpSummedWithLocal[,1])
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000    0.7981    2.5550    7.7090    7.3420 1793.0000

# excess of trans is robust to requiring a local eQTL to be present



# fairest comparison is likely for genes that have at least one local and one distant
localAndDistant <- sapply(unique(allPeaksODPad$gene), function(x){
   subPeaks <- allPeaksODPad[allPeaksODPad$gene == x,]
   length(which(subPeaks$cis)) > 0 & nrow(subPeaks) > length(which(subPeaks$cis))
})
summary(localAndDistant)
# TRUE for 2846 genes, false for 2797

varExpSummedWithBoth <- t(sapply(names(localAndDistant)[localAndDistant], function(x){
    c(sum(allPeaksODPad$"var.exp"[allPeaksODPad$cis & allPeaksODPad$gene == x]),
    sum(allPeaksODPad$"var.exp"[(!allPeaksODPad$cis) & allPeaksODPad$gene == x]))
}))
colnames(varExpSummedWithBoth) <- c("local", "distant")
summary(varExpSummedWithBoth)
#Min.   :0.0004321   Min.   :0.0004247
#1st Qu.:0.0228002   1st Qu.:0.0848599
#Median :0.0528947   Median :0.1498827
#Mean   :0.1137831   Mean   :0.1815511
#3rd Qu.:0.1373104   3rd Qu.:0.2455859
#Max.   :0.8993735   Max.   :0.8893136

summary(varExpSummedWithBoth[,2] / varExpSummedWithBoth[,1])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.0005    0.8425    2.6225    7.8124    7.4602 1792.5327

t.test(varExpSummedWithBoth[,1], varExpSummedWithBoth[,2], paired=TRUE)

# plot with violins:
pdf("effectSizesLocalVDistant_violins_cisAndTrans.pdf", height=5, width=5)
plotDat <- melt(varExpSummedWithBoth)
names(plotDat) <- c("gene", "cisOrTrans", "varExplained")
p <- ggplot(plotDat, aes(x=cisOrTrans, y=varExplained)) +
geom_violin(trim=TRUE, fill='#A4A4A4', color="black") +
labs(title="genes with at least one local and one distant eQTL", x="eQTL type", y = "summed variance explained") +
theme_bw()
p

dev.off()

# how many eQTL are on a different chromosome than their gene?
summary(apply(allPeaksODPad, 1, function(x){
    x["chr"] != geneAnnotation[x["gene"],"chr"]
}))
#   Mode   FALSE    TRUE    NA's
#logical    5094   31404       0

# how many genes have a trans eQTL?
hasTrans <- sapply(unique(allPeaksODPad$gene), function(x){
    subPeaks <- allPeaksODPad[allPeaksODPad$gene == x,]
    length(which(!subPeaks$cis)) > 0
})
summary(hasTrans)
# 38 FALSE, 5605 TRUE






#########################################
# build GRanges, find overlaps
# need to decide what definition we want to use

# now, the overlaps

genesInHotspots = findOverlaps(hotspotGRanges, geneGRanges, ignore.strand=TRUE)
# 755 genes in hotspots defined by 95% confidence intervals, and with perfect LD extension

# to make things easier on the eyes, make a list of overlapping genes for each hotspot, with clear names
genesInHotspotList <- lapply(1:length(hotspotGRanges), function(i){
    list(
    geneGRanges[subjectHits(genesInHotspots[queryHits(genesInHotspots) == i,])]$gene,
    allNames[geneGRanges[subjectHits(genesInHotspots[queryHits(genesInHotspots) == i,])]$gene])
})
names(genesInHotspotList) <- names(hotspotGRanges)

# convenient for seeing hotspot gene content:
zz <- file("genesInHotspotList.txt", open="wt")
sink(zz)
genesInHotspotList
sink()

#save(genesInHotspotList, file="R_genesInHotspotList_95PadWLD_170131.RData")
# load in future runs:
#load("R_genesInHotspotList_95PadWLD_170131.RData")


hotspotPeakMarkerGRanges <- GRanges(seqnames = hotspotChrs, ranges = IRanges(
    start=sapply(rownames(hotspot.boot.intervals), function(tH){getPositionOfLinkedMarker(marker.LD, strsplit(tH, ":")[[1]][1], as.integer(strsplit(strsplit(tH, ":")[[1]][2], "_")[[1]][1]), leftOrRight="left")}),
    end=sapply(rownames(hotspot.boot.intervals), function(tH){getPositionOfLinkedMarker(marker.LD, strsplit(tH, ":")[[1]][1], as.integer(strsplit(strsplit(tH, ":")[[1]][2], "_")[[1]][1]), leftOrRight="right")})
))

# genes that overlap the peak hotspot marker
hotspotPeakGeneHits <- nearest(hotspotPeakMarkerGRanges, geneGRangesExtended, select="all")
hotspotPeakGenes <- lapply(unique(queryHits(hotspotPeakGeneHits)), function(x){allNames[geneGRangesExtended[subjectHits(hotspotPeakGeneHits)[queryHits(hotspotPeakGeneHits) == x]]$gene]})
names(hotspotPeakGenes) <- rownames(hotspot.boot.intervals)

# find gene closest to boot peak
hotspotBootPeakGRanges <- GRanges(seqnames = hotspotChrs, ranges = IRanges(start=sapply(rownames(hotspot.boot.intervals), function(tH){getPositionOfLinkedMarker(marker.LD, strsplit(tH, ":")[[1]][1], hotspot.boot.intervals[tH,3], leftOrRight="left")}), end=sapply(rownames(hotspot.boot.intervals), function(tH){getPositionOfLinkedMarker(marker.LD, strsplit(tH, ":")[[1]][1], hotspot.boot.intervals[tH,3], leftOrRight="right")})))
names(hotspotBootPeakGRanges) <- rownames(hotspot.boot.intervals)
bootPeakGeneHits <- nearest(hotspotBootPeakGRanges, geneGRangesExtended, select="all")
bootPeakGenes <- lapply(unique(queryHits(bootPeakGeneHits)), function(x){allNames[geneGRangesExtended[subjectHits(bootPeakGeneHits)[queryHits(bootPeakGeneHits) == x]]$gene]})
names(bootPeakGenes) <- rownames(hotspot.boot.intervals)

length(unlist(bootPeakGenes))
# 186 genes



#################################
# hotspot region plots

# the columns in hotspot.boot.intervals are quantiles of the available bootstraps
# c(.025, .05, .95, .975)
# i.e., 95% confidence are the outermost columns

# what is the difference between 90 and 95% CI width
# i.e. is it worth 5% extra error for much smaller regions?
pdf("hotspotCI_arithmetic.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(hotspot.boot.intervals[,5] - hotspot.boot.intervals[,1], hotspot.boot.intervals[,4] - hotspot.boot.intervals[,2], xlab="95% confidence interval", ylab="90% confidence interval")
abline(0, 1, col="grey", lty=2)
hist((hotspot.boot.intervals[,5] - hotspot.boot.intervals[,1]) - (hotspot.boot.intervals[,4] - hotspot.boot.intervals[,2]), xlab="difference between 95% and 90% confidence interval", main="", breaks=40)
dev.off()
# the real 90-95 CIs are very similar to each other
# go with 95 throughout


# also want to add the distance to the end of the chromosome

chrLengthsFile <- read.table("sacCer3ChromLengths.txt", stringsAsFactors=FALSE, head=FALSE, sep="\t")
chrLengths <- chrLengthsFile[,2]
names(chrLengths) <- chrLengthsFile[,1]

# Josh now bootstraps within 40 markers either side of peak
load("boot.windows.RData")

plotBufferSize = 20e3

# generate a multi-page pdf
# one plot per region
pdf("bootstrapDistributionsWithGeneModels.pdf", width=10, height=5)
thisChr <- ""
gtrack <- GenomeAxisTrack(labelPos="above")
peakPad=0
# if we use the hotspot list order, we get the hotspots in the order they were called
for (i in hotspot.list[[1]]){
#for (i in hotspot.list[[1]][78]){

#thisHist <- hist(hotspot.list[[2]][[i]], breaks=100, plot=FALSE)
    thisHist <- hist(hotspot.list[[2]][i,], breaks=100, plot=FALSE)
    
    # only fetch from UCSC if new chromosome
    if (thisChr != strsplit(i, ":")[[1]][1]){
        thisChr = strsplit(i, ":")[[1]][1]
    
        # pull the gene annotation (for the whole chromosome) and the ideogram
        ensGenes = UcscTrack(genome="sacCer3", chromosome=thisChr, track ="ensGene", from=0, to=10e6, trackType="GeneRegionTrack", name="EnsGenes", rstarts="exonStarts", rends="exonEnds", gene="name", strand="strand", symbol="name2", id="name", transcript="name")
        displayPars(ensGenes) <- list(stacking="squish", showTitle=FALSE, shape="arrow", transcriptAnnotation="gene", fill="grey")
        identifier(ensGenes) <- allNames[sapply(identifier(ensGenes), function(x){strsplit(x, "\\.")[[1]][1]})]
        iTrack <- IdeogramTrack(genome="sacCer3", chromosome=thisChr)
    }
    
    thisAllPeak <- as.integer(strsplit(strsplit(i, ":")[[1]][2], "_")[[1]][1])

    leftBootBoundaryLeft <- getPositionOfLinkedMarker(marker.LD, thisChr, boot.windows[i, 1], leftOrRight="left")
    rightBootBoundaryLeft <- getPositionOfLinkedMarker(marker.LD, thisChr, boot.windows[i, 2], leftOrRight="left")
    leftBootBoundaryRight <- getPositionOfLinkedMarker(marker.LD, thisChr, boot.windows[i, 1], leftOrRight="right")
    rightBootBoundaryRight <- getPositionOfLinkedMarker(marker.LD, thisChr, boot.windows[i, 2], leftOrRight="right")
    
    thisFrom = max(c(-1000, min(c(thisAllPeak - plotBufferSize/2, leftBootBoundaryLeft-1000)))) # outer max to prevent negative plot cooordinates
    thisTo = min(c(max(c(thisAllPeak + plotBufferSize/2, rightBootBoundaryRight+1000)), chrLengths[thisChr] + 1000)) # outer min to truncate at chr end
    
    variantsToPlot = allVariants[allVariants[,"chr"] == thisChr & allVariants[,"pos"] >= thisFrom & allVariants[,"pos"] <= thisTo, ]
    variantCols <- rep("grey", nrow(variantsToPlot))
    variantCols[variantsToPlot[,5] == "HIGH"] <- "red"
    variantCols[variantsToPlot[,5] == "MODERATE"] <- "orange"
    variantTrack <- AnnotationTrack(genome="sacCer3", chromosome=thisChr, start=variantsToPlot[,"pos"]-1, end=variantsToPlot[,"pos"]+1)
    displayPars(variantTrack) <- list(showTitle=FALSE, background.title="white", col.line="transparent", fill=variantCols, lty=0, stacking="dense")
    
    ht90 <- HighlightTrack(trackList = ensGenes, genome = "sacCer3", chromosome = thisChr,
        start=c(
            getPositionOfLinkedMarker(marker.LD, thisChr, hotspot.boot.intervals[i,1], leftOrRight="left"),
            getPositionOfLinkedMarker(marker.LD, thisChr, hotspot.boot.intervals[i,2], leftOrRight="left"),
            getPositionOfLinkedMarker(marker.LD, thisChr, hotspot.boot.intervals[i,3]-peakPad, leftOrRight="left"),
            getPositionOfLinkedMarker(marker.LD, thisChr, thisAllPeak-peakPad, leftOrRight="left"),
            leftBootBoundaryLeft,
            rightBootBoundaryLeft
            ),
        end=c(
            getPositionOfLinkedMarker(marker.LD, thisChr, hotspot.boot.intervals[i,5], leftOrRight="right"),
            getPositionOfLinkedMarker(marker.LD, thisChr, hotspot.boot.intervals[i,4], leftOrRight="right"),
            getPositionOfLinkedMarker(marker.LD, thisChr, hotspot.boot.intervals[i,3]+peakPad, leftOrRight="right"),
            getPositionOfLinkedMarker(marker.LD, thisChr, thisAllPeak+peakPad, leftOrRight="right"),
            leftBootBoundaryRight,
            rightBootBoundaryRight
            )
    )
    displayPars(ht90) <- list(col=c("lightskyblue1", "lightskyblue",  "darkblue", "red", "grey", "grey"), fill=c("lightskyblue1", "lightskyblue", "darkblue", "red", "grey", "grey"), lty=c(0, 0, 1, 1, 1, 1), lwd=c(1,1, 2, 1, 1, 1))
    histTrack <- DataTrack(data = thisHist$counts, start = thisHist$breaks[-length(thisHist$breaks)], end = thisHist$breaks[-1], name="bootstrap frequency", chromosome=thisChr, genome="sacCer3")
    plotTracks(list(iTrack, histTrack, gtrack, variantTrack, ht90), chromosome = thisChr, from=thisFrom, to = thisTo, type="histogram", main=paste0(i, "; ", length(which(B[i,] != 0)), " genes affected"), cex.main=1.5, fontcolor.title="black", col.axis="black", background.title="white")
    
}
dev.off()



##########################################
# network plots

# with all genes and cis SELECTIVELY blocked out
load("R_thetaOutAllHotspotsAllGenes0.1_0.1_170415.RData")

# variable is called yySolve for historical reasons
yySolve <- thetaOutAllHotspotsAllGenes$Theta$yy

# want to give "local" genes a different shape if they have an eQTL here (these were protected from cis correction during network construction):
# they may represent spurious, non-mediating eQTL because there are so many local eQTL
load("R_genesInHotspotList_95PadWLD_170131.RData")
thres = 10^-(5.5)

# genes within the attempted bootstrap region
bootWindowGRanges = GRanges(seqnames = hotspotChrs, ranges = IRanges(start=boot.windows[,1], end=boot.windows[,2]))
names(bootWindowGRanges) <- rownames(boot.windows)

genesInHotspotsBootDist = findOverlaps(bootWindowGRanges, geneGRanges, ignore.strand=TRUE)
# 4720 genes in 95% with perfect LD extension
# a window this wide is a little silly, but still useful to know when a local eQTL is even remotely close

genesInHotspotsBootDistList <- lapply(1:length(bootWindowGRanges), function(i){
    geneGRanges[subjectHits(genesInHotspotsBootDist[queryHits(genesInHotspotsBootDist) == i,])]$gene
})
names(genesInHotspotsBootDistList) <- names(bootWindowGRanges)


# multi-page pdf, with one page per hotspot
maxGenesForGraph = 100
pdf(paste0("hotspotNetworks_BGeneSelection_top",maxGenesForGraph,"genes.pdf", sep=""), width=21, height=7)
par(mfrow=c(1, 3))
for (thisHotSpot in hotspot.list[[1]]){
    print(thisHotSpot)
    genesToSelect = list()
    # do separately for top total, up, down
    # select on B != zero
    # need to also consider B.cis so that cis effects don't get ignored
    # probably easiests to construct a vector from B and replace the cis effects with those from B.cis, but only for genes in the genesInHotspotsBootDistList window
    # using genesInHotspotsBootDistList (which is the entire 80 marker window used for bootstraps) adds a lot of local eQTL, and brekas up the trans networks
    # only add cis effects for genes in the confidence interval
    theseEffects <- B[thisHotSpot,]
    cisToBeReplaced <- names(theseEffects)[names(theseEffects) %in% genesInHotspotList[[thisHotSpot]][[1]] & theseEffects == 0]
    #cisToBeReplaced <- names(theseEffects)[names(theseEffects) %in% genesInHotspotsBootDistList[[thisHotSpot]]]
    if(length(cisToBeReplaced) > 0){
        theseEffects[cisToBeReplaced] <- sapply(cisToBeReplaced, function(x){B.cis[which(abs(B.cis[,x]) == max(abs(B.cis[,x])))[1], x]})
    }
    if (length(which(theseEffects > 0)) <= maxGenesForGraph){
        genesToSelect[[1]] <- names(which((theseEffects) > 0))
    }else{
        genesToSelect[[1]] <- names(sort(theseEffects, decreasing=TRUE))[1:maxGenesForGraph]
    }
    if (length(which(theseEffects < 0)) <= maxGenesForGraph){
        genesToSelect[[2]] <- names(which((theseEffects) < 0))
    }else{
        genesToSelect[[2]] <- names(sort(theseEffects, decreasing=FALSE))[1:maxGenesForGraph]
    }
    genesToSelect[[3]] <- c(genesToSelect[[1]], genesToSelect[[2]])
    
    
    plotMains = c(paste0("top ", length(genesToSelect[[1]]), " up genes"), paste0("top ", length(genesToSelect[[2]]), " down genes"), paste0("top ", length(genesToSelect[[3]]), " up and down genes"))
    
    for (i in 1:length(genesToSelect)){
        if (length(genesToSelect[[i]]) > 1){
            tryTest <- try({
                gHotSpot <- graph.adjacency(yySolve[genesToSelect[[i]][genesToSelect[[i]] %in% rownames(yySolve)], genesToSelect[[i]][genesToSelect[[i]] %in% rownames(yySolve)]], mode="undirected", diag=FALSE, weighted=TRUE)
                E(gHotSpot)$width <- abs(E(gHotSpot)$weight) * 100
                E(gHotSpot)$edge.color[E(gHotSpot)$weight > 0] <- "#FF000033"
                E(gHotSpot)$edge.color[E(gHotSpot)$weight < 0] <- "#0000FF33"
                gHotSpot <- igraph::delete.edges(gHotSpot, E(gHotSpot)[abs(weight) < thres ])
                vertexCol = rep("lightblue", length(V(gHotSpot)))
                vertexCol[sign(theseEffects[names(V(gHotSpot))]) < 0] <- "yellow"
                vertexFrameCol <- rep("black", length(V(gHotSpot)))
                vertexFrameCol[names(V(gHotSpot)) %in% genesInHotspotsBootDistList[[thisHotSpot]]] <- "orange"
                vertexFrameCol[names(V(gHotSpot)) %in% genesInHotspotList[[thisHotSpot]][[1]]] <- "red"
                plot(gHotSpot, vertex.color=vertexCol, vertex.frame.color=vertexFrameCol, vertex.size=10*abs(theseEffects[names(V(gHotSpot))]), edge.color=E(gHotSpot)$edge.color, vertex.label=allNames[names(V(gHotSpot))], vertex.label.family="sans", vertex.label.cex=0.8, vertex.label.dist=0.2, layout=layout_with_fr(gHotSpot), main=paste0(thisHotSpot, " ", plotMains[i]))
            })
            if(class(tryTest) == "try-error"){
                plot.new()
            }
        }
        else{
            plot(1,1, type="n")
        }
    }
}
dev.off()





###########################
# TFBS enrichment per hotspot
# NOTE that we had to do IMP2' -> IMP2
TFBSData <- read.table("regulationRelationshipsFromSGD_160501_MOD.txt", head=TRUE, sep="\t", quote="\"", stringsAsFactors=FALSE)
colnames(TFBSData) <- c("TFName", "TFSystematic", "targetName", "targetSystematic", "evidenceType", "evidenceCode", "condition", "construct", "strain", "evidenceObservation", "refID", "SGD")

maxGenesForGraph = 50
TFBSEnrichments <- lapply(hotspot.list[[1]], function(thisHotSpot){
    genesToSelect = list()
    # do separately for top total, up, down
    # select on lasso != zero
    if (length(which(B[thisHotSpot,] > 0)) <= maxGenesForGraph){
        genesToSelect[[1]] <- names(which((B[thisHotSpot,]) > 0))
    }else{
        genesToSelect[[1]] <- names(sort(B[thisHotSpot,], decreasing=TRUE))[1:maxGenesForGraph]
    }
    if (length(which(B[thisHotSpot,] < 0)) <= maxGenesForGraph){
        genesToSelect[[2]] <- names(which((B[thisHotSpot,]) < 0))
    }else{
        genesToSelect[[2]] <- names(sort(B[thisHotSpot,], decreasing=FALSE))[1:maxGenesForGraph]
    }
    genesToSelect[[3]] <- c(genesToSelect[[1]], genesToSelect[[2]])
    
    plotMains = c(paste0("top ", length(genesToSelect[[1]]), " up genes"), paste0("top ", length(genesToSelect[[2]]), " down genes"), paste0("top ", length(genesToSelect[[3]]), " up and down genes"))
    
    ret1 <- lapply(1:length(genesToSelect), function(i){
        if (length(genesToSelect[[i]]) > 1){
            #print(genesToSelect[[i]])
            t(sapply(unique(TFBSData[,"TFSystematic"]), function(thisTF){
                #print(thisTF)
                subTFMat = TFBSData[TFBSData[,"TFSystematic"] == thisTF,]
                #print(genesToSelect[[i]] %in% subTFMat[,"targetSystematic"])
                testMat = cbind(
                c(length(which(genesToSelect[[i]] %in% subTFMat[,"targetSystematic"])),
                length(which(!genesToSelect[[i]] %in% subTFMat[,"targetSystematic"]))),
                c(length(which(colnames(B)[which(!colnames(B) %in% genesToSelect[[i]])] %in% subTFMat[,"targetSystematic"])),
                length(which(!colnames(B)[which(!colnames(B) %in% genesToSelect[[i]])] %in% subTFMat[,"targetSystematic"])))
                )
                res <- fisher.test(testMat)
                res1 <- c(res$p, res$est, as.numeric(testMat))
                names(res1) <- c("p", "odds", "regulatedAndTarget", "regulatedAndNotTarget", "notRegulatedAndTarget", "notRegulatedAndNotTarget")
                res1
            }))
        }
    })
    names(ret1) <- c("up", "down", "both")
    ret1
})
names(TFBSEnrichments) <- hotspot.list[[1]]

#save(TFBSEnrichments, file=paste0("R_TFBSEnrichments_top", maxGenesForGraph,"Genes_170420.RData"))


# visualize the enrichments
# which enrichments should we believe?
# there are 189 TFs in the annotation file
# and 102 hotspots
# times 2 for up & down direction of effect
# 38556 tests => Bonferroni 1.296815e-06
# use that in the plots
# in addition, always label (but not color in) the most significant TFBS, even if it's not significant
pdf(paste0("TFBS_SGD_enrichments_scatterplot_top", maxGenesForGraph, "Genes.pdf"), width=12, height=6)
par(mfrow=c(1,2))
mainLabels = c("up", "down")
for (thisHS in hotspot.list[[1]]){
    #for (thisHS in "chrXII:657022_T/C"){
    for (i in 1:2){
        print(c(thisHS, i))
        plotCols = rep("#00000044", nrow(TFBSEnrichments[[thisHS]][[i]]))
        namesArg <- allNames[rownames(TFBSEnrichments[[thisHS]][[i]])]
        plotCols[TFBSEnrichments[[thisHS]][[i]][,1] < 1.3e-6] <- "red"
        namesArg[TFBSEnrichments[[thisHS]][[i]][,1] >= 1.3e-6 & TFBSEnrichments[[thisHS]][[i]][,1] != min(TFBSEnrichments[[thisHS]][[i]][,1])] <- ""
        toPlot1 <- TFBSEnrichments[[thisHS]][[i]][,2]
        toPlot1[is.infinite(toPlot1)] <- NA
        toPlot2 <- -log10(TFBSEnrichments[[thisHS]][[i]][,1])
        toPlot2[is.infinite(toPlot2)] <- NA
        plot(toPlot1, toPlot2, col=plotCols, xlab="odds", ylab="-log10(p)", main=paste0(thisHS, " ", mainLabels[i]), pch=19)
        text(toPlot1, toPlot2, labels=namesArg, pos=1)
    }
}
dev.off()





########################################
#make an overview table of the hotspots for the paper


# best 5 GO terms
# these are pulled from GO results computed on MSI, will need to altered on other locations
# see "RESULTS_FOLDER" in file path
best5GO <- lapply(rownames(hotspot.boot.intervals), function(x){
    thisFile <- gsub(":", "_", x)
    thisFile <- gsub("/", "_", thisFile)
    lapply(c("topGenes50", "topGenes100"), function(thisGeneNumber){
        lapply(c("BP", "MF", "CC"), function(thisGOCat){
            lapply(c("Up", "Down"), function(thisUpDown){
                thisFile <- paste0("RESULTS_FOLDER/topGO/", thisGeneNumber, "/", thisGOCat, "/GOGraph_hotspot_", thisFile, "_Fisher", thisUpDown, "_GOResultTable.txt", sep="")
                thisGOResult <- read.table(thisFile, stringsAsFactors=FALSE, head=TRUE, sep="\t", quote="\"")
                # use a within-hotspot Bonferroni for each category (i.e. lenient overall)
                thisGOResult <- thisGOResult[thisGOResult[,"pValue"] < (0.05 / nrow(thisGOResult)),]
                if (nrow(thisGOResult) > 5){thisGOResult <- thisGOResult[1:5,]}
                thisGOResult
            })
        })
    })
})
names(best5GO) <- rownames(hotspot.boot.intervals)

# flag hotspots at the end of their chromosomes
# make two flags: is either the boot or the peak marker the last marker in the linkage map?
# or are they only 5kb away?
# first is a strict chr end definition, 2nd is more strict in terms of throwing hotspots out later on
isEndMarker <- t(sapply(rownames(hotspot.boot.intervals), function(x){
    maxAllowedDist = 5e3
    firstMarker <- rownames(marker.LD[[hotspotChrs[x]]])[1]
    lastMarker <- rownames(marker.LD[[hotspotChrs[x]]])[nrow(marker.LD[[hotspotChrs[x]]])]
    ret <- (x == firstMarker) | (x == lastMarker)
    bootMarkerLeft <- getClosestMarker(marker.LD, as.character(seqnames(hotspotBootPeakGRanges[x])), start(hotspotBootPeakGRanges[x]))
    bootMarkerRight <- getClosestMarker(marker.LD, as.character(seqnames(hotspotBootPeakGRanges[x])), end(hotspotBootPeakGRanges[x]))
    ret <- ret | ((firstMarker == bootMarkerLeft) | (lastMarker == bootMarkerRight))
    thisMarkerPosition <- as.numeric(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])
    lastMarkerPosition <- as.numeric(strsplit(strsplit(lastMarker, ":")[[1]][2], "_")[[1]][1])
    firstMarkerPosition <- as.numeric(strsplit(strsplit(firstMarker, ":")[[1]][2], "_")[[1]][1])
    bootMarkerLeftPosition <- as.numeric(strsplit(strsplit(bootMarkerLeft, ":")[[1]][2], "_")[[1]][1])
    bootMarkerRightPosition <- as.numeric(strsplit(strsplit(bootMarkerRight, ":")[[1]][2], "_")[[1]][1])
    ret1 <- abs(thisMarkerPosition - lastMarkerPosition) < maxAllowedDist | abs(thisMarkerPosition - firstMarkerPosition) < maxAllowedDist | abs(bootMarkerLeftPosition - firstMarkerPosition) < maxAllowedDist | abs(bootMarkerRightPosition - lastMarkerPosition) < maxAllowedDist
    ret2 <- c(ret, ret1)
    names(ret2) <- c("isEndMarker", "isCloseToEnd")
    ret2
}))
#save(isEndMarker, file="R_isEndMarker_170214.RData")

tableForPaper <- data.frame(rownames(hotspot.boot.intervals),
    hotspotChrs,
    hotspot.boot.intervals[,3],
    leftHotspotFlank, rightHotspotFlank,
    rightHotspotFlank-leftHotspotFlank,
    sapply(genesInHotspotList, function(x){length(x[[1]])}),
    countOverlaps(hotspotPeakMarkerGRanges, eQTLGRangesPad),
    sapply(rownames(hotspot.boot.intervals), function(i){length(which(B[i,] != 0))}),
    isEndMarker,
    sapply(bootPeakGenes, function(x){paste0(x, collapse=";")}),
    sapply(hotspotPeakGenes, function(x){paste0(x, collapse=";")}),
    sapply(genesInHotspotList, function(x){paste0(x[[2]][x[[1]] %in% allPeaksODPad$gene[allPeaksODPad$cis]], collapse=";")}),
    sapply(genesInHotspotList, function(x){paste0(x[[2]][x[[2]] %in% allVariants[,6][allVariants[,5] == "HIGH"]], collapse=";")}),
    sapply(genesInHotspotList, function(x){paste0(x[[2]][x[[2]] %in% allVariants[,6][allVariants[,5] == "MODERATE"]], collapse=";")}),
    sapply(genesInHotspotList, function(x){paste0(x[[2]], collapse=";")}),
    sapply(best5GO, function(x){paste0(x[[2]][[1]][[1]][,1], collapse=";")}), sapply(best5GO, function(x){paste0(x[[2]][[1]][[1]][,2], collapse=";")}),
    sapply(best5GO, function(x){paste0(x[[2]][[1]][[2]][,1], collapse=";")}), sapply(best5GO, function(x){paste0(x[[2]][[1]][[2]][,2], collapse=";")}),
    sapply(best5GO, function(x){paste0(x[[2]][[2]][[1]][,1], collapse=";")}), sapply(best5GO, function(x){paste0(x[[2]][[2]][[1]][,2], collapse=";")}),
    sapply(best5GO, function(x){paste0(x[[2]][[2]][[2]][,1], collapse=";")}), sapply(best5GO, function(x){paste0(x[[2]][[2]][[2]][,2], collapse=";")}),
    sapply(best5GO, function(x){paste0(x[[2]][[3]][[1]][,1], collapse=";")}), sapply(best5GO, function(x){paste0(x[[2]][[3]][[1]][,2], collapse=";")}),
    sapply(best5GO, function(x){paste0(x[[2]][[3]][[2]][,1], collapse=";")}), sapply(best5GO, function(x){paste0(x[[2]][[3]][[2]][,2], collapse=";")}),
# make sure to have the version based on 100 genes loaded, to make consistent with the GO results
    sapply(rownames(hotspot.boot.intervals), function(x){paste0(allNames[rownames(TFBSEnrichments[[x]][["up"]])[TFBSEnrichments[[x]][["up"]][,1] < 1.3e-6]], collapse=";")}),
    sapply(rownames(hotspot.boot.intervals), function(x){allNames[rownames(TFBSEnrichments[[x]][["up"]])[order(TFBSEnrichments[[x]][["up"]][,1], decreasing=FALSE)][1]]}),
    sapply(rownames(hotspot.boot.intervals), function(x){paste0(allNames[rownames(TFBSEnrichments[[x]][["down"]])[TFBSEnrichments[[x]][["down"]][,1] < 1.3e-6]], collapse=";")}),
    sapply(rownames(hotspot.boot.intervals), function(x){allNames[rownames(TFBSEnrichments[[x]][["down"]])[order(TFBSEnrichments[[x]][["down"]][,1], decreasing=FALSE)][1]]})
, stringsAsFactors=FALSE)

colnames(tableForPaper) <- c("hotspotMarker", "chromosome", "bootstrapPeak", "bootstrapIntervalLeft", "bootstrapIntervalRight", "boostrapIntervalWidth", "numberGenesInHotspot", "numberEQTLInHotspot", "numberNonzeroEffects", "telomeric", "telomere5k", "bestGeneBootstraps", "bestGenePeak", "localeQTL", "genesWithHighImpactVariants", "genesWithModerateImpactVariants", "allGenesInInterval", "GOID_BP_upInRM", "GOTerms_BP_upInRM", "GOID_BP_downInRM", "GOTerms_BP_downInRM", "GOID_CC_upInRM", "GOTerms_CC_upInRM", "GOID_CC_downInRM", "GOTerms_CC_downInRM", "GOID_MF_upInRM", "GOTerms_MF_upInRM", "GOID_MF_downInRM", "GOTerms_MF_downInRM", "TFBS_sig_upInRM", "TFBS_best_upInRM", "TFBS_sig_downInRM", "TFBS_best_downInRM")

#write.table(tableForPaper, file="tableForPaper_topGenes100_170420.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")


#####################################
# distribution of number of genes in hotspot
pdf("numberGenesInHotspotBarplot.pdf", width=9, height=8)
par(mfrow=c(2,1))
barplot(table(sapply(genesInHotspotList, function(x){length(x[[1]])})), xlab="number of genes in hotspot", col="blue", ylab="number of hotspots", main="all hotspots")
barplot(table(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])})), xlab="number of genes in hotspot", col="blue", ylab="number of hotspots", main="non-telomeric hotspots")
dev.off()


# number of hotspots with < 5 genes
# not telomeres
length(which(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}) <= 3))
sum(sapply(genesInHotspotList[!isEndMarker[,2]][which(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}) <= 3)], function(x){length(x[[1]])}))
# 26 hotspots, 58 genes

# three hotspots with exactly one gene:
genesInHotspotList[!isEndMarker[,2]][which(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}) == 1)]

# how many genes total?
sum(sapply(genesInHotspotList, function(x){length(x[[1]])}))
# 755
# how many unique?
length(unique(unlist(lapply(genesInHotspotList, function(x){x[[1]]}))))
# 755 as well

sum(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}))
# 722




###################################
# GO ON THE HOTSPOT GENES

myGene2GO = read.table("gene_association_MOD_fromGO_160229.txt", stringsAsFactors = FALSE, sep="\t", header=FALSE)
# make that gene to GO list format they want:
gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
names(gene2GOList) = unique(myGene2GO[,1])

plotGOToTree <- function(GOdat, GOres, sigThres = 0.0005){
    # only plot if there are any significant GO terms (SEE ABOVE for "significance"; I am somewhat lenient here):
    # we need these extra lines because very small p-values are reported as a text string "< X", rather than a numeric
    toTest <- as.numeric(GenTable(GOdat, pVal = GOres)[1,6])
    if(is.na(toTest)){toTest <- 0.000000000000000000000000000001}
    if (toTest < sigThres){
        showSigOfNodes(GOdat, score(GOres), firstSigNodes = length(score(GOres)[score(GOres) < sigThres]), useInfo = "all")
    }else{
        plot(1,1)
        legend("topright", legend="no significant GO categories", cex=0.4, box.lty=0)
    }
}


# restrict to hotspots with three or fewer genes, and throw out telomeric hotspots
# can try different definitions of which genes to consider by changing the inclusion criteria here
testGenes <- factor(as.integer(geneAnnotation[,"geneID"] %in% unlist(genesInHotspotList[!isEndMarker[,2]][which(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}) <= 3)])))
names(testGenes) <- geneAnnotation[,"geneID"]

# CAREFUL WITH THE FILE NAMES

#pdf("topGO/hotspotGenes/hotspotGeneGO_max3Genes_noTelomeres.pdf", width=15, height=15)
for (thisOntology in c("BP", "MF", "CC")){
    GOData <- new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
    GOresult <- runTest(GOData, algorithm="classic", statistic="fisher")
    plotGOToTree(GOData, GOresult, sigThres = 0.05 / length(score(GOresult)))
    
    write.table(GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult))), file=paste0("topGO/hotspotGenes/hotspotGOresult_", thisOntology, "_max3Genes_noTelomeres.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}
dev.off()


## randomize these enrichment tests
# let's focus on the test with ≤3 genes / hotspot
# can pick random trios of neighboring genes and run on those
# could restrict to avoid chr ends - let's for now say we avoid the first/last 10 genes on each chromosome
# important to pick adjacent genes to account for grouping of genes by function
# then rerun the test, and record pValue for each GO category
# can base this off code written during BoG 2016:

# want a table that has the same order of GOs as the actual one:
BPResult = read.table("topGO/hotspotGenes/hotspotGOresult_BP_max3Genes_noTelomeres.txt", head=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
MFResult = read.table("topGO/hotspotGenes/hotspotGOresult_MF_max3Genes_noTelomeres.txt", head=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
CCResult = read.table("topGO/hotspotGenes/hotspotGOresult_CC_max3Genes_noTelomeres.txt", head=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

pdf("hotspotGeneGOPpValueHists.pdf", width=9, height=6)
par(mfrow=c(2,3))
hist(BPResult$result1, breaks=100, main="BP")
hist(MFResult$result1, breaks=100, main="MF")
hist(CCResult$result1, breaks=100, main="CC")
hist(BPResult$result1[BPResult$result1 < 1], breaks=100, main="BP")
hist(MFResult$result1[MFResult$result1 < 1], breaks=100, main="MF")
hist(CCResult$result1[CCResult$result1 < 1], breaks=100, main="CC")
dev.off()

hotspotGOResult <- list(BPResult, MFResult, CCResult)
for (i in 1:3){
    rownames(hotspotGOResult[[i]]) <- hotspotGOResult[[i]][,1]
}
hotspotGOStats <- lapply(hotspotGOResult, function(x){
    y <- x$result1
    ret <- rep(NA, 3)
    try({ret[1] <- 1-qvalue(y[y < 1])$pi0})
    ret[2] <- length(y[y < 0.05/nrow(x)])
    ret[3] <- length(y[y < 0.05])
    names(ret) <- c("pi1", "bonferroniTerms", "nominalTerms")
    ret
})
hotspotGOStats
# pi1s are 0.67; 0.75; 0
# OK, clearly some signal AFTER taking out the (extremely many) p==1


GOTermsInOrder <- list(BPResult$GO.ID, MFResult$GO.ID, CCResult$GO.ID)
GOCats <- c("BP", "MF", "CC")

# throw our mitochondrial genes (no hotspots there...)
geneAnno4Perms <- geneAnnotation[geneAnnotation$chr != "chrM",]
# exclude first/last 10 genes from each chromosome (crude way of avoiding telomeric regions)
last10GenesPerChr <- lapply(unique(geneAnno4Perms$chr), function(x){
    excludeGenes <- 10
    c(rep(TRUE, excludeGenes), rep(FALSE, length(which(geneAnno4Perms$chr == x)) - 2*excludeGenes), rep(TRUE, excludeGenes))
})
geneAnno4Perms <- geneAnno4Perms[!unlist(last10GenesPerChr),]

# the exact distribution of gene stretches we want to sample:
genesToSample <- sapply(genesInHotspotList[!isEndMarker[,2]][which(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}) <= 3)], function(y){length(y[[1]])})

nPerms = 1000
permGO <- lapply(1:length(GOCats), function(thisOntology){
    ret <- sapply(1:nPerms, function(j){
        print(paste0(GOCats[thisOntology], " randomization: ", j))
        
        # construct the random set
        # there is a small chance that this will select some genes multiple times
        # ignore for now, would become problematic if there were more & bigger regions we'd sample
        sampledGenes <- unlist(lapply(genesToSample, function(x){
            # pick one chr at random (note that this will weigh the chr by gene number because there is no unique(), which is fine)
            theseChrGenes <- geneAnno4Perms[geneAnno4Perms$chr == sample(geneAnno4Perms$chr, 1), "geneID"]
            # pick a gene, staying away far enough from the chr end to accomodate the number we want to include
            thisSampledAnchorGene <- sample(1:(length(theseChrGenes) - x), 1)
            # pick the anchor gene and the needed number of additional genes to its right
            theseChrGenes[thisSampledAnchorGene:(thisSampledAnchorGene + x - 1)]
        }))
        #print(length(sampledGenes))
        permGenes <- as.factor(as.integer(geneAnno4Perms$geneID %in% sampledGenes))
        names(permGenes) <- geneAnno4Perms$geneID
        #print(permGenes)
        
        permGOData <- new("topGOdata", ontology=GOCats[thisOntology], allGenes = permGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
        permGOresult <- runTest(permGOData, algorithm="classic", statistic="fisher")
        permTable <- GenTable(permGOData, permGOresult, topNodes = length(score(permGOresult)))
        rownames(permTable) <- permTable[,"GO.ID"]
        return(as.numeric(permTable[GOTermsInOrder[[thisOntology]], "result1"]))
    })
    rownames(ret) <- GOTermsInOrder[[thisOntology]]
    ret
})
names(permGO) <- GOCats
#save(permGO, file="topGO/hotspotGenes/R_permGO_170206.RData")

# get distributions of pi1 & other stats from the permutations
permGOStats <- lapply(permGO, function(x){
    t(apply(x, 2, function(y){
        ret <- rep(NA, 3)
        try({ret[1] <- 1-qvalue(y[y < 1])$pi0})
        ret[2] <- length(y[y < 0.05/nrow(x)])
        ret[3] <- length(y[y < 0.05])
        names(ret) <- c("pi1", "bonferroniTerms", "nominalTerms")
        ret
    }))
})
summary(permGOStats[[1]])

pdf("hotspotGeneGOPpermResults.pdf", width=10, height=10)
par(mfrow=c(3,3))
for (i in 1:3){
    for(j in colnames(permGOStats[[i]])){
        hist(permGOStats[[i]][,j], breaks=100,
        xlim=c(min(c(permGOStats[[i]][,j], hotspotGOStats[[i]][j]), na.rm=TRUE), max(c(permGOStats[[i]][,j], hotspotGOStats[[i]][j]), na.rm=TRUE)),
        main=paste(GOCats[i], j, sep=" "))
        abline(v = hotspotGOStats[[i]][j], lty=1, col="red")
    }
}
dev.off()
# the pi1 statistic doesn't work here - cutting off the 1s (which is necessary for these GO tests) leaves the distribution very jumpy
# for the # sig genes, BP is significant, esp in Bonferroni
# for the nominal pValues:
sapply(1:3, function(i){
    length(which(permGOStats[[i]][,"nominalTerms"] >= hotspotGOStats[[i]]["nominalTerms"])) / nrow(permGOStats[[i]])
})
#  0.001 0.060 0.695
# BP fine, MF just not, CC not at all

# OK, so the distribution of the GO enrichment statistic suggests signal in BP at least
# get permutation p-values per term:
permGOResultPerTerm <- lapply(1:3, function(x){
    sapply(rownames(permGO[[x]]), function(y){
        length(which(permGO[[x]][y,] <= hotspotGOResult[[x]][y,"result1"])) / ncol(permGO[[x]])
    })
})
# in each, how many are < 0.05?
sapply(permGOResultPerTerm, function(x){length(which(x < 0.05))})
# 352  99  18
sapply(permGOResultPerTerm, function(x){length(which(x < 0.005))})
# 149  53   1
# out of
sapply(permGOResultPerTerm, length)
# 3408 1122  810
# expectation:
# don't really have an expectation for the permutation-based p-values other than just the fraction of tests?
sapply(permGOResultPerTerm, length) * 0.005
#  17.04  5.61  4.05



# find false discovery rates at a range of p-values:
# i.e. how many terms are significant at a given threshold in the permutations vs observed
# do we need the mean or the max of the permutations?
pRange <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
sapply(1:3, function(x){
    t(sapply(pRange, function(thisP){
        mean(apply(permGO[[x]], 2, function(y){length(which(y < thisP))}))  / length(which(hotspotGOResult[[x]]$result1 < thisP))
    }))
})
# mean
#[,1]         [,2]     [,3]
#[1,] 0.271977612 0.3859180328 1.274429
#[2,] 0.079651852 0.0832857143      Inf
#[3,] 0.049393162 0.0526875000      Inf <= ~5% FDR at p=0.005
#[4,] 0.016666667 0.0125161290      Inf
#[5,] 0.009012346 0.0079655172      Inf
#[6,] 0.001878788 0.0007407407      Inf
#[7,] 0.001101695 0.0003703704      Inf
#[8,] 0.000325000 0.0000800000      Inf

#max
#[1,] 0.9253731 0.96721311 4.142857
#[2,] 0.7925926 0.88571429      Inf
#[3,] 0.8119658 0.84375000      Inf
#[4,] 0.6190476 0.35483871      Inf
#[5,] 0.5679012 0.27586207      Inf
#[6,] 0.4848485 0.11111111      Inf
#[7,] 0.3559322 0.07407407      Inf
#[8,] 0.1000000 0.04000000      Inf
# so, presumably, the mean

# how many terms are at p < 0.005?
# observed
sapply(1:3, function(x){length(which(hotspotGOResult[[x]]$result1 < 0.005))})
# 117   32  0
# expected
sapply(1:3, function(x){mean(apply(permGO[[x]], 2, function(y){length(which(y < 0.005))}))})
# 5.779 1.686 1.406


# and finally, a plot
# how about fold change vs significance, with FDR gene colored in?
# looks fine. Note that we always label the top five terms, as well as one term "cellular response to nutrient levels" because we present that one in the paper
pdf("hotspotGOEnrichments_170918.pdf", width=12, height=8)
par(mfcol=c(2,3))
for (i in 1:3){
    plot(log2(hotspotGOResult[[i]]$Significant / hotspotGOResult[[i]]$Expected), -log10(hotspotGOResult[[i]]$result1), xlab="log2(enrichment)", ylab="-log10(pValue)", pch=19, col="#00000022", cex=log(hotspotGOResult[[i]]$Annotated, base=20))
    points(
            log2(hotspotGOResult[[i]]$Significant / hotspotGOResult[[i]]$Expected)[hotspotGOResult[[i]]$result1 <= 0.005],
            -log10(hotspotGOResult[[i]]$result1)[hotspotGOResult[[i]]$result1 <= 0.005]
            , col="orange", cex=log(hotspotGOResult[[i]]$Annotated, base=20)[hotspotGOResult[[i]]$result1 <= 0.005])
    points(
            log2(hotspotGOResult[[i]]$Significant / hotspotGOResult[[i]]$Expected)[hotspotGOResult[[i]]$result1 <= 0.05/nrow(hotspotGOResult[[i]])],
            -log10(hotspotGOResult[[i]]$result1)[hotspotGOResult[[i]]$result1 <= 0.05/nrow(hotspotGOResult[[i]])]
            , col="red", cex=log(hotspotGOResult[[i]]$Annotated, base=20)[hotspotGOResult[[i]]$result1 <= 0.05/nrow(hotspotGOResult[[i]])])
    points(
            log2(hotspotGOResult[[i]]$Significant / hotspotGOResult[[i]]$Expected)[permGOResultPerTerm[[i]] <= 0.01 & hotspotGOResult[[i]]$result1 > 0.005],
            -log10(hotspotGOResult[[i]]$result1)[permGOResultPerTerm[[i]] <= 0.01 & hotspotGOResult[[i]]$result1 > 0.005]
            , col="blue", cex=log(hotspotGOResult[[i]]$Annotated, base=20)[permGOResultPerTerm[[i]] <= 0.01 & hotspotGOResult[[i]]$result1 > 0.005])
            #plot(log10(hotspotGOResult[[i]]$Annotated), -log10(hotspotGOResult[[i]]$result1), xlab="log10(category size)", ylab="-log10(pValue)", pch=19, col="#00000022")
            theseLabels <- c(1:5, which(hotspotGOResult[[i]]$Term == "cellular response to nutrient levels"))
            text(log2(hotspotGOResult[[i]]$Significant / hotspotGOResult[[i]]$Expected)[theseLabels], jitter(-log10(hotspotGOResult[[i]]$result1)[theseLabels], factor=1.6^(max(-log10(hotspotGOResult[[i]]$result1)[theseLabels]))), labels=hotspotGOResult[[i]]$Term[theseLabels], pos=2, cex=0.5)
     
     plot(hotspotGOResult[[i]]$Expected, hotspotGOResult[[i]]$Significant, xlab="Expected", ylab="Observed", pch=19, col="#00000022", cex=log(hotspotGOResult[[i]]$Annotated, base=20))
     
     points(hotspotGOResult[[i]]$Expected[hotspotGOResult[[i]]$result1 <= 0.005], hotspotGOResult[[i]]$Significant[hotspotGOResult[[i]]$result1 <= 0.005], xlab="Expected", ylab="Observed", col="orange", cex=log(hotspotGOResult[[i]]$Annotated, base=20)[hotspotGOResult[[i]]$result1 <= 0.005])
     points(hotspotGOResult[[i]]$Expected[hotspotGOResult[[i]]$result1 <= 0.05/nrow(hotspotGOResult[[i]])], hotspotGOResult[[i]]$Significant[hotspotGOResult[[i]]$result1 <= 0.05/nrow(hotspotGOResult[[i]])], xlab="Expected", ylab="Observed", col="red", cex=log(hotspotGOResult[[i]]$Annotated, base=20)[hotspotGOResult[[i]]$result1 <= 0.05/nrow(hotspotGOResult[[i]])])
     points(hotspotGOResult[[i]]$Expected[permGOResultPerTerm[[i]] <= 0.01 & hotspotGOResult[[i]]$result1 > 0.005], hotspotGOResult[[i]]$Significant[permGOResultPerTerm[[i]] <= 0.01 & hotspotGOResult[[i]]$result1 > 0.005], xlab="Expected", ylab="Observed", col="blue", cex=log(hotspotGOResult[[i]]$Annotated[permGOResultPerTerm[[i]] <= 0.01 & hotspotGOResult[[i]]$result1 > 0.005], base=20))
     
     abline(0, 1, lty=2, col="grey")
     text(hotspotGOResult[[i]]$Expected[theseLabels], jitter(hotspotGOResult[[i]]$Significant[theseLabels], factor=0.7^(mean(hotspotGOResult[[i]]$Significant[theseLabels]))), labels=hotspotGOResult[[i]]$Term[theseLabels], pos=4, cex=0.5)
}
dev.off()










################################################
# other external characteristics of hotspot genes

load("R_externalGeneData_170501.RData")

testGenes <- geneAnnotation[,"geneID"] %in% unlist(genesInHotspotList[!isEndMarker[,2]][which(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}) <= 3)])
names(testGenes) <- geneAnnotation[,"geneID"]
# 58 genes in here


pdf("hotspotGenesVsExternalGeneAnnotations.pdf", width=9, height=9)
par(mfrow=c(3,3))
for(i in c(1, 3:8)){
    boxplot(list(externalGeneData[names(testGenes[testGenes]), i], externalGeneData[names(testGenes[!testGenes]), i]), names=c("in hotspot", "not in hotspot"), ylab=names(externalGeneData)[i])
    legend("topright", box.lty=0, legend = paste0("p = ", round(t.test(externalGeneData[names(testGenes[testGenes]), i], externalGeneData[names(testGenes[!testGenes]), i])$p.value, 2)))
}
dev.off()

# nothing going on at all in these tests
# better to do via multiple logistic regression:

externalGeneDataWithHotspotGenes <- cbind(externalGeneData, as.numeric(testGenes[rownames(externalGeneData)]))
colnames(externalGeneDataWithHotspotGenes)[ncol(externalGeneDataWithHotspotGenes)] <- "inHotSpot"

summary(glm(inHotSpot ~ expression + essential + dNdS + PPI + GxGStrict + isTF + hasHumanHomolog + hasParalog, data=externalGeneDataWithHotspotGenes, family=binomial(link="logit")))
#Coefficients:
#Estimate Std. Error z value Pr(>|z|)
#(Intercept)         -5.0550059  1.1526541  -4.386 1.16e-05 ***
#expression           0.1235415  0.1487337   0.831   0.4062
#essentialTRUE       -1.8197674  1.0872497  -1.674   0.0942 .
#dNdS                -7.4769771  4.3900914  -1.703   0.0885 .
#PPI                  0.0004215  0.0026175   0.161   0.8721
#GxGStrict            0.0024931  0.0016693   1.494   0.1353
#isTFTRUE             2.5639101  0.5142486   4.986 6.17e-07 ***
#hasHumanHomologTRUE -1.2941383  0.5825589  -2.221   0.0263 *
#hasParalogTRUE       0.7576005  0.4727560   1.603   0.1090

# so there is your TF enrichment again :-)
# and marginally less likely to be essential or to have a human homolog

library(car)
Anova(glm(inHotSpot ~ expression + essential + dNdS + PPI + GxGStrict + isTF + hasHumanHomolog + hasParalog, data=externalGeneDataWithHotspotGenes, family=binomial(link="logit")))
#Response: inHotSpot
#LR Chisq Df Pr(>Chisq)
#expression        0.6807  1    0.40935
#essential         4.2255  1    0.03982 *
#dNdS              3.3521  1    0.06712 .
#PPI               0.0221  1    0.88182
#GxGStrict         1.9258  1    0.16522
#isTF             18.7080  1  1.523e-05 ***
#hasHumanHomolog   5.9904  1    0.01438 *
#hasParalog        2.4428  1    0.11807










################################################
# enrichment of local eQTL in hotspots?

# this is the table with corrected cis effects (one per gene, and with ASE data)
load("R_localeQTL_170222.RData")
rownames(localeQTL) <- localeQTL$gene

# are the tightly localized hotspots more likely to contain (strong) local eQTL than random regions?
localEffectsInHotspots <- lapply(names(genesInHotspotList), function(x){
    theseGenes <- genesInHotspotList[[x]][[1]]
    ret <- rep(NA, length(theseGenes))
    names(ret) <- theseGenes
    ret[names(ret) %in% localeQTL$gene] <- sapply(names(ret)[names(ret) %in% localeQTL$gene], function(i){localeQTL[localeQTL$gene == i,"lm.coeff"]})
    ret
})
names(localEffectsInHotspots) <- names(genesInHotspotList)

genesToSampleAll <- sapply(genesInHotspotList, function(y){length(y[[1]])})

geneAnno4Perms <- geneAnnotation[geneAnnotation$chr != "chrM",]
# exclude first/last 10 genes
last10GenesPerChr <- lapply(unique(geneAnno4Perms$chr), function(x){
    excludeGenes <- 10
    c(rep(TRUE, excludeGenes), rep(FALSE, length(which(geneAnno4Perms$chr == x)) - 2*excludeGenes), rep(TRUE, excludeGenes))
})
geneAnno4Perms <- geneAnno4Perms[!unlist(last10GenesPerChr),]

nPerms = 1000
permLocalEffectsHotspots <- lapply(1:nPerms, function(j){
    print(j)
    # construct the random set
    # there is a small chance that this will select some genes multiple times
    # ignore for now, would become problematic if there were more & bigger regions we'd sample
    sampledGenes <- lapply(genesToSampleAll, function(x){
        # pick one chr at random (note that this will weigh the chr by gene number because there is no unique(), which is fine)
        theseChrGenes <- geneAnno4Perms[geneAnno4Perms$chr == sample(geneAnno4Perms$chr, 1), "geneID"]
        # pick a gene, staying away far enough from the chr end to accomodate the number we want to include
        thisSampledAnchorGene <- sample(1:(length(theseChrGenes) - x), 1)
        # pick the anchor gene and the needed number of additional genes to its right
        theseChrGenes[thisSampledAnchorGene:(thisSampledAnchorGene + x - 1)]
    })
    # and get their local effects (if any)
    ret1 <- lapply(sampledGenes, function(y){
        ret <- rep(NA, length(y))
        names(ret) <- y
        ret[names(ret) %in% localeQTL$gene] <- sapply(names(ret)[names(ret) %in% localeQTL$gene], function(i){localeQTL[localeQTL$gene == i,"lm.coeff"]})
        ret
    })
    ret1
})
#save(permLocalEffectsHotspots, file="R_permLocalEffectsHotspots_170504.RData")
#load("R_permLocalEffectsHotspots_170504.RData")

# statistics to look at:
# how many hotspots have a local eQTL?
# are the local eQTL stronger than random?
# restrict to finely localized?

load("R_isEndMarker_170214.RData")
fineHotspots <- (!isEndMarker[,2]) & sapply(genesInHotspotList, function(x){length(x[[1]])}) <= 3

sum(sapply(localEffectsInHotspots[fineHotspots], function(x){length(which(!is.na(x)))}) > 0)
# max 3 genes: 21/26 have a local eQTL

# get the strongest local eQTL per hotspot
maxLocaleQTLPerHotspot <- sapply(localEffectsInHotspots[fineHotspots], function(x){ret<-NA; try({ret <- max(abs(x), na.rm=TRUE)}); ret})
summary(maxLocaleQTLPerHotspot, na.rm=TRUE)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.1106  0.1933  0.2613  0.3909  0.5766  0.8880       5


# now look at the random sets
maxLocaleQTLPerHotspotPerm <- sapply(permLocalEffectsHotspots, function(i){
    sapply(i[fineHotspots], function(x){ret<-NA; try({ret <- max(abs(x), na.rm=TRUE)}); ret})
})

# what about effect sizes?
#pdf("hotspotGenesVsLocalPerms_nonTelomere_3GeneHotspots.pdf", width=9, height=9)
pdf("hotspotGenesVsLocalPerms_nonTelomere_Hotspots.pdf", width=12, height=9)
par(mfrow=c(2,3))

hist(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){length(which(!is.na(x)))}), breaks=100, xlab="number of hotspots with local eQTL", main="number with local")
abline(v = sum(sapply(localEffectsInHotspots[fineHotspots], function(x){length(which(!is.na(x)))}) > 0), col="red", lty=1, lwd=2)
legend("topright", box.lty=0, legend=paste("p = ", length(which(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){length(which(!is.na(x)))}) >= sum(sapply(localEffectsInHotspots[fineHotspots], function(x){length(which(!is.na(x)))}) > 0))) / 1000, sep=""))

hist(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){length(which(!is.na(x))) / length(which(fineHotspots))}), breaks=100, xlab="fraction of hotspots with local eQTL", main="fraction with local", xlim=c(0, 1))
abline(v = sum(sapply(localEffectsInHotspots[fineHotspots], function(x){length(which(!is.na(x)))}) > 0) / length(which(fineHotspots)), col="red", lty=1, lwd=2)
legend("topright", box.lty=0, legend=paste("p = ", length(which(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){length(which(!is.na(x)))}) >= sum(sapply(localEffectsInHotspots[fineHotspots], function(x){length(which(!is.na(x)))}) > 0))) / 1000, sep=""))

hist(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){mean(x, na.rm=TRUE)}), breaks=100, xlab="mean strongest local eQTL effect per hotspot", main="mean of max local effect", xlim=c(min(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){mean(x, na.rm=TRUE)})), mean(maxLocaleQTLPerHotspot, na.rm=TRUE)))
abline(v = mean(maxLocaleQTLPerHotspot, na.rm=TRUE), col="red", lty=1, lwd=2)
legend("topright", box.lty=0, legend=paste("p = ", length(which(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){mean(x, na.rm=TRUE)}) >= mean(maxLocaleQTLPerHotspot, na.rm=TRUE))) / 1000, sep=""))

hist(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){median(x, na.rm=TRUE)}), breaks=100, xlab="median strongest local eQTL effect per hotspot", main="median of max local effect")
abline(v = median(maxLocaleQTLPerHotspot, na.rm=TRUE), col="red", lty=1, lwd=2)
legend("topright", box.lty=0, legend=paste("p = ", length(which(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){median(x, na.rm=TRUE)}) >= median(maxLocaleQTLPerHotspot, na.rm=TRUE))) / 1000, sep=""))

hist(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){max(x, na.rm=TRUE)}), breaks=100, xlab="maximum strongest local eQTL effect per hotspot", main="max of max local effect")
abline(v = max(maxLocaleQTLPerHotspot, na.rm=TRUE), col="red", lty=1, lwd=2)
legend("topright", box.lty=0, legend=paste("p = ", length(which(apply(maxLocaleQTLPerHotspotPerm, 2, function(x){max(x, na.rm=TRUE)}) >= max(maxLocaleQTLPerHotspot, na.rm=TRUE))) / 1000, sep=""))

dev.off()




# enrichment of ASE in hotspots genes?

load("R_localeQTL_170222.RData")
rownames(localeQTL) <- localeQTL$gene

testGenes <- geneAnnotation[,"geneID"] %in% unlist(genesInHotspotList[!isEndMarker[,2]][which(sapply(genesInHotspotList[!isEndMarker[,2]], function(x){length(x[[1]])}) <= 3)])
names(testGenes) <- geneAnnotation[,"geneID"]
# 58 genes in here

# do the hotspot genes have stronger ASE?
t.test(abs(localeQTL[names(testGenes)[which(testGenes)], "ASE"]), abs(localeQTL[names(testGenes)[which(!testGenes)], "ASE"]))
#t = 1.4236, df = 25.139, p-value = 0.1669
#mean of x mean of y
#0.3780978 0.2321898
# nope!
# although the mean is higher...
# wilcox: p=0.17


pdf("hotspotGenes_vs_ASE.pdf")
boxplot(list(abs(localeQTL[names(testGenes)[which(testGenes)], "ASE"]), abs(localeQTL[names(testGenes)[which(!testGenes)], "ASE"])), names=c("in hotspot", "not in hotspot"), ylab="abs(log2(ASE fold change))")
points(jitter(rep(1, length(localeQTL[names(testGenes)[which(testGenes)], "ASE"])), factor=10), abs(localeQTL[names(testGenes)[which(testGenes)], "ASE"]), col="blue", pch=19)
points(jitter(rep(2, length(localeQTL[names(testGenes)[which(!testGenes)], "ASE"])), factor=10), abs(localeQTL[names(testGenes)[which(!testGenes)], "ASE"]), col="#0000FF22")
dev.off()
# more, but not enough for significance

# are they more likely to *have* ASE?
wilcox.test(abs(localeQTL[names(testGenes)[which(testGenes)], "sigDeciderBonferroni"]), abs(localeQTL[names(testGenes)[which(!testGenes)], "sigDeciderBonferroni"]))
# p = 0.12

fisher.test(cbind(
    c(nrow(localeQTL[rownames(localeQTL) %in% names(testGenes[testGenes]) & localeQTL$sigDeciderBonferroni > 0 & !is.na(localeQTL$sigDeciderBonferroni),]),
    nrow(localeQTL[rownames(localeQTL) %in% names(testGenes[testGenes]) & localeQTL$sigDeciderBonferroni == 0 & !is.na(localeQTL$sigDeciderBonferroni),])),
    c(nrow(localeQTL[(!rownames(localeQTL) %in% names(testGenes[testGenes])) & localeQTL$sigDeciderBonferroni > 0 & !is.na(localeQTL$sigDeciderBonferroni),]),
    nrow(localeQTL[(!rownames(localeQTL) %in% names(testGenes[testGenes])) & localeQTL$sigDeciderBonferroni == 0 & !is.na(localeQTL$sigDeciderBonferroni),]))
))
# p = 0.16, OR=1.8

# nope, nothing significant here





###################
# fraction of distant linkages in hotspots

# simply do this based on overlaps for now
# what fraction of these overlaps a hotspot?

# distant
summary(overlapsAny(eQTLGRangesPad[!eQTLGRangesPad$cis], hotspotGRanges))
# 95%CI: 2948 FALSE 30581 TRUE => 9% do not overlap a hotspot

# local
summary(overlapsAny(eQTLGRangesPad[eQTLGRangesPad$cis], hotspotGRanges))
# 95% CI: 1790 FALSE, 1179 TRUE 60% are NOT "in" hotspots


# what fraction of the genome is covered by the hotspots?
sum(width(reduce(hotspotGRanges)))
# 1,080,947 ~1Mb, 8.3% of 12Mb
# only 8% if genome is in hotspots, but 91% of trans eQTL overlap a hotspot
# note that 'reduce' collapses them into a nonredundant set - result is the same with or without reduction





########################
# number of genes affected per hotspot
# using the B matrix

BGenesPerHotspot <- apply(B, 1, function(x){length(which(x != 0))})
summary(BGenesPerHotspot)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#26.0   168.0   425.0   660.0   709.5  4594.0

BGenesPerHotspot[which(BGenesPerHotspot == min(BGenesPerHotspot))]
# chrIII:166390_G/C
BGenesPerHotspot[which(BGenesPerHotspot == max(BGenesPerHotspot))]
# chrXIV:466588_T/G

# how many have > half the genes?
BGenesPerHotspot[which(BGenesPerHotspot >= ncol(B)/2)]
#chrXII:657022_T/C chrXIV:267161_G/A chrXIV:372376_G/A chrXIV:466588_T/G
#3640              3169              4172              4594






#####################
# HAP1 targets

HAP1Targets <- read.table("HAP1_targets_SGD_170623MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)

HAP1hotspotTargets <- B["chrXII:657022_T/C",]

HAP1topDownBY <- sort(HAP1hotspotTargets, decreasing=TRUE)[1:50]

length(which(HAP1Targets[,4] %in% names(HAP1hotspotTargets)))
# 69 of the targets are in the data

fisher.test(
    cbind(
        c(length(which((names(HAP1hotspotTargets) %in% names(HAP1topDownBY)) & (names(HAP1hotspotTargets) %in% HAP1Targets[,4]))),
            length(which((!names(HAP1hotspotTargets) %in% names(HAP1topDownBY)) & (names(HAP1hotspotTargets) %in% HAP1Targets[,4])))),
        c(length(which((names(HAP1hotspotTargets) %in% names(HAP1topDownBY)) & (!names(HAP1hotspotTargets) %in% HAP1Targets[,4]))),
            length(which((!names(HAP1hotspotTargets) %in% names(HAP1topDownBY)) & (!names(HAP1hotspotTargets) %in% HAP1Targets[,4]))))
    )
)

# 26   24
# 43 5536
# p < 2.2e-16, odds = 138

# better to just test for different effects for targets / non-targets?
t.test(HAP1hotspotTargets[names(HAP1hotspotTargets) %in% HAP1Targets[,4]], HAP1hotspotTargets[!names(HAP1hotspotTargets) %in% HAP1Targets[,4]])
# p = 2.189e-12, means 0.43001801 0.05866876



#######################
# CBF1 targeted analyses

# GO enrichment on ALL targets, not just the top 50/100

myGene2GO = read.table("gene_association_MOD_fromGO_160229.txt", stringsAsFactors = FALSE, sep="\t", header=FALSE)
# make that gene to GO list format they want:
gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
names(gene2GOList) = unique(myGene2GO[,1])

plotGOToTree <- function(GOdat, GOres, sigThres = 0.0005){
    # only plot if there are any significant GO terms (SEE ABOVE for "significance"; I am somewhat lenient here):
    # we need these extra lines because very small p-values are reported as a text string "< X", rather than a numeric
    toTest <- as.numeric(GenTable(GOdat, pVal = GOres)[1,6])
    if(is.na(toTest)){toTest <- 0.000000000000000000000000000001}
    if (toTest < sigThres){
        showSigOfNodes(GOdat, score(GOres), firstSigNodes = length(score(GOres)[score(GOres) < sigThres]), useInfo = "all")
    }else{
        plot(1,1)
        legend("topright", legend="no significant GO categories", cex=0.4, box.lty=0)
    }
}


testGenes <- factor(as.integer(B["chrX:548174_T/C", ] != 0))
names(testGenes) <- colnames(B)

# CAREFUL WITH THE FILE NAMES
# run > as UP, < as DOWN, != 0 as BOTH
pdf("topGO/CBF1/hotspotGeneGO_bootPeakGenes_BOTH.pdf", width=15, height=15)
for (thisOntology in c("BP", "MF", "CC")){
    GOData <- new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
    GOresult <- runTest(GOData, algorithm="classic", statistic="fisher")
    plotGOToTree(GOData, GOresult, sigThres = 0.005)
    
    write.table(GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult))), file=paste0("topGO/CBF1/hotspotGOresult_", thisOntology, "_bootPeakGenes_BOTH.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}
dev.off()


# now for TFBS
CBF1targets <- TFBSData[TFBSData[,1] == "CBF1", "targetSystematic"]

# alternatively, use the MacIssac 2006 predictions directly:
MacIsaac2006 <- readLines("orfs_by_factor_p0.001_cons2.txt")
MacIsaac2006 <- lapply(MacIsaac2006, function(x){strsplit(x, "\t")[[1]]})
names(MacIsaac2006) <- sapply(MacIsaac2006, function(x){x[1]})
MacIsaac2006 <- lapply(MacIsaac2006, function(x){x[2:length(x)]})
length(MacIsaac2006$CBF1)
# this is only 169 targets, compared to 210 from the total SGD data

CBF1targets <- MacIsaac2006$CBF1

# UP
fisher.test(
    cbind(
        c(length(which(B["chrX:548174_T/C", ] > 0 & colnames(B) %in% CBF1targets)),
            length(which(B["chrX:548174_T/C", ] > 0 & !(colnames(B) %in% CBF1targets)))),
        c(length(which(B["chrX:548174_T/C", ] <= 0 & colnames(B) %in% CBF1targets)),
            length(which(B["chrX:548174_T/C", ] <= 0 & !(colnames(B) %in% CBF1targets))))
    )
)
# SGD targets: p = 0.16, OR =1.8
# MacIssac targets: p = 0.06, OR=2.3

# DOWN
fisher.test(
cbind(
    c(length(which(B["chrX:548174_T/C", ] < 0 & colnames(B) %in% CBF1targets)),
        length(which(B["chrX:548174_T/C", ] < 0 & !(colnames(B) %in% CBF1targets)))),
    c(length(which(B["chrX:548174_T/C", ] >= 0 & colnames(B) %in% CBF1targets)),
        length(which(B["chrX:548174_T/C", ] >= 0 & !(colnames(B) %in% CBF1targets))))
    )
)
# SGD targets: p 8e-12, OR 5.3
# MacIssac targets: p = 2e-10, OR=5.4

# BOTH
fisher.test(
cbind(
    c(length(which(B["chrX:548174_T/C", ] != 0 & colnames(B) %in% CBF1targets)),
        length(which(B["chrX:548174_T/C", ] != 0 & !(colnames(B) %in% CBF1targets)))),
    c(length(which(B["chrX:548174_T/C", ] == 0 & colnames(B) %in% CBF1targets)),
        length(which(B["chrX:548174_T/C", ] == 0 & !(colnames(B) %in% CBF1targets))))
    )
)
# SGD targets: p = 3e-11, OR = 4.3
# MacIssac targets: p = 1e-10, OR=4.6


# which are the strongest effects at CBF1 that are NOT CBF1 targets?
length(B["chrX:548174_T/C", ][B["chrX:548174_T/C", ] != 0 & (colnames(B) %in% CBF1targets)])
# 37 ARE targets (32 with MacIssac)
length(B["chrX:548174_T/C", ][B["chrX:548174_T/C", ] != 0 & !(colnames(B) %in% CBF1targets)])
# 270 are NOT targets (275 with MacIsaac)

sort(B["chrX:548174_T/C", ][B["chrX:548174_T/C", ] != 0 & (colnames(B) %in% CBF1targets)])
# strongest effect (negative) that is also CBF1 target is CAP1, other targets like this are SNC2 (Vesicle membrane receptor protein (v-SNARE)) and MSH3 (Mismatch repair protein)
# strongest effect (positive) that is also CBF1 target is RPL11A, EPL1 (Subunit of NuA4, an essential histone H4/H2A acetyltransferase complex), NDI1 (NADH:ubiquinone oxidoreductase); these are all smaller effect than the negative genes above

sort(B["chrX:548174_T/C", ][B["chrX:548174_T/C", ] != 0 & !(colnames(B) %in% CBF1targets)])
# strongest effect (negative) is SLT2 (Serine/threonine MAP kinase), followed by CRH1 (a chitin transglycosylase), followed by DFG5 (Putative mannosidase); all three are cell wall; all three have larger absolute effects than the genes in the positive direction (with the exception of that local gene, which may be a random, colocalized cis eQTL)
# strongest effect (positive) is YJR061W (a local gene), followed by YIL092, uncharacterized, followed by HFM1 (Meiosis specific DNA helicase)





#####################
# 10/09/2017
# pull in all GO results and print as ONE table

# best 5 GO terms
# NOTE RESULTSFOLDER
allHotspotGO <- lapply(rownames(hotspot.boot.intervals), function(x){
    thisFile <- gsub(":", "_", x)
    thisFile <- gsub("/", "_", thisFile)
    lapply(c("topGenes100"), function(thisGeneNumber){
        lapply(c("BP", "MF", "CC"), function(thisGOCat){
            lapply(c("Up", "Down"), function(thisUpDown){
                thisFile1 <- paste0("RESULTSFOLDER/topGO/", thisGeneNumber, "/", thisGOCat, "/GOGraph_hotspot_", thisFile, "_Fisher", thisUpDown, "_GOResultTable.txt", sep="")
                thisGOResult <- read.table(thisFile1, stringsAsFactors=FALSE, head=TRUE, sep="\t", quote="\"")
                thisGOResult <- cbind(x, as.character(thisGOCat), as.character(thisUpDown), thisGOResult, stringsAsFactors=FALSE)
                thisGOResult[thisGOResult == "Up"] <- "RM"
                thisGOResult[thisGOResult == "Down"] <- "BY"
                thisGOResult
            })
        })
    })
})
names(allHotspotGO) <- rownames(hotspot.boot.intervals)

# print into one giant table
# this file gets to be ~100Mb, even when not printing the "both" direction
# remove the text descriptions - people will need to parse this themselves
write.table(t(c("hotspot", "GO_subtree", "directionOfHotspotEffect", "GO_ID", "numberOfGenesInTerm", "significantGenesInTerm", "expectedGenesInTerm", "pValue")), file="allGOResults4paper_171009.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
for (i in (1:length(allHotspotGO))){
    for (j in 1:length(allHotspotGO[[i]][[1]])){
        for (k in 1:length(allHotspotGO[[i]][[1]][[j]])){
            write.table(allHotspotGO[[i]][[1]][[j]][[k]][,-c(5)], file="allGOResults4paper_171009.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
        }
    }
}




#####################
# 10/09/2017
# pull in all TFBS results and print as ONE table

load("R_TFBSEnrichments_top100Genes_170420.RData")

printDirections <- c("RM", "BY")

write.table(t(c("hotspot", "directionOfHotspotEffect", "regulatorGene", "pValue", "odds", "regulatedAndTarget", "regulatedAndNotTarget", "notRegulatedAndTarget", "notRegulatedAndNotTarget")), file="TFBS_allResults4paper_171009.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
for (i in (1:length(TFBSEnrichments))){
    for (j in 1:2){
            write.table(cbind(names(TFBSEnrichments)[i], printDirections[j], TFBSEnrichments[[i]][[j]]), file="TFBS_allResults4paper_171009.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=FALSE, append=TRUE)
    }
}





