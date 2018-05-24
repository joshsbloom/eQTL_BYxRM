library(GenomicRanges)
library(qvalue)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#####################
# load many objects (not all may actually be needed for these analyses for the revision)


serverPrefix = ""

load(paste0(serverPrefix, "gbatch_fact.RData", sep=""))
load(paste0(serverPrefix, "ODcov.RData", sep=""))
load(paste0(serverPrefix, "log2_t.tpm.matrix.RData", sep=""))

# a cleaned-up version for plots
phenoBatchOD = apply(t.tpm.matrix, 2, function(y){scale(lm(y ~ gbatch.fact + OD.cov)$res)})
rownames(phenoBatchOD) <- sapply(rownames(t.tpm.matrix), function(x){strsplit(x, "-")[[1]][1]})

load(paste0(serverPrefix, "gdata_42k.RData", sep=""))

hotspot.boot.intervals <- read.table(paste0(serverPrefix, "MV.12417_mod.txt", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE)
rownames(hotspot.boot.intervals) <- hotspot.boot.intervals[,1]
hotspot.boot.intervals <- hotspot.boot.intervals[,2:6]

load(paste0(serverPrefix, "mv.bootstrap.pos.RData", sep=""))
hotspot.list <- list(rownames(hotspot.boot.intervals), mv.boots, hotspot.boot.intervals)
names(hotspot.list[[2]]) <- hotspot.list[[1]]
rownames(hotspot.list[[3]]) <- hotspot.list[[1]]

load(paste0(serverPrefix, "B.Forward.RData", sep=""))
load(paste0(serverPrefix, "B.cis.Forward.RData", sep=""))

geneAnnotation = read.table(paste0(serverPrefix, "ensemblGenes_ensembl83_160307_MOD.txt", sep=""), sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]

allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

load("marker.LD.RData")

getLDBlock <- function(markerDat, thisMarker, LDCut = 1){
    thisChr <- strsplit(thisMarker, ":")[[1]][1]
    allLDmarkers <- colnames(markerDat[[thisChr]])[markerDat[[thisChr]][thisMarker,] >= LDCut]
    return(c(allLDmarkers[1], allLDmarkers[length(allLDmarkers)]))
}

getClosestMarker <- function(markerDat, thisChr, thisPosition){
    thesePositions <- sapply(colnames(markerDat[[thisChr]]), function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})
    colnames(markerDat[[thisChr]])[which.min(abs(thesePositions - thisPosition))][1]
}

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

# replace original object:
geneGRanges <- geneGRangesExtended

# eQTL:
load(paste0(serverPrefix, "R_allPeaksODPad_161213.RData", sep=""))
eQTLGRangesPad <- GRanges(seqnames = allPeaksODPad[,"chr"], ranges = IRanges(start=sapply(allPeaksODPad[,"CI.l"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])}), end=sapply(allPeaksODPad[,"CI.r"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})))
eQTLGRangesPad$gene <- allPeaksODPad$gene
eQTLGRangesPad$cis <- allPeaksODPad$cis

eQTLGRangesPad_withGeneNames <- GRanges(seqnames = paste(allPeaksODPad[,"chr"], allPeaksODPad[,"gene"],  sep="_"), ranges = IRanges(start=sapply(allPeaksODPad[,"CI.l"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])}), end=sapply(allPeaksODPad[,"CI.r"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})))
eQTLGRangesPad_withGeneNames$gene <- allPeaksODPad$gene
eQTLGRangesPad_withGeneNames$cis <- allPeaksODPad$cis


# load the actual epistasis results:
load("R_lm2D_between_additive.RData")
load("R_lm2D_from_full.RData")

# replace factors in object
lm2D_from_full %>% mutate_if(is.factor, as.character) -> lm2D_from_full
lm2D_between_additive %>% mutate_if(is.factor, as.character) -> lm2D_between_additive




############################
# analyses start here

# for the full 2D:
# how many pairs have 0, 1, 2 overlap a detected additive QTL?

thisPad <- 10e6
overlapsAdditive <- t(apply(lm2D_from_full, 1, function(x){
    # make GRanges for the two markers
    pairMarkers <- lapply(c("m1_pmarker", "m2_pmarker"), function(i){
        markerSplit <- strsplit(x[i], ":")[[1]]
        markerChr <- markerSplit[1]
        markerPos <- as.integer(strsplit(markerSplit[2], "_")[[1]][1])
        markerGRange <- GRanges(seqnames = paste(markerChr, x["trait"], sep="_"), ranges = IRanges(start = markerPos - thisPad, end= markerPos + thisPad))
        markerGRange
    })
    sapply(pairMarkers, function(x){
        overlapsAny(x, eQTLGRangesPad_withGeneNames)
    })
}))
dim(overlapsAdditive)
# 387 pairs, i.e. 774 markers

table(rowSums(overlapsAdditive))
# 10k padding
#0   1   2
#54 108 225
# => (54 * 2 + 108)/774 = 28% have no additive QTL; 72% do

# 25k padding (more conservative)
#0   1   2
#40  81 266

# 10M padding (guarantees no other QTL on this chr; this is most conservative for "interesting" examples)
#0   1   2
#11  40 336

# pi1 for the p-values
full2DAdditivePs <- c(lm2D_from_full$m1_p, lm2D_from_full$m2_p)
1 - qvalue(full2DAdditivePs)$pi0
# 0.8421334
# at least 84% have additive signal

# focus on those where none overlap; anything good?
noAdditive <- lm2D_from_full[rowSums(overlapsAdditive) == 0 & lm2D_from_full$m1_p > 0.05 & lm2D_from_full$m2_p > 0.05,]
# 9 such with the strongest filter of 10M pad
noAdditive[order(noAdditive$interaction_variance_explained, decreasing=TRUE),]
# all are trans-trans
save(noAdditive, file="R_noAdditive_180410.RData")

allPeaksODPad[allPeaksODPad$gene == "YML116W",]


# which pairs overlap a hotspot?
# use a less aggressive pad – every chr has a hotspot, so the huge pad would be true for all
thisPadForHotspot <- 20000
overlapsHotspot <- t(apply(lm2D_from_full, 1, function(x){
    # make GRanges for the two markers
    pairMarkers <- lapply(c("m1_pmarker", "m2_pmarker"), function(i){
        markerSplit <- strsplit(x[i], ":")[[1]]
        markerChr <- markerSplit[1]
        markerPos <- as.integer(strsplit(markerSplit[2], "_")[[1]][1])
        markerGRange <- GRanges(seqnames = markerChr, ranges = IRanges(start = markerPos - thisPadForHotspot, end= markerPos + thisPadForHotspot))
        markerGRange
    })
    sapply(pairMarkers, function(x){
        overlapsAny(x, hotspotGRanges)
    })
}))
table(rowSums(overlapsHotspot))
#0   1   2
#39 141 207

# are the non-additive pairs truly alone, or do they perhaps correspond to hotspots after all, if only epistatically?
noAdditiveNoHotspot <- lm2D_from_full[rowSums(overlapsAdditive) == 0 & lm2D_from_full$m1_p > 0.05 & lm2D_from_full$m2_p > 0.05 & rowSums(overlapsHotspot) == 0,]
noAdditiveNoHotspot[order(noAdditiveNoHotspot$interaction_variance_explained, decreasing=TRUE),]
# three in total, two share the same chrIX marker chrIX:83087_T/C
# all three are on chrVII: chrVII:275701_C/T, chrVII:329109_A/G, chrVII:344309_T/G

#trait        m1_pmarker       m1_beta m1_variance_explained      m1_p
#177 YOR300W chrVII:275701_C/T -0.0266125540          4.461644e-04 0.3893026 <- dubious ORF on chr15, overlaps BUD7; off in my FP data
#152 YAL028W chrVII:329109_A/G  0.0363199107          1.141276e-03 0.2418441 <- FRT2, tail-anchored ER membrane protein of unknown function
#162 YNR033W chrVII:344309_T/G  0.0008529821          1.071796e-06 0.9781029 <- ABZ1, Para-aminobenzoate (PABA) synthase
#m1_is_cis        m2_pmarker       m2_beta m2_variance_explained      m2_p
#177     FALSE chrXVI:592674_T/C -0.0232510050          8.253172e-04 0.4519961
#152     FALSE   chrIX:83087_T/C  0.0303644157          7.331300e-04 0.3278494
#162     FALSE   chrIX:83087_T/C -0.0008437293          9.823003e-06 0.9783410
#m2_is_cis interaction_beta interaction_variance_explained interaction_p
#177     FALSE       -0.1890525                     0.03572262  1.354594e-09
#152     FALSE        0.1647656                     0.02715988  1.329715e-07
#162     FALSE        0.1613093                     0.02603751  2.493301e-07

# plot the expression levels for these three:
plotOnePair <- function(gene, marker1, marker2){
    plotList <- list(
        phenoBatchOD[rownames(gdata)[gdata[,marker1] == -1 & gdata[,marker2] == -1], gene],
        phenoBatchOD[rownames(gdata)[gdata[,marker1] == -1 & gdata[,marker2] == 1], gene],
        phenoBatchOD[rownames(gdata)[gdata[,marker1] == 1 & gdata[,marker2] == -1], gene],
        phenoBatchOD[rownames(gdata)[gdata[,marker1] == 1 & gdata[,marker2] == 1], gene]
    )
    ggDat <- data.frame(expression = unlist(plotList), genotype=rep(c("BB", "BR", "RB", "RR"), times = sapply(plotList, length)))
    p <- ggplot(ggDat, aes(x = genotype, y = expression)) + geom_boxplot(color="darkblue", outlier.shape=NA) + theme_light() + labs(x="Genotype", y="Normalized expression level", title=paste(gene, allNames[gene], sep=" - "), subtitle=paste(marker1, marker2, sep=" and ")) +
    geom_jitter(width=0.2, alpha=0.4)
    print(p)
}


pdf("epistasis_noAdditiveNoHotspotExamples.pdf")
for (i in order(noAdditiveNoHotspot$interaction_variance_explained, decreasing=TRUE)){
    plotOnePair(noAdditiveNoHotspot[i, "trait"], noAdditiveNoHotspot[i, "m1_pmarker"], noAdditiveNoHotspot[i, "m2_pmarker"])
}
dev.off()
# in each of the three noAdditive noHotspot, the parents are similar, both hets are higher or lower

# plot all the nine non-additive
pdf("epistasis_noAdditiveExamples.pdf")
for (i in order(noAdditive$interaction_variance_explained, decreasing=TRUE)){
    plotOnePair(noAdditive[i, "trait"], noAdditive[i, "m1_pmarker"], noAdditive[i, "m2_pmarker"])
}
dev.off()
# none of these are visually more striking than the three above; in fact YOR300W is the strongest of these overall

# plot the ten biggest 2D effects, irrespective of how they relate to additive effects
pdf("epistasis_strongest2DExamples.pdf")
for (i in order(lm2D_from_full$interaction_variance_explained, decreasing=TRUE)[1:10]){
    plotOnePair(lm2D_from_full[i, "trait"], lm2D_from_full[i, "m1_pmarker"], lm2D_from_full[i, "m2_pmarker"])
}
dev.off()

lm2D_from_full[order(lm2D_from_full$interaction_variance_explained, decreasing=TRUE)[1:10],c("trait", "m1_pmarker", "m2_pmarker")]
#trait        m1_pmarker         m2_pmarker
#165 YGR027W-B chrVII:529347_A/T  chrXII:818534_G/T <- Ty-GAG retrotransposon; itself with a region where my RNASeq suggests a BY-specific, unannotated ncRNA upstream of a tRNA. An insertion of this gene in BY? It's mostly gone in RR
#46    YDR461W chrIII:204750_G/A  chrIV:1385256_C/A <- MFA1, between itself and the mating type locus. Nice! it's only on in BB. why is not on in BR, when the cell is 'a', but MFA1 is RM?
#65    YJR004C chrIII:204750_G/A chrVIII:103046_A/G <- SAG1, Alpha-agglutinin of alpha-cells: mating type & GPA1; on only in RB
#143   YFL054C   chrVI:30187_G/A     chrX:28380_C/T <- AQY3, might be chr ends? OFF in RB
#132   YOR133W chrIV:1248530_C/G   chrXV:582005_C/T <- EFT1, has close paralog EFT2 ("also encoded by"), interaction is paralog with itself, off in BR. alignments were almost impossible in Albert14
#144   YFL055W   chrVI:30187_G/A     chrX:28380_C/T <- AGP3, next to AQY3, see above, same markers, same pattern
#206 YHL009W-B chrVIII:83441_C/T  chrXVI:434157_C/T <- retrotransposon, itself with another hypothetical Ty-GAG, OFF in RR. Same as first one above
#209   YJL225C   chrIX:28874_G/C     chrX:28380_C/T <- telomeres
#78    YCL019W  chrIII:81425_C/A  chrXIV:560313_A/T <- Ty-GAG, next to a region without any aligned RNA reads and no annotations. another insertion? off in RR
#86  YDR210W-B  chrIV:878249_A/C  chrVII:405432_C/T <- Ty-Gag, itself with region next to tRNA




##########################################################
# for 2D:
# where are they in the genome?
# bin them and additive loci, make histogram-style plot
chromosomeDividers <- c(0, 230218, 1043402, 1360022, 2891955, 3468829, 3738990, 4829930, 5392573, 5832461, 6578212, 7245028, 8323205, 9247636, 10031969, 11123260, 12071326)
getGcoords = function ( chr , pos, spacing=0, sgd.table="sacCer3ChromLenghts.txt" ) {
    offind = as.vector(cumsum(read.table(sgd.table, header=FALSE, sep="\t")[,2] + spacing))
    offind=    offind[-length(offind)]
    offind= c(0, offind)
    names(offind) = as.character(read.table(sgd.table, header=FALSE, sep="\t")[,1])
    chr.off=as.numeric(sapply(chr, function(x) {offind[[x]]}))
    return(chr.off+pos)
}
sepBetweenChr <- 1e5

# need to make plotable positions:
# eQTL:
eQTLGCoords <- sapply(allPeaksODPad$pmarker, function(x){
    thisSplit <- strsplit(x, ":")[[1]]
    print(thisSplit)
    getGcoords(thisSplit[1], as.numeric(strsplit(thisSplit[2], "_")[[1]][1]))
})

pairs2DGcoords <- sapply(c(lm2D_from_full$m1_pmarker, lm2D_from_full$m2_pmarker), function(x){
    thisSplit <- strsplit(x, ":")[[1]]
    getGcoords(thisSplit[1], as.numeric(strsplit(thisSplit[2], "_")[[1]][1]))
})


doubleHistPlot <- function(xVals, yVals1, yVals2, thisYLim = NA, thisMain="", thisYLabel="", yFract=TRUE){
    if(is.na(thisYLim)){thisYLim <- c(-max(yVals2), max(yVals1))}
    
    plot(xVals, yVals1, type="h", main=thisMain, xaxt="n", xlab="genome position (chromosome name)", col="black", ylab=thisYLabel, yaxt="n",
    ylim=thisYLim)
    points(xVals, -yVals2, type="h")
    for(i in 2:(length(chromosomeDividers)-1)){
        abline(v=chromosomeDividers[i], col="grey", lty=2)
    }
    axis(1, at=sapply(1:16, function(i){chromosomeDividers[i] + ((chromosomeDividers[i+1] - chromosomeDividers[i])/2)}), labels=as.roman(1:16), tick=FALSE)
    
    if(yFract){axis(2, at = seq(from=-1, to=1, by=0.025), labels=abs(seq(from=-1, to=1, by=0.025)))}
    else{axis(2, at = pretty((-max(yVals2)):max(yVals1) ), labels=pretty((-max(yVals2)):max(yVals1) ))}
    legend("topleft", legend="eQTL markers", box.lty=0)
    legend("bottomleft", legend="pair markers", box.lty=0)
    
}

pdf("epistasis_2D_locationHistogram.pdf", width=9, height=5)

eQTLHist <- hist(eQTLGCoords, breaks=600, plot=FALSE)
pairQTLHist <- hist(pairs2DGcoords, breaks=eQTLHist$breaks, plot=FALSE)

doubleHistPlot(eQTLHist$mids, eQTLHist$counts/length(unique(allPeaksODPad$gene)), pairQTLHist$counts/length(unique(lm2D_from_full$trait)), thisMain = "all eQTL vs 2D pair markers", thisYLabel="fraction of genes with an eQTL or epistatic pair")

doubleHistPlot(eQTLHist$mids, eQTLHist$counts/length(eQTLGCoords), pairQTLHist$counts/(length(pairs2DGcoords)/2), thisMain = "all eQTL vs 2D pair markers", thisYLabel="fraction of eQTL or epistatic pairs")

doubleHistPlot(eQTLHist$mids, eQTLHist$counts, pairQTLHist$counts, thisMain = "all eQTL vs 2D pair markers", thisYLabel="number of eQTL or epistatic pairs", yFract=FALSE)

dev.off()
# clear agreement: HAP1, MKT1, IRA2, GPA1
# but also differences: mating type seems overenriched compared to additive?
# several additive hotspots that are NOT epistaistic hotspots: 12 left, 11, 7
# little evidence for epistatic-only hotspots

# what is that thing on chrIV that has no hotspot?
head(sort(table(c(lm2D_from_full$m1_pmarker, lm2D_from_full$m2_pmarker)), decreasing=TRUE), n=10)
#chrVIII:103046_A/G  chrXII:657022_T/C  chrXIV:376313_C/T  chrXII:648649_C/T
#29                 22                 22                 18
#chrXIV:373628_AT/A  chrIII:204750_G/A  chrXIV:449640_A/G   chrIV:997624_A/G
#18                 13                 12                 11
#chrXII:640445_G/A chrXIII:913429_G/A
#9                  9

# the one on IV just misses chrIV:975156_A/G (YAP6?) with 20k padding



################################
# which hotspots do they overlap?

hotspotsFor2D <- t(apply(lm2D_from_full, 1, function(x){
    # make GRanges for the two markers
    pairMarkers <- lapply(c("m1_pmarker", "m2_pmarker"), function(i){
        markerSplit <- strsplit(x[i], ":")[[1]]
        markerChr <- markerSplit[1]
        markerPos <- as.integer(strsplit(markerSplit[2], "_")[[1]][1])
        markerGRange <- GRanges(seqnames = markerChr, ranges = IRanges(start = markerPos - thisPadForHotspot, end= markerPos + thisPadForHotspot))
        markerGRange
    })
    sapply(pairMarkers, function(x){
        findOverlaps(x, hotspotGRanges)
    })
}))
# how many each?
summary(t(sapply(hotspotsFor2D, function(x){c(length(x[[1]]), length(x[[2]]))})))
# max is two

# which additive hotspots are most frequently overlapped by the pairs?
hotspotsThatOverlap2D <- table(unlist(sapply(unlist(hotspotsFor2D), subjectHits)))
names(hotspotsThatOverlap2D) <- rownames(hotspot.boot.intervals)[as.numeric(names(hotspotsThatOverlap2D))]
t(t(head(sort(hotspotsThatOverlap2D, decreasing=TRUE), n=10)))
#chrXII:657022_T/C    79    HAP1
#chrXIV:372376_G/A    66    KRE33
#chrXIV:466588_T/G    49    MKT1
#chrVIII:104423_C/A   45    GAP1
#chrXV:150803_G/A     28    HAL9
#chrIII:201185_C/A    22    mating type
#chrXV:171150_T/C     20    IRA2
#chrXIII:33836_G/A    16    PHO84?
#chrIII:143912_A/G    14    ???
#chrXVI:26510_T/C     12    telomeric

# which pairs of hotspots are the most frequent?
overlapsWithHotspotsPairwise <- matrix(0, nrow=nrow(hotspot.boot.intervals), ncol=nrow(hotspot.boot.intervals))
rownames(overlapsWithHotspotsPairwise) <- rownames(hotspot.boot.intervals)
colnames(overlapsWithHotspotsPairwise) <- rownames(hotspot.boot.intervals)
for(i in 1:length(hotspotsFor2D)){
    marker1HS <- subjectHits(hotspotsFor2D[[i]][[1]])
    marker2HS <- subjectHits(hotspotsFor2D[[i]][[2]])
    # not every marker overlaps a hotspot
    if(length(marker1HS) > 0 & length(marker2HS) > 0){
    # there can be more than one overlap between an interaction and a hotspot...
    for (i1 in 1:length(marker1HS)){
            for (i2 in 1:length(marker2HS)){
    # the interactions have no direction, so only add in one half of the matrix
                overlapsWithHotspotsPairwise[min(c(marker1HS[i1], marker2HS[i2])), max(c(marker1HS[i1], marker2HS[i2]))] <- overlapsWithHotspotsPairwise[min(c(marker1HS[i1], marker2HS[i2])), max(c(marker1HS[i1], marker2HS[i2]))] + 1
        }
    }
    }
}
overlapsWithHotspotsPairwise[rowSums(overlapsWithHotspotsPairwise) > 6, colSums(overlapsWithHotspotsPairwise) > 6]

write.table(overlapsWithHotspotsPairwise, file="overlapsWithHotspotsPairwise_180422.txt", sep="\t", quote=FALSE)


# for between-additive:
hotspotsForBA <- t(apply(lm2D_between_additive, 1, function(x){
    # make GRanges for the two markers
    pairMarkers <- lapply(c("m1_pmarker", "m2_pmarker"), function(i){
        markerSplit <- strsplit(x[i], ":")[[1]]
        markerChr <- markerSplit[1]
        markerPos <- as.integer(strsplit(markerSplit[2], "_")[[1]][1])
        markerGRange <- GRanges(seqnames = markerChr, ranges = IRanges(start = markerPos - thisPadForHotspot, end= markerPos + thisPadForHotspot))
        markerGRange
    })
    sapply(pairMarkers, function(x){
        findOverlaps(x, hotspotGRanges)
    })
}))
# how many each?
summary(t(sapply(hotspotsForBA, function(x){c(length(x[[1]]), length(x[[2]]))})))
# max is two

# which additive hotspots are most frequently overlapped by the pairs?
hotspotsThatOverlapBA <- table(unlist(sapply(unlist(hotspotsForBA), subjectHits)))
names(hotspotsThatOverlapBA) <- rownames(hotspot.boot.intervals)[as.numeric(names(hotspotsThatOverlapBA))]
hotspotsThatOverlapBA

# correlate # pairs with # genes affected
genesAffectedHS <- apply(B, 1, function(x){length(which(x != 0))})
pairsPerHS <- rep(0, length(genesAffectedHS))
names(pairsPerHS) <- names(genesAffectedHS)
pairsPerHS[names(hotspotsThatOverlapBA)] <- hotspotsThatOverlapBA

cor.test(pairsPerHS, genesAffectedHS)
# r=0.78, p < 2.2e-16
pdf("epistasis_pairsVsGenesAffectedPerHS.pdf")
plot(pairsPerHS, genesAffectedHS)
dev.off()


###############################################
# cis x trans

# start with HAP1 targets; if see polarization, run search for polarization at other loci to determine more "active" trans allele
# find all hotspots where one partner is HAP1, not cis, and the other is cis

# CAREFUL with variable names, this is not great coding
hotspotsFor2Dadditive <- t(apply(lm2D_between_additive, 1, function(x){
    # make GRanges for the two markers
    pairMarkers <- lapply(c("m1_pmarker", "m2_pmarker"), function(i){
        markerSplit <- strsplit(x[i], ":")[[1]]
        markerChr <- markerSplit[1]
        markerPos <- as.integer(strsplit(markerSplit[2], "_")[[1]][1])
        markerGRange <- GRanges(seqnames = markerChr, ranges = IRanges(start = markerPos - thisPadForHotspot, end= markerPos + thisPadForHotspot))
        markerGRange
    })
    sapply(pairMarkers, function(x){
        findOverlaps(x, hotspotGRanges)
    })
}))
# full2D <- lm2D_from_full
# hsF2D <- hotspotsFor2D

#lm2D_from_full <- lm2D_between_additive
#hotspotsFor2D <- hotspotsFor2Dadditive

thisHS <- "chrXII:657022_T/C" #HAP1
#thisHS <- "chrXIV:372376_G/A" # KRE33
#thisHS <- "chrXIV:466588_T/G" # MKT1
#thisHS <- "chrVIII:104423_C/A" # GAP1
#thisHS <- "chrXV:171150_T/C" # IRA2
#thisHS <- "chrIII:201185_C/A" # mating type => only two cases, one of which is MFA1
thisHSIndex <- which(rownames(hotspot.boot.intervals) == thisHS)
inThisHS <- t(sapply(hotspotsFor2D, function(x){
    c(thisHSIndex %in% subjectHits(x[[1]]), thisHSIndex %in% subjectHits(x[[2]]))
}))

thisHS2D <- cbind(lm2D_from_full, inThisHS)[rowSums(inThisHS) > 0,]
thisOneIsTheHS <- sapply(apply(thisHS2D[,15:16], 1, which), function(x){c("m1_is_cis", "m2_is_cis")[x]})
thisOneIsNotTheHS <- sapply(apply(!thisHS2D[,15:16], 1, which), function(x){c("m1_is_cis", "m2_is_cis")[x]})
# HS not cis AND other one is actually cis
thisHS2D <- thisHS2D[!sapply(1:length(thisOneIsTheHS), function(i){thisHS2D[i, thisOneIsTheHS[i]]}) & sapply(1:length(thisOneIsNotTheHS), function(i){thisHS2D[i, thisOneIsNotTheHS[i]]}),]
# for HAP1, this has 15 rows

# make the plots for these
pdf("epistasis_HAP1_by_local.pdf")
for (i in order(thisHS2D$interaction_variance_explained, decreasing=TRUE)){
    plotOnePair(thisHS2D[i, "trait"], thisHS2D[i, "m1_pmarker"], thisHS2D[i, "m2_pmarker"])
}
dev.off()
# YHR054C has one cis allele that shuts it off, the on one reacts to HAP1

# looking through the plots, it does NOT seem like these genes are known HAP1 targets??
# do more systematically:
TFBSData <- read.table(paste0(serverPrefix, "/falbert/eQTL1k/annotations/regulationRelationshipsFromSGD_160501_MOD.txt", sep=""), head=TRUE, sep="\t", quote="\"", stringsAsFactors=FALSE)
colnames(TFBSData) <- c("TFName", "TFSystematic", "targetName", "targetSystematic", "evidenceType", "evidenceCode", "condition", "construct", "strain", "evidenceObservation", "refID", "SGD")
thisHS2D$trait %in% TFBSData$targetSystematic
# all true; they are targets of SOMETHING
thisHS2D$trait %in% TFBSData$targetSystematic[TFBSData$TFName == "HAP1"]
# 75 direct targets
# only one in this set
thisHS2D[thisHS2D$trait %in% TFBSData$targetSystematic[TFBSData$TFName == "HAP1"],]
# YGR049W – SCM4
# shows expected pattern of higher with RM HAP1
# RM has higher cis; this effect is stronger with RM HAP1; all three betas are positive
# does SCM4 have two HAP1 sites (so it will be turned on either way), only one of which has a variant?
# or just one that is weaker but still responds?
# Parts & Stegle mapped this to inferred HAP1 activity
# no for the proximal C/T SNP, but the RM SNP creates a (two?) Hsf1 BS
# distal C/G: BY, but not RM, has a Mot3p BS

# how many of the cis partners are significant by themselves?


# now we need the segregants with the BY and RM allele at the HS:
HSBYSegs <- rownames(gdata)[gdata[,thisHS] == -1]
HSRMSegs <- rownames(gdata)[gdata[,thisHS] == 1]

cisEffectsSplitByHS <- data.frame(t(sapply(1:nrow(thisHS2D), function(i){
    thisOneIsTheCisMarker <- thisHS2D[i, c("m1_pmarker", "m2_pmarker")[apply(!thisHS2D[i,15:16], 1, which)]]
    HSBYCor <- cor(phenoBatchOD[HSBYSegs, thisHS2D[i, "trait"]], gdata[HSBYSegs, thisOneIsTheCisMarker])
    HSRMCor <- cor(phenoBatchOD[HSRMSegs, thisHS2D[i, "trait"]], gdata[HSRMSegs, thisOneIsTheCisMarker])
    c(HSBYCor, HSRMCor, thisOneIsTheCisMarker)
})), stringsAsFactors=FALSE)
rownames(cisEffectsSplitByHS) <- thisHS2D[,"trait"]
cisEffectsSplitByHS[,1] <- as.numeric(cisEffectsSplitByHS[,1])
cisEffectsSplitByHS[,2] <- as.numeric(cisEffectsSplitByHS[,2])
t.test(abs(cisEffectsSplitByHS[,1]), abs(cisEffectsSplitByHS[,2]))
# p = 0.86
# nope
# KRE33: p = 0.6
# MKT1: 0.3
# GAP1: 0.6
# IRA2: 0.7
cor.test(cisEffectsSplitByHS[,1], cisEffectsSplitByHS[,2])
# 0.95, p=8e-8
# no systematic difference for the cis effect sizes depending on HAP1
# and for none of the others either

# No evidence for boosted cis effects by active trans

# just checked for all-pairwise:
# HAP1 0.5, KR33 0.8, MKT1 0.6, GAP1 0.8, IRA2 0.9

# is there any enrichment of cis among these trans targets?
# among interactions?

summary(rowSums(lm2D_from_full[,c("m1_is_cis", "m2_is_cis")])> 0)
# 158 / 387 have a cis. 41%

# cis enrichment?
summary(rowSums(lm2D_between_additive[,c("m1_is_cis", "m2_is_cis")])> 0)
# 300 / 1464 => 20%

# out of...
summary(allPeaksODPad$cis)
# 2969/36,498 => 8%

# for each gene with an interaction, draw pairs at random
cisPerms <- sapply(1:100, function(j){
    t(sapply(unique(lm2D_between_additive$trait), function(x){
    trueInteractions <- length(which(lm2D_between_additive$trait == x))
    trueCisInteractions <- length(which(lm2D_between_additive$trait == x & rowSums(lm2D_between_additive[,c("m1_is_cis", "m2_is_cis")]) > 0))
    randomCis <- sapply(1:trueInteractions, function(i){
        theseTwo <- sample(1:length(which(allPeaksODPad$gene == x)), 2)
        length(which(allPeaksODPad[which(allPeaksODPad$gene == x),"cis"][theseTwo]))
    })
    c(trueInteractions, trueCisInteractions, sum(randomCis))
}))[,3]
})
# some genes do have multiple pairs!
summary(apply(cisPerms, 2, sum))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#164.0   184.0   192.5   192.1   201.2   223.0
# so, 0/100 reaches 300, not even close

# yes, there is enrichment of cis


###############################################
# how many of the 2D are also found in the between additive?
# NOT elegant code

# make a reference out of the between-additive:

tableToGRanges <- function(dat, thisPad=100000){
    m1Splits <- lapply(dat[,"m1_pmarker"], function(x){strsplit(x, ":")[[1]]})
    m1Chrs <- sapply(m1Splits, function(x){x[1]})
    m1Pos <- sapply(m1Splits, function(x){as.integer(strsplit(x[2], "_")[[1]][1])})
    m1GRange <- GRanges(seqnames = paste(dat[,"trait"], m1Chrs, sep="_"), ranges = IRanges(start = m1Pos - thisPad, end= m1Pos + thisPad))

    m2Splits <- lapply(dat[,"m2_pmarker"], function(x){strsplit(x, ":")})
    m2Chrs <- sapply(m2Splits, function(x){x[[1]][1]})
    m2Pos <- sapply(m2Splits, function(x){as.integer(strsplit(x[[1]][2], "_")[[1]][1])})
    m2GRange <- GRanges(seqnames = paste(dat[,"trait"], m2Chrs, sep="_"), ranges = IRanges(start = m2Pos - thisPad, end= m2Pos + thisPad))

    list(m1GRange, m2GRange)
}

full2DGRanges <- tableToGRanges(lm2D_from_full)
betweenAdditiveGRanges <- tableToGRanges(lm2D_between_additive)

full2DInBetweenAdditive <- sapply(1:length(full2DGRanges[[1]]), function(i){
    m1Hits1 <- subjectHits(findOverlaps(full2DGRanges[[1]][i], betweenAdditiveGRanges[[1]]))
    m1Hits2 <- subjectHits(findOverlaps(full2DGRanges[[1]][i], betweenAdditiveGRanges[[2]]))
    m2Hits1 <- subjectHits(findOverlaps(full2DGRanges[[2]][i], betweenAdditiveGRanges[[1]]))
    m2Hits2 <- subjectHits(findOverlaps(full2DGRanges[[2]][i], betweenAdditiveGRanges[[2]]))

    print(i)
    return(length(intersect(m1Hits1, m2Hits2)) > 0 | length(intersect(m1Hits2, m2Hits1)) > 0)
})

summary(full2DInBetweenAdditive)
# 259 / 387
# with 50k pad: 272
# 100k pad: 276

# this could be faster, but keeps the code consistent with above
full2DInBetweenAdditiveAtLeastOne <- sapply(1:length(full2DGRanges[[1]]), function(i){
    m1Hits1 <- subjectHits(findOverlaps(full2DGRanges[[1]][i], betweenAdditiveGRanges[[1]]))
    m1Hits2 <- subjectHits(findOverlaps(full2DGRanges[[1]][i], betweenAdditiveGRanges[[2]]))
    m2Hits1 <- subjectHits(findOverlaps(full2DGRanges[[2]][i], betweenAdditiveGRanges[[1]]))
    m2Hits2 <- subjectHits(findOverlaps(full2DGRanges[[2]][i], betweenAdditiveGRanges[[2]]))
    
    print(i)
    return(length(m1Hits1) > 0 | length(m1Hits2) > 0 | length(m2Hits1) > 0 | length(m2Hits2) > 0)
})

summary(full2DInBetweenAdditiveAtLeastOne)
# 25 pad: 301
# 50k pad: 303; 78%


