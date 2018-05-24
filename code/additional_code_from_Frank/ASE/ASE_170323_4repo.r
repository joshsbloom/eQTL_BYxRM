library(GenomicRanges)
library(qvalue)
library(smatr)

##################
# compare to ASE

# data for Albert & Noori ASE data
load("R_jointASEData.RData")
# eQTL table
load("R_allPeaksODPad_161213.RData")
# peak-based effect sizes for eQTL, calculated on uncorrected data to preserve "natural" effect sizes
load("peakModel.RData")

# for pulling LODs at individual genes and genome positions
load("scanoneLODS_OD_stranded.RData") # file is ~1Gb
load("marker.LD.RData") # file is 130Mb

# need to retstrict marker.LD to the set that is present in scanoneLODS
marker.LD.subset <- lapply(marker.LD, function(x){x[rownames(x) %in% colnames(scanoneLODS.OD[[2]]), colnames(x) %in% colnames(scanoneLODS.OD[[2]])]})

# load gene info to make vectors relating systematic to common names
geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

localeQTL <- allPeaksODPad[allPeaksODPad$cis,]
ncol(localeQTL)
# 2969

# if there are more than 1 local eQTL, keep the one closer to the gene
# for each gene, how often is it there:
# table(localeQTL$gene)

# run the filter
localeQTL <- cbind(localeQTL, table(localeQTL$gene)[localeQTL$gene])
colnames(localeQTL)[ncol(localeQTL)] <- "Freq"
throwOutDecider <- rownames(localeQTL[localeQTL$Freq > 1,])
throwOutDeciderBool <- sapply(throwOutDecider, function(x){
    theseeQTL <- localeQTL[localeQTL$gene == localeQTL[x,]$gene, ]
    theseDists <- abs(theseeQTL$marker.gcoord - theseeQTL$gene.gcoord)
    abs(localeQTL[x,]$marker.gcoord - localeQTL[x,]$gene.gcoord) != min(theseDists)
})
localeQTL <- localeQTL[!rownames(localeQTL) %in% names(throwOutDeciderBool[throwOutDeciderBool]),]
# 2969 going in, 2884 coming out


# the peak table above uses lm.coeffs from scaled data
# add lm.coeffs from peakModel that was fit using the covariates in the same model
# these should be on the proper, "natural" scale for comparison to ASE
naturalLMCoeff <- sapply(localeQTL$gene, function(x){
    thisPM <- peakModel[[x]]
    thisRowName <- rownames(localeQTL)[which(localeQTL$gene == x)]
    print(thisRowName)
    thisPM[thisRowName, "lm.coeff"]
})
# sanity check
cor.test(naturalLMCoeff, localeQTL$"lm.coeff")
# r = 0.69, rho = 0.94

pdf("ASE_lmCoeffsScaledVsNatural.pdf")
plot(naturalLMCoeff, localeQTL$"lm.coeff", xlab="natural scale", ylab="scaled", col="#00000022")
abline(0, 1, col="grey", lwd=2, lty=2)
dev.off()
# interesting - the "natural" ones are much more compressed (i.e., mostly smaller), but have a few big outliers

localeQTL <- cbind(localeQTL, naturalLMCoeff)

# width of the local eQTL:
summary(
    sapply(localeQTL$CI.r, function(x){as.numeric(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})
    - sapply(localeQTL$CI.l, function(x){as.numeric(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})
)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0    5054   16208   76858   45540 1338215



############
# some ASE data wrangling

qVals <- apply(jointASEData[,c("NooriB", "AlbertB")], 2, function(x){
    ret <- rep(1, length(x))
    ret[x < 1] <- qvalue(x[x< 1])$q
    ret
})

# make an average ASE between the two ASE datasets
jointFC <- apply(jointASEData, 1, function(x){mean(c(x[1], x[2]))})

# make a decider for significance (Bonferroni)
sigDecider <- rep(0, nrow(jointASEData))
names(sigDecider) <- rownames(jointASEData)
sigDecider[jointASEData[,"NooriB"] < 0.05/nrow(jointASEData) | jointASEData[,"AlbertB"] < 0.05/nrow(jointASEData)] <- 1
sigDecider[jointASEData[,"NooriB"] < 0.05/nrow(jointASEData) & jointASEData[,"AlbertB"] < 0.05/nrow(jointASEData)] <- 2

# FDR-based threshold (not used in paper)
sigDeciderFDR <- rep(0, nrow(jointASEData))
names(sigDeciderFDR) <- rownames(jointASEData)
sigDeciderFDR[qVals[,"NooriB"] < 0.05 | qVals[,"AlbertB"] < 0.05] <- 1
sigDeciderFDR[qVals[,"NooriB"] < 0.05 & qVals[,"AlbertB"] < 0.05] <- 2


# append the ASE info to the eQTL data

# !!!****NOTE that we need to invert the ASE sign to make consistent with eQTL data***!!!
# previous ASE data was polarized as "positive = higher in BY"
# current eQTL are other way around

localeQTL <- cbind(localeQTL, -jointFC[localeQTL$gene], sigDecider[localeQTL$gene], sigDeciderFDR[localeQTL$gene])
colnames(localeQTL)[(ncol(localeQTL)-2): ncol(localeQTL)] <- c("ASE", "sigDeciderBonferroni", "sigDeciderFDR")

#save(localeQTL, file="R_localeQTL_170222.RData")



############
# load the previously saved table here, as the basis for the remaining analyses
load("R_localeQTL_170222.RData")
rownames(localeQTL) <- localeQTL$gene

# for every local eQTL, pull the ASE signal

pdf("ASE_vs_eQTLFC.pdf", width=12, height=12)
plotCol <- rep("#00000022", nrow(localeQTL))
plotCol[localeQTL$gene %in% names(sigDecider)[sigDecider == 1]] <- "#FF000022"
plotCol[localeQTL$gene %in% names(sigDecider)[sigDecider == 2]] <- "red"
names(plotCol) <- localeQTL$gene

par(mfrow=c(2,2))

plot(localeQTL$ASE, localeQTL$"naturalLMCoeff" * 2, pch=19, col=plotCol, ylab="local eQTL fold change", xlab="log2(ASE signal)", main="significant eQTL")
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")

plot(localeQTL$ASE, localeQTL$"naturalLMCoeff" * 2, pch=19, col=plotCol, ylab="local eQTL fold change", xlab="log2(ASE signal)", main="significant eQTL", xlim=c(-2, 2), ylim=c(-2, 2))
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")

plot(-jointFC[localeQTL$gene], localeQTL$"naturalLMCoeff" * 2, pch=19, col=plotCol, ylab="local eQTL fold change", xlab="log2(ASE signal)", type="n", main="significant ASE and eQTL")
for(i in localeQTL$gene){
    if (i %in% rownames(jointASEData) & sigDecider[i] > 0){
        lines(-jointASEData[i, c("NooriFC", "AlbertFC")], rep(localeQTL[localeQTL$gene == i, "naturalLMCoeff"] * 2, 2), col=plotCol[i])
    }
}
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")

plot(-jointFC[localeQTL$gene], localeQTL$"naturalLMCoeff" * 2, pch=19, col=plotCol, ylab="local eQTL fold change", xlab="log2(ASE signal)", type="n", main="significant ASE and eQTL")
for(i in localeQTL$gene){
    if (i %in% rownames(jointASEData)){
        lines(-jointASEData[i, c("NooriFC", "AlbertFC")], rep(localeQTL[localeQTL$gene == i, "naturalLMCoeff"] * 2, 2), col=plotCol[i])
    }
}
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")

dev.off()

# SMA analysis
xForLM <- -jointFC[localeQTL$gene]
yForLM <- localeQTL$"naturalLMCoeff" * 2
sma(yForLM ~ xForLM, na.action="na.omit")
#slope=0.95 (0.91-0.98)


# overall numbers for agreement:
# how many genes have ASE?
length(which(sigDecider == 2))
# Bonferroni: 121 both, 598 either

# how many genes have a local eQTL and are included in the ASE dataset?
nrow(localeQTL[!is.na(localeQTL$sigDeciderBonferroni),])
# 1974

# how many local eQTL have ASE in at least one (both) datasets?
nrow(localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni > 0,])
# Bonferroni: > 0: 451, 2: 100

# how many local eQTL have NO ASE in at least one (both) datasets?
nrow(localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == 0,])
# Bonf: 1523


# are sig ASE and local eQTL enriched:
# the universe in these tests is the genes with measured ASE
fisher.test(cbind(
    c(length(names(sigDecider)[sigDecider > 0 & names(sigDecider) %in% localeQTL$gene]),
        length(names(sigDecider)[sigDecider == 0 & names(sigDecider) %in% localeQTL$gene])),
    c(length(names(sigDecider)[sigDecider > 0 & !(names(sigDecider) %in% localeQTL$gene)]),
        length(names(sigDecider)[sigDecider == 0 & !(names(sigDecider) %in% localeQTL$gene)])))
)
# Bonf: p < 2.2e-16, odds = 2.5

fisher.test(cbind(
    c(length(names(sigDecider)[sigDecider == 2 & names(sigDecider) %in% localeQTL$gene]),
        length(names(sigDecider)[sigDecider < 2 & names(sigDecider) %in% localeQTL$gene])),
    c(length(names(sigDecider)[sigDecider == 2 & !(names(sigDecider) %in% localeQTL$gene)]),
        length(names(sigDecider)[sigDecider < 2 & !(names(sigDecider) %in% localeQTL$gene)])))
)
# p = 2.5e-8, odds = 3.4


# do eQTL with/without ASE differ in effect size?
wilcox.test(
abs(2 * localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == 1, "naturalLMCoeff"]),
abs(2 * localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == 2, "naturalLMCoeff"])
)
# 0 vs >=0, 0 vs 1, 0 vs 2: p<2.2e-16
# 2 vs 1: p = 5e-7


pdf("eQTLFCByASESignificance.pdf")
boxplot(list(
    abs(2 * localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == 0, "naturalLMCoeff"]),
    abs(2 * localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == 1, "naturalLMCoeff"]),
    abs(2 * localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == 2, "naturalLMCoeff"])
), names=0:2, xlab="number of significant ASE datasets", ylab="absolute eQTL log2(fold change)", main="", ylim=c(0,3))
for(i in 0:2){
    points(jitter(rep(i+1, length(which(!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == i))), amount=0.3), abs(2 * localeQTL[!is.na(localeQTL$sigDeciderBonferroni) & localeQTL$sigDeciderBonferroni == i, "naturalLMCoeff"]), col="#0000FF22", pch=19, cex=0.5)
}
dev.off()





#################################
# plot all ASE with eQTL r, irrespective of eQTL sig


# use Josh's recomputed local FCs at all genes (irrespective of local eQTL significance):
load("cisModel.effects.RData")

# compare these to the ones in the eQTL table (which have background QTL)
pdf("ASE_cisEffectComparison_170324.pdf")
plot(localeQTL$naturalLMCoeff, cisModel.effects[localeQTL$gene], col="#00000022")
abline(0, 1, lwd=2, lty=2, col="grey")
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
dev.off()
# there's a handful of cases where the new cis effects are more extreme than the local eQTL
# probably because closeness to a trans eQTL is not considered here


# plot all genes, irrespective of eQTL or ASE sig (but can color in)
allFCs <- cbind(-jointFC, 2*cisModel.effects[names(jointFC)])
colnames(allFCs) <- c("ASE_fc", "eQTL_fc")


plotColsASESig <- rep("green", nrow(allFCs)) # to verify the data selection works. THere should be no green on the plot
plotColsASESig[sigDecider[rownames(allFCs)] == 2] <- "red"
plotColsASESig[sigDecider[rownames(allFCs)] == 1] <- "#FF000022"
names(plotColsASESig) <- rownames(allFCs)

pdf("ASE_all_FCs.pdf", width=12, height=12)
par(mfrow=c(2,2))
plot(allFCs, col="#00000022", xlab="log2(ASE fold change)", ylab="log2(eQTL fold change)", pch=19, main="all genes")
abline(0, 1, lwd=2, lty=2, col="grey")
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")

plot(allFCs[sigDecider[rownames(allFCs)] > 0,], col=plotColsASESig[sigDecider[rownames(allFCs)] > 0], xlab="log2(ASE fold change)", ylab="log2(eQTL fold change)", pch=19, main="genes with significant ASE")
points(allFCs[sigDecider[rownames(allFCs)] > 0 & !rownames(allFCs) %in% localeQTL$gene,])
abline(0, 1, lwd=2, lty=2, col="grey")
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
legend("topleft", box.lty=0,
col=c("red", "#FF000022", "black"),
legend=c("2 significant ASE datasets", "1 significant ASE dataset", "no local eQTL"),
pch=c(19, 19, 1))

textDecider <- (sigDecider[rownames(allFCs)] == 2 & !rownames(allFCs) %in% localeQTL$gene & abs(allFCs[,1]) >= 1) #| (sigDecider[rownames(allFCs)] >0 & !rownames(allFCs) %in% localeQTL$gene & sign(jointASEData[,"NooriFC"])[rownames(allFCs)] == sign(jointASEData[,"AlbertFC"])[rownames(allFCs)] & abs(allFCs[,1]) > 0.5)

text(allFCs[textDecider,], labels = allNames[rownames(allFCs)[textDecider]], pos=2)


plot(allFCs[sigDecider[rownames(allFCs)] > 0,], col=plotColsASESig[sigDecider[rownames(allFCs)] > 0], xlab="log2(ASE fold change)", ylab="log2(eQTL fold change)", pch=19, main="genes with significant ASE")
abline(0, 1, lwd=2, lty=2, col="grey")
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")

for(i in rownames(allFCs)[sigDecider[rownames(allFCs)] > 0]){
    lines(-jointASEData[i,c("AlbertFC", "NooriFC")], rep(allFCs[i, 2], 2), col=plotColsASESig[i])
}
points(allFCs[sigDecider[rownames(allFCs)] > 0 & !rownames(allFCs) %in% localeQTL$gene,])

dev.off()




###################
# prepare table for paper

paralogInfo <- read.table("paralogs_ygob_2012.txt", stringsAsFactors=FALSE, head=TRUE, sep="\t", quote="\"")


# extend the "allFC" object with eQTL info
# need to look up local LODs
# at every gene, pull the closest marker and calculate a fold change
getClosestMarker <- function(markerDat, thisChr, thisPosition){
    thesePositions <- sapply(colnames(markerDat[[thisChr]]), function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})
    colnames(markerDat[[thisChr]])[which.min(abs(thesePositions - thisPosition))][1]
}

allLocalLODs <- sapply(rownames(allFCs)[geneAnnotation[rownames(allFCs),"chr"] != "chrM"], function(thisGene){
    ret <- NA
    print(thisGene)
    thisClosestMarker <- getClosestMarker(marker.LD.subset, geneAnnotation[thisGene, "chr"], mean(as.numeric(geneAnnotation[thisGene, c("start", "end")])))
    print(thisClosestMarker)
    try({
        ret <- scanoneLODS.OD[[2]][thisGene,thisClosestMarker]
    })
    ret
})



allFCsWithSig <- data.frame(allNames[rownames(allFCs)], allFCs, sigDecider[rownames(allFCs)], allLocalLODs[rownames(allFCs)], jointASEData[rownames(allFCs),1:4])
colnames(allFCsWithSig)[1:5] <- c("geneName", "ASE_fc", "eQTL_fc", "numberSig_ASESets", "local_eQTL_LOD")
# set NA LODs to -1 to avoid having to test for !is.na
allFCsWithSig$eQTL_LOD[is.na(allFCsWithSig$local_eQTL_LOD)] <- -1

# also attach paralogy info
allFCsWithSig <- cbind(allFCsWithSig,
as.numeric(sapply(rownames(allFCsWithSig), function(x){paralogInfo[paralogInfo[,"Gene1"] == x | paralogInfo[,"Gene2"] == x, "AA_percent_id"]})),
    sapply(rownames(allFCsWithSig), function(x){
        ret=NA
        if(x %in% paralogInfo[,"Gene1"]){ret=allNames[paralogInfo[paralogInfo[,"Gene1"] == x, "Gene2"]]}
        if(x %in% paralogInfo[,"Gene2"]){ret=allNames[paralogInfo[paralogInfo[,"Gene2"] == x, "Gene1"]]}
        ret
    })
)
colnames(allFCsWithSig)[(ncol(allFCsWithSig)-1):ncol(allFCsWithSig)] <- c("paralogIdentity", "paralogGene")

save(allFCsWithSig, file="R_allFCsWithSig_170324.RData")

# table for paper
write.table(allFCsWithSig, file="ASE_table4paper_170324.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")




############################################
# variant density in promoter

# load variant numbers per gene
load("R_promoterVariantsPerGene_171020.RData")

pdf("ASE_promoterVariantsVseQTL.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(promoterVariantsPerGene[rownames(allFCs)], abs(allFCs[,1]), xlab="number of promoter variants", ylab="abs(log2(ASE fold change))", col="#00000022")
plot(promoterVariantsPerGene[rownames(allFCs)], abs(allFCs[,2]), xlab="number of promoter variants", ylab="abs(log2(eQTL fold change))", col="#00000022")
dev.off()

cor.test(promoterVariantsPerGene[rownames(allFCs)], abs(allFCs[,1]), method="s")
# rho=0.04, p=0.014

cor.test(promoterVariantsPerGene[rownames(allFCs)], abs(allFCs[,2]), method="s")
# rho=0.23, p < 2.2e-16


