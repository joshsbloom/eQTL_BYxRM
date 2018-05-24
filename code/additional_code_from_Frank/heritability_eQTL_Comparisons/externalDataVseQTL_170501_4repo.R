library(car)
library(topGO)

geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

load("R_externalGeneData_170501.RData")
load("log2_t.tpm.matrix.RData")
load("R_allPeaksODPad_161213.RData")
load("vcA.unscaled.RData")

h2 <- vcA.OD.unscaled[,1] / rowSums(vcA.OD.unscaled)
VTotal <- rowSums(vcA.OD.unscaled)


localeQTL <- allPeaksODPad[allPeaksODPad$cis,]
localeQTL <- cbind(localeQTL, table(localeQTL$gene)[localeQTL$gene])
colnames(localeQTL)[ncol(localeQTL)] <- "Freq"
throwOutDecider <- rownames(localeQTL[localeQTL$Freq > 1,])
throwOutDeciderBool <- sapply(throwOutDecider, function(x){
    theseeQTL <- localeQTL[localeQTL$gene == localeQTL[x,]$gene, ]
    theseDists <- abs(theseeQTL$marker.gcoord - theseeQTL$gene.gcoord)
    abs(localeQTL[x,]$marker.gcoord - localeQTL[x,]$gene.gcoord) != min(theseDists)
})
localeQTL <- localeQTL[!rownames(localeQTL) %in% names(throwOutDeciderBool[throwOutDeciderBool]),]



# distribution of peaks per gene
# cannot just "table" since that would miss genes with no eQTL
# also cannot use the external gene annotation since that would include genes we didn't find as expressed
eQTLPerGene <- sapply(colnames(t.tpm.matrix), function(x){
    length(which(allPeaksODPad$gene == x))
})
summary(eQTLPerGene)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.000   4.000   6.000   6.381   8.000  21.000


pdf("eQTLPerGene.pdf", width=8, height=5)
barplot(as.numeric(table(eQTLPerGene)), ylab="number of genes", xlab="eQTL per gene", names.arg=names(table(eQTLPerGene)))
dev.off()

# how many have none?
length(which(eQTLPerGene == 0))
# 77

# distribution for > 0
summary(eQTLPerGene[which(eQTLPerGene > 0)])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   4.000   6.000   6.468   8.000  21.000


# which are the genes with most eQTL?
eQTLPerGene[which(eQTLPerGene == max(eQTLPerGene))]
#YGL038C YIL165C YIL164C YLL053C YLL052C
#21      21      21      21      21

allNames[names(eQTLPerGene[which(eQTLPerGene == max(eQTLPerGene))])]
#   YGL038C   YIL165C   YIL164C   YLL053C   YLL052C
# "OCH1" "YIL165C"    "NIT1" "YLL053C"    "AQY2"

externalGeneDataWithQTLData <- cbind(externalGeneData, rownames(externalGeneData) %in% localeQTL$gene, sapply(rownames(externalGeneData), function(x){ret=NA; if(x %in% localeQTL$gene){ret <- localeQTL[localeQTL$gene == x, "lm.coeff"]}; return(ret)}), eQTLPerGene[rownames(externalGeneData)], h2[rownames(externalGeneData)], VTotal[rownames(externalGeneData)], vcA.OD.unscaled[,1][rownames(externalGeneData)], vcA.OD.unscaled[,2][rownames(externalGeneData)])
colnames(externalGeneDataWithQTLData)[12:18] <- c("hasLocaleQTL", "localeQTL_lmCoeff", "numbereQTL", "h2", "VTotal", "A", "E")


summary(lm(h2 ~ VTotal + expression + essential + dNdS + PPI + GxGStrict + isTF + hasHumanHomolog + hasParalog, data=externalGeneDataWithQTLData))
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)          2.145e-01  1.523e-02  14.088  < 2e-16 ***
#VTotal               7.476e-02  4.655e-03  16.059  < 2e-16 ***
#expression           2.089e-02  1.992e-03  10.485  < 2e-16 ***
#essentialTRUE       -7.923e-02  8.358e-03  -9.479  < 2e-16 ***
#dNdS                -2.304e-01  4.931e-02  -4.672 3.13e-06 ***
#PPI                 -1.188e-04  4.045e-05  -2.937 0.003342 **
#GxGStrict           -9.173e-05  2.595e-05  -3.535 0.000414 ***
#isTFTRUE             1.166e-04  1.552e-02   0.008 0.994008
#hasHumanHomologTRUE -6.856e-03  6.305e-03  -1.087 0.276964
#hasParalogTRUE       2.160e-02  7.561e-03   2.857 0.004309 **

Anova(lm(h2 ~ VTotal + expression + essential + dNdS + PPI + GxGStrict + isTF + hasHumanHomolog + hasParalog, data=externalGeneDataWithQTLData))
#Anova Table (Type II tests)))
#Response: h2
#Sum Sq   Df  F value    Pr(>F)
#VTotal           5.776    1 257.9074 < 2.2e-16 ***
#expression       2.462    1 109.9357 < 2.2e-16 ***
#essential        2.013    1  89.8594 < 2.2e-16 ***
#dNdS             0.489    1  21.8282 3.126e-06 ***
#PPI              0.193    1   8.6259 0.0033419 **
#GxGStrict        0.280    1  12.4974 0.0004143 ***
#isTF             0.000    1   0.0001 0.9940080
#hasHumanHomolog  0.026    1   1.1824 0.2769637
#hasParalog       0.183    1   8.1623 0.0043093 **
#Residuals       61.321 2738

# note identical p-values (so both are type II anovas), and also gives F, SS, residual




# how does h2 relate to eQTL number etc?
pdf("h2_vs_Others.pdf", width=16, height=8)
par(mfcol=c(2,4))
plot(vcA.OD.unscaled[,"E"], vcA.OD.unscaled[,"A"], pch=19, col="#00000022")
plot(vcA.OD.unscaled[,"E"], vcA.OD.unscaled[,"A"], pch=19, col="#00000022", xlim=c(0, 1), ylim=c(0, 1))

plot(externalGeneDataWithQTLData$VTotal, externalGeneDataWithQTLData$h2, pch=19, col="#00000022")
plot(externalGeneDataWithQTLData$VTotal, externalGeneDataWithQTLData$h2, pch=19, col="#00000022", xlim=c(0, 1))

plot(externalGeneDataWithQTLData$VTotal, externalGeneDataWithQTLData$numbereQTL, pch=19, col="#00000022")
plot(externalGeneDataWithQTLData$VTotal, externalGeneDataWithQTLData$numbereQTL, pch=19, col="#00000022", xlim=c(0, 1))

plot(externalGeneDataWithQTLData$h2, externalGeneDataWithQTLData$numbereQTL, pch=19, col="#00000022")
points()
dev.off()

cor.test(externalGeneDataWithQTLData$h2, externalGeneDataWithQTLData$numbereQTL, method="s")
# r = 0.56, rho = 0.66, p < 2.2e-16
# overall positive

cor.test(externalGeneDataWithQTLData$h2[which(externalGeneDataWithQTLData$h2 >= 0.6)], externalGeneDataWithQTLData$numbereQTL[which(externalGeneDataWithQTLData$h2 >= 0.6)], method="p")
# but for higher h2, actually negative!
# r = -0.47, rho=-0.45, p < 2.2e-16
# this is 423 genes, i.e. 7.4% of 5720

cor.test(externalGeneDataWithQTLData$h2[which(externalGeneDataWithQTLData$h2 < 0.6)], externalGeneDataWithQTLData$numbereQTL[which(externalGeneDataWithQTLData$h2 < 0.6)], method="p")
# r = 0.70, p < 2.2e-16

# are there any with high h2 and no eQTL?
rownames(externalGeneDataWithQTLData)[which(externalGeneDataWithQTLData$h2 > 0.9 & externalGeneDataWithQTLData$numbereQTL == 1),]
# none with zero
# 7 with just one eQTL: "YDL227C"   "YDR038C"   "YER109C"   "YGL053W"   "YHL048C-A" "YJL217W" "YKR104W"
 allPeaksODPad[allPeaksODPad$gene %in% rownames(externalGeneDataWithQTLData)[which(externalGeneDataWithQTLData$h2 > 0.9 & externalGeneDataWithQTLData$numbereQTL == 1)],]
# only one classified as cis
# but in fact, all seem to be close to their gene:
#gene           pmarker pcind          r      LOD
#chrIV.263     YDL227C   chrIV:46090_C/G    19 -0.9608724 564.2244 <= HO. classified as cis
#chrIV.436     YDR038C  chrIV:538866_C/T   537 -0.9125001 392.8538 <= ENA5 (i.e. complex locus)!; marker & gene coords 1898888   1887444; long stretch with no markers; no reads here in footprints; not included in ASE data
#chrV.245      YER109C   chrV:374819_A/C   353 -0.9372776 463.1851 <= FLO8; marker just outside of 200bp downstream interval; mRNA but no footprints; only in BY => missing in RM? does not look like a healthy gene; S288C known to contain premature stop (blocked reading frame), but why nothing in RM?; not in ASE data
#chrV.347    YHL048C-A    chrV:10202_A/T     1 -0.9444324 488.9910 <= at chr start; not in ASE data; small putative ORF; marker is really far away, probably at chr end
#chrVII.1259   YGL053W chrVII:405567_G/T   425 -0.9411609 476.7897 <= PRM8, next to OLE1 & a Ty element; no reads go here in RM & gene looks broken in BY; marker is 3kb downstream, perhaps due to mapping problems here?; not in ASE data
#chrX.598      YJL217W    chrX:26108_C/T     1 -0.9460741 495.3960 <= REE1, pretty close to chr end; outside of linkage map; not in ASE data
#chrXI.740     YKR104W  chrXI:643662_G/A   820 -0.9478398 502.5125 <= chr end; joint protein with neighbor YKR103W in other strains; no FP reads, not in ASE data; outside of linkage map

# are some of the very high h2 cases genes that aren't present in both strains?

# other genes with very high h2

allPeaksODPad[allPeaksODPad$gene %in% rownames(externalGeneDataWithQTLData)[which(externalGeneDataWithQTLData$h2 > 0.95)],]
# highest h2 is 0.9542
# identical for YCL066W (HMLALPHA1) & YCR040W (MATALPHA1)
# total of five genes with h2 > 0.95
h2[which(h2 > 0.95)]
#YCL066W   YCR040W   YDL227C   YLR040C   YPL187W
#0.9542172 0.9542172 0.9517206 0.9520537 0.9532949
# first two, s. above
# YDL227C = HO
# YLR040C = AFB1 "A-factor barrier"
# YPL187W  = MF(ALPHA)1
# clearly mating type is highly heritable

# there is a correlation between # eQTL and h2
# but the highest h2 do NOT have many eQTL:
summary(externalGeneDataWithQTLData$numbereQTL[which(externalGeneDataWithQTLData$h2 > 0.9)])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.000   1.750   2.500   3.036   4.000   8.000
# 28 genes


pdf("eQTLNumberVsExternal_170403.pdf", width=9, height=9)
par(mfrow=c(3, 3))
for (j in c("numbereQTL", "h2",  "VTotal")){
    for(i in c("expression", "dNdS", "PPI", "GxGLenient", "GxGStrict")){
        plot(externalGeneDataWithQTLData[, i], externalGeneDataWithQTLData[,j], col="#00000022", pch=19, ylab=j, xlab=i)
        legend("topright", box.lty=0, legend=c(
        paste0("cor = ", round(cor.test(externalGeneDataWithQTLData[,j], externalGeneDataWithQTLData[,i], method="s")$est, 2)),
        paste0("p = ", cor.test(externalGeneDataWithQTLData[,j], externalGeneDataWithQTLData[,i], method="s")$p.value))
        )
        if (i == "PPI"){
            plot(externalGeneDataWithQTLData[, i], externalGeneDataWithQTLData[,j], col="#00000022", pch=19, ylab=j, xlab=i, xlim=c(0, 200))
        }
    }
    for(i in c("essential", "isTF", "hasHumanHomolog")){
        boxplot(list(externalGeneDataWithQTLData[,j][which(externalGeneDataWithQTLData[,i])], externalGeneDataWithQTLData[,j][which(!externalGeneDataWithQTLData[,i])]), names=c(i, paste0("not ", i)), ylab=j)
        legend("topright", box.lty=0, legend = paste0("p = ", t.test(externalGeneDataWithQTLData[,j][which(externalGeneDataWithQTLData[,i])], externalGeneDataWithQTLData[,j][which(!externalGeneDataWithQTLData[,i])])$p.value))
    }
}
dev.off()

# for each gene, find the eQTL with the highest VE
highesteQTLVE <- sapply(rownames(externalGeneDataWithQTLData), function(x){
    ret <- allPeaksODPad[allPeaksODPad$gene == x & allPeaksODPad$"var.exp" == max(allPeaksODPad[allPeaksODPad$gene == x, "var.exp"]),"var.exp"]
    if(length(ret) == 0){ret <- NA}
    return(ret)
})

pdf("h2_vs_biggesteQTL.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(externalGeneDataWithQTLData$h2, highesteQTLVE, col="#00000022", pch=19, xlab="heritability", ylab="variance from strongest eQTL")
abline(0, 1, col="grey", lwd=2, lty=2)
plot(externalGeneDataWithQTLData$h2, highesteQTLVE/externalGeneDataWithQTLData$h2, col="#00000022", pch=19, xlab="heritability", ylab="fraction of h2 from strongest eQTL", ylim=c(0, 1.2))
abline(h=1, col="grey", lwd=2, lty=2)
dev.off()

cor.test(externalGeneDataWithQTLData$h2, highesteQTLVE/externalGeneDataWithQTLData$h2)
# r=-0.04, p =  0.006 not the strongest relationship
cor.test(externalGeneDataWithQTLData$h2[externalGeneDataWithQTLData$h2 >=0.6], (highesteQTLVE/externalGeneDataWithQTLData$h2)[externalGeneDataWithQTLData$h2 >=0.6])
# r = 0.56, p < 2.2e16

# at higher heritability, more of it tends to come from a single eQTL
cor.test(externalGeneDataWithQTLData$h2[externalGeneDataWithQTLData$h2 < 0.6], (highesteQTLVE/externalGeneDataWithQTLData$h2)[externalGeneDataWithQTLData$h2 < 0.6])
# r = -0.04, p = 0.001






#################################
# GO on genes with high/low h2

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

testGenes <- externalGeneDataWithQTLData$h2
names(testGenes) <- rownames(externalGeneDataWithQTLData)

testGenes <- testGenes[!is.na(testGenes)]

# there has to be some function to build the topGOdata object
dummyFun <- function(x){x > 0.5}
# confirmed that the threshold in this function doesn't matter for Ttest
# the "expected" and "significant" columns in the result file change, but the p-values and ranking remain the same

# by default, the score is treated as "increasing", i.e. low is better
# CAREFUL WITH THE FILE NAMES!!!
pdf("topGO/h2/h2_increasing.pdf", width=15, height=15)
for (thisOntology in c("BP", "MF", "CC")){
    try({
        GOData <- new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSelectionFun=dummyFun, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
        GOresult <- runTest(GOData, algorithm="classic", statistic="t", scoreOrder = "increasing")
        plotGOToTree(GOData, GOresult, sigThres = 1.5e-5)
        
        write.table(GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult))), file=paste0("topGO/h2/h2GOresult_", thisOntology, "_increasing.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    })
}
dev.off()

# the increasing and decreasing lists are not the same
# they appear to be in reverse order
# so one is enrichment for genes with many (vs all) eQTL, the other is for genes with few eQTL
# to tell which is which, look at the top category in the "increasing" file: GO:0022613 (ribonucleoprotein comlex biogenesis) and plot the scores:
# best decreasing: GO:0055114, oxidation-reduction process

testGOData <- new("topGOdata", ontology="BP", allGenes = testGenes, geneSelectionFun=dummyFun, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
pdf("GO:0055114.pdf")
showGroupDensity(testGOData, "GO:0055114", ranks = FALSE)
dev.off()
# The top GO from the DEcreasing category (GO:0055114) has HIGHER scores.
# the top GO from the INcreasing category has LOWER scores than the universe.


# can also run a version where we look only at the highest/lowest h2 genes and do Fisher's exact test

top5Percent <- function(x){x > quantile(x, probs=0.95)}
bottom5Percent <- function(x){x < quantile(x, probs=0.05)}

highH2 <- function(x){x >= 0.9} # 28 genes
lowH2 <- function(x){x <= 0.01} # 39 genes


# by default, the score is treated as "increasing"
# CAREFUL WITH THE FILE NAMES - swap them out according to which gene set and test is being done
#pdf("topGO/h2/h2_atLeast_0.9.pdf", width=15, height=15)
pdf("topGO/h2/h2_atMost_0.01.pdf", width=15, height=15)
#pdf("topGO/h2/h2_bottom5Percent.pdf", width=15, height=15)
for (thisOntology in c("BP", "MF", "CC")){
    try({
        GOData <- new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSelectionFun=lowH2, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
        #GOData <- new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSelectionFun=bottom5Percent, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
        GOresult <- runTest(GOData, algorithm="classic", statistic="fisher")
        plotGOToTree(GOData, GOresult, sigThres = 1.5e-5)
        
        write.table(GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult))), file=paste0("topGO/h2/h2_", thisOntology, "_atMost_0.01.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
        #write.table(GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult))), file=paste0("topGO/h2/h2_", thisOntology, "_bottom5Percent.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    })
}
dev.off()



