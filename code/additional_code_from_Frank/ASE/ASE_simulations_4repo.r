library(RColorBrewer)
library(VGAM)
library(smatr)

geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]

allNamesInv <- names(allNames)
names(allNamesInv) <- allNames


load("R_jointASEData")
# previously generated in initial ASE comparison:
load("R_localeQTL_170222.RData")
rownames(localeQTL) <- localeQTL$gene

AlbertCounts <- read.table("Albert_2014_hybridCounts.txt", sep="\t", head=TRUE, stringsAsFactors=FALSE)
rownames(AlbertCounts) <- AlbertCounts[,1]
AlbertCounts <- AlbertCounts[,2:3]
summary(rowSums(AlbertCounts))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#20.0    75.0   188.0   535.4   447.0 33170.0

# AlbertCounts has two genes that are not in jointASEdata
# force the two to have the same genes:
AlbertCounts <- AlbertCounts[rownames(jointASEData),]


# Noori data
load("R_NooriCountByGeneDownSampled_170316.RData")
NooriCountByGeneDownSampled <- NooriCountByGeneDownSampled[rownames(jointASEData),]



######################################################
# need to incorporate overdispersion into the simulations
# for that, need to find an overdispersion parameter that matches the observed data

# ASEdat has 3 columns: observed count for one allele, total coverage, "true" fold change
omDbetabinom=function(targetRho, ASEdat) {
    -sum(dbetabinom(ASEdat$AlleleCount, ASEdat$Coverage, ASEdat$Prob, targetRho, log=TRUE))
}


ASEdat4OptimAlbert <- data.frame(
    AlbertCounts[, 1],
    rowSums(AlbertCounts),
    (jointASEData[,"NooriFC"] + jointASEData[,"AlbertFC"] / 2)
)
ASEdat4OptimAlbert[,3] <- 2^ASEdat4OptimAlbert[,3] / ( 1 + 2^ASEdat4OptimAlbert[,3])
names(ASEdat4OptimAlbert) <- c("AlleleCount", "Coverage", "Prob")

optimResultAlbert <- optim(0.5, omDbetabinom, ASEdat = ASEdat4OptimAlbert, method='Brent', lower=1/1e8, upper=1)
optimResultAlbert$par
# 0.005365516


ASEdat4OptimNoori <- data.frame(
    NooriCountByGeneDownSampled[, 1],
    rowSums(NooriCountByGeneDownSampled),
    (jointASEData[,"NooriFC"] + jointASEData[,"AlbertFC"] / 2)
)
ASEdat4OptimNoori[,3] <- 2^ASEdat4OptimNoori[,3] / ( 1 + 2^ASEdat4OptimNoori[,3])
names(ASEdat4OptimNoori) <- c("AlleleCount", "Coverage", "Prob")

optimResultNoori <- optim(0.5, omDbetabinom, ASEdat = ASEdat4OptimNoori, method='Brent', lower=1/1e8, upper=1)
optimResultNoori$par
# 0.00414922
# similar, but not the same

# in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837449/
# they essentially follow our procedure, but treat ALL genes as hvaing expected prob=0.5
# this was also done in Pritchard's WASP: http://www.nature.com/nmeth/journal/v12/n11/extref/nmeth.3582-S1.pdf & http://www.nature.com/nmeth/journal/v12/n11/full/nmeth.3582.html

# so let's try it too:
ASEdat4OptimAll0.5_Albert <- ASEdat4OptimAlbert
ASEdat4OptimAll0.5_Albert$Prob <- 0.5
optimResultAlbert_all0.5 <- optim(0.5, omDbetabinom, ASEdat = ASEdat4OptimAll0.5_Albert, method='Brent', lower=1/1e8, upper=1)
optimResultAlbert_all0.5$par
# 0.01284717
# big, but less than when using the eQTL FCs

ASEdat4OptimAll0.5_Noori <- ASEdat4OptimNoori
ASEdat4OptimAll0.5_Noori$Prob <- 0.5
optimResultNoori_all0.5 <- optim(0.5, omDbetabinom, ASEdat = ASEdat4OptimAll0.5_Noori, method='Brent', lower=1/1e8, upper=1)
optimResultNoori_all0.5$par
# 0.009470975


################################################
# simulate counts without overdispersion

# function to generate paired count data
simulateCountsBinom <- function(thisFC, coverage){
        thisX <- (thisFC * coverage) / (1 + thisFC)
        thisProb = thisX/coverage
        thisRandomCount = rbinom(1, coverage, thisProb)
        c(thisRandomCount, coverage - thisRandomCount)
}

# 2^(localeQTL$"naturalLMCoeff" * 2)[1]
simCountsBinomial <- t(sapply(1:10000, function(i){simulateCountsBinom(0.5, 100)}))

# how well does this recapitulate the observed spread in the ASE data?
oneRandomASESampleAcrossGenesBinom <- function(theseGenes){
    t(sapply(theseGenes, function(x){
        thisFC <- 2^((jointASEData[x, "NooriFC"] + jointASEData[x, "AlbertFC"]) / 2)
        thisCov <- rowSums(AlbertCounts[x,])
        if (is.na(thisCov) | is.na(thisFC)){return(c(NA, NA))}
        simCounts <- simulateCountsBinom(thisFC, thisCov)
        simCounts
    }))
}
oneRandomASESampleAcrossGenesBinom(rownames(jointASEData)[1:3])

# make 100 sets and get the respective FCs
binomialFCs <- sapply(1:100, function(i){
    print(i)
    thisSet <- oneRandomASESampleAcrossGenesBinom(rownames(jointASEData))
    log2(thisSet[,1] / thisSet[,2])
})



################################################
# simulate counts with overdispersion

# function to generate paired count data
# default is for Albert data
simulateCountsBetaBinom <- function(thisFC, coverage, rho=0.005365516){
    thisX <- (thisFC * coverage) / (1 + thisFC)
    thisProb = thisX/coverage
    thisRandomCount = rbetabinom(1, coverage, thisProb, rho)
    c(thisRandomCount, coverage - thisRandomCount)
}
simCountsBetaBinomialAlbert <- t(sapply(1:10000, function(i){simulateCountsBetaBinom(0.5, 100, 0.01284717)}))
simCountsBetaBinomialNoori <- t(sapply(1:10000, function(i){simulateCountsBetaBinom(0.5, 100, 0.009470975)}))

oneRandomASESampleAcrossGenesBetaBinomAlbert <- function(theseGenes){
    t(sapply(theseGenes, function(x){
        thisFC <- 2^((jointASEData[x, "NooriFC"] + jointASEData[x, "AlbertFC"]) / 2)
        thisCov <- rowSums(AlbertCounts[x,])
        if (is.na(thisCov) | is.na(thisFC)){return(c(NA, NA))}
        simCounts <- simulateCountsBetaBinom(thisFC, thisCov, rho=0.01284717)
        simCounts
    }))
}
oneRandomASESampleAcrossGenesBetaBinomNoori <- function(theseGenes){
    t(sapply(theseGenes, function(x){
        thisFC <- 2^((jointASEData[x, "NooriFC"] + jointASEData[x, "AlbertFC"]) / 2)
        thisCov <- sum(NooriCountByGeneDownSampled[x,])
        if (is.na(thisCov) | is.na(thisFC)){return(c(NA, NA))}
        simCounts <- simulateCountsBetaBinom(thisFC, thisCov, rho=0.009470975)
        simCounts
    }))
}

# make 50 sets and get the respective FCs
betaBinomialFCsAlbert <- sapply(1:50, function(i){
    print(i)
    thisSet <- oneRandomASESampleAcrossGenesBetaBinomAlbert(rownames(jointASEData))
    log2(thisSet[,1] / thisSet[,2])
})
betaBinomialFCsNoori <- sapply(1:50, function(i){
    print(i)
    thisSet <- oneRandomASESampleAcrossGenesBetaBinomNoori(rownames(jointASEData))
    log2(thisSet[,1] / thisSet[,2])
})


# compare the three sets of simulated data
pdf("ASE_binomialVsbetaBinomial.pdf")
boxplot(list(simCountsBinomial[,1] / simCountsBinomial[,2], simCountsBetaBinomialAlbert[,1] / simCountsBetaBinomialAlbert[,2], simCountsBetaBinomialNoori[,1] / simCountsBetaBinomialNoori[,2]), names=c("binomial", "betabinomial Albert", "betabinomial Noori"))
dev.off()

summary(simCountsBinomial[,1] / simCountsBinomial[,2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1905  0.4286  0.4925  0.5066  0.5625  1.0410
sd(simCountsBinomial[,1] / simCountsBinomial[,2])
# 0.1089637

summary(simCountsBetaBinomialAlbert[,1] / simCountsBetaBinomialAlbert[,2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.08696 0.38890 0.49250 0.51710 0.61290 1.50000
sd(simCountsBetaBinomialAlbert[,1] / simCountsBetaBinomialAlbert[,2])
# 0.1722746

summary(simCountsBetaBinomialNoori[,1] / simCountsBetaBinomialNoori[,2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1236  0.4085  0.4925  0.5152  0.6129  1.5640
sd(simCountsBetaBinomialNoori[,1] / simCountsBetaBinomialNoori[,2])
# 0.156884

# spread is wider in betabinom, and a tiny little bit wider in Albert vs Noori



# how do the two sets of FCs compare to the observed spread between datasets?
# what is the real correlation between Noori & Albert?
cor.test(jointASEData[,"NooriFC"], jointASEData[,"AlbertFC"])
# r = 0.3909369


# what are correlations between two typical binomial datasets (we'd expect these to be better than Noori/Albert, i.e. too good)
binomialCorrelations <- sapply(1:ncol(binomialFCs), function(i){
    sapply(1:ncol(binomialFCs), function(j){
        if(i >= j){return(NA)}
        cor.test(binomialFCs[,i], binomialFCs[,j])$est
    })
})
summary(as.numeric(binomialCorrelations))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.475   0.528   0.538   0.538   0.549   0.588    5535

# we'd expect the betaBinomials to be worse
# we can directly run a bunch of simulated Albert vs Noori comparisons
betaBinomialCorrelations <- sapply(1:ncol(betaBinomialFCsAlbert), function(i){
    sapply(1:ncol(betaBinomialFCsNoori), function(j){
        if(i >= j){return(NA)}
        cor.test(betaBinomialFCsAlbert[,i], betaBinomialFCsNoori[,j])$est
    })
})
summary(as.numeric(betaBinomialCorrelations))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.3713  0.4058  0.4153  0.4147  0.4240  0.4584    1285
# this is pretty close, and seems to mimick the real construction of this distribution well


# what are the observed fold changes in the local eQTL
testFCs <- quantile(2^(abs(localeQTL$"naturalLMCoeff") * 2), probs=c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
# natural scale:
#0%      25%      50%      75%      90%      95%      99%
#1.020001 1.059405 1.099968 1.187402 1.415266 1.825191 5.841893
# wow, median is only 1.1 fold


# to get power for a different combinations of FC, simulate 1k at a given depth and test
nSims <- 1000
covRangeAlbert <- c(round(quantile(rowSums(AlbertCounts), probs=c(0.1, 0.5, 0.9))), 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
covRangeNoori <- c(round(quantile(rowSums(NooriCountByGeneDownSampled), probs=c(0.1, 0.5, 0.9))), 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
# pick one (Noori had higher depth)
covRange <- covRangeAlbert
#covRange <- covRangeNoori

# beta binomials
# careful with the variable assignments!
simpVals <- lapply(covRange, function(thisCov){
    sapply(testFCs, function(x){
        sapply(1:nSims, function(i){
            simCounts <- simulateCountsBetaBinom(x, thisCov, rho=0.009470975)
            binom.test(simCounts[1], sum(simCounts))$p.value
        })
    })
})
simFCs <- lapply(covRange, function(thisCov){
    sapply(testFCs, function(x){
        sapply(1:nSims, function(i){
            simCounts <- simulateCountsBetaBinom(x, thisCov, rho=0.009470975)
            simCounts[1]/simCounts[2]
        })
    })
})


simPower <- lapply(simpVals, function(y){
    apply(y, 2, function(x){
        c(length(which(x < 0.05)), length(which(x < (0.05/nrow(AlbertCounts)))))/length(x)
    })
})

simDirCorrect <- lapply(simFCs, function(y){
    apply(y, 2, function(x){
        length(which(x > 1))/length(x)
    })
})

#pdf("ASE_powerSim.pdf", width=12, height=4)
#pdf("ASE_powerSim_betaBinomialAlbert.pdf", width=12, height=4)
#pdf("ASE_powerSim_betaBinomialNoori.pdf", width=12, height=4)

par(mfrow=c(1,3))
plotCols = colorRampPalette(brewer.pal(9, "YlOrRd"))(length(covRange))
plotCols[c(1, 3)] <- "grey"
plotCols[2] <- "blue"
lwds=rep(1, length(covRange))
lwds[1:3] <- 3
plotTypes <- rep("l", length(covRange))
plotTypes[1:3] <- "l"
ltys <- rep(2, length(covRange))
ltys[1:3] <- 1

FCHist <- hist((abs(localeQTL$"naturalLMCoeff") * 2), plot=FALSE, breaks=30)

plot(log2(testFCs), simPower[[1]][1,], type="n", ylim=c(0,1), xlim=c(0, log2(max(testFCs))), ylab="power", xlab="log2(fold change)", main="nominal power")
polygon(c(FCHist$mids, rev(FCHist$mids)), c(FCHist$density/max(FCHist$density), rep(0, length(FCHist$mids))), col="lightgreen", border=NA)
for (i in 1:length(simPower)){
    points(log2(testFCs), simPower[[i]][1,], col=plotCols[i], type=plotTypes[i], lwd=lwds[i], lty=ltys[i])
}

plot(log2(testFCs), simPower[[1]][2,], type="n", ylim=c(0,1), xlim=c(0, log2(max(testFCs))), ylab="power", xlab="log2(fold change)", main="Bonferroni power")
polygon(c(FCHist$mids, rev(FCHist$mids)), c(FCHist$density/max(FCHist$density), rep(0, length(FCHist$mids))), col="lightgreen", border=NA)
for (i in 1:length(simPower)){
    points(log2(testFCs), simPower[[i]][2,], col=plotCols[i], type=plotTypes[i], lwd=lwds[i], lty=ltys[i])
}

plot(log2(testFCs), simPower[[1]][2,], type="n", ylim=c(0,1), xlim=c(0, log2(max(testFCs))), ylab="fraction correct", xlab="log2(fold change)", main="directional agreement")
polygon(c(FCHist$mids, rev(FCHist$mids)), c(FCHist$density/max(FCHist$density), rep(0, length(FCHist$mids))), col="lightgreen", border=NA)
for (i in 1:length(simPower)){
    points(log2(testFCs), simDirCorrect[[i]], col=plotCols[i], type=plotTypes[i], lwd=lwds[i], lty=ltys[i])
}

legend("bottomright", box.lty=0, col=plotCols, lwd=lwds, lty=ltys, legend=c("10th percentile observed", "median observed", "90th percentile observed", covRange[4:length(covRange)]))

dev.off()

# for the median FC, power with Bonferroni reaches ≥80% at 1k coverage
# but the direction of effect is already correctly estimated 50% of time at 10X
# ≥90% at 100X







###################################
# calculate power for each local eQTL, with coverage specific for that gene

# first, restrict eQTL to det that are in ASE data
localeQTLForASE <- localeQTL[localeQTL$gene %in% rownames(jointASEData),]

nSimByGene <- 100
powerByGeneAlbert <- t(apply(localeQTLForASE, 1, function(x){
    print(x["gene"])
    thisFC <- 2^(as.numeric(x["naturalLMCoeff"])*2)
    thisCov <- rowSums(AlbertCounts)[x["gene"]]
    if (is.na(thisCov) | is.na(thisFC)){return(c(NA, NA))}
    thesePs <- sapply(1:nSimByGene, function(i){
        simCounts <- simulateCountsBetaBinom(thisFC, thisCov, rho=0.01284717)
        binom.test(simCounts[1], sum(simCounts))$p.value
    })
    c(length(which(thesePs < 0.05)), length(which(thesePs < (0.05/nrow(AlbertCounts)))))/length(thesePs)
}))

powerByGeneNoori <- t(apply(localeQTLForASE, 1, function(x){
    print(x["gene"])
    thisFC <- 2^(as.numeric(x["naturalLMCoeff"])*2)
    thisCov <- rowSums(NooriCountByGeneDownSampled)[x["gene"]]
    if (is.na(thisCov) | is.na(thisFC)){return(c(NA, NA))}
    thesePs <- sapply(1:nSimByGene, function(i){
        simCounts <- simulateCountsBetaBinom(thisFC, thisCov, rho=0.009470975)
        binom.test(simCounts[1], sum(simCounts))$p.value
    })
    c(length(which(thesePs < 0.05)), length(which(thesePs < (0.05/nrow(NooriCountByGeneDownSampled)))))/length(thesePs)
}))


# how many can be detected with Bonferroni?
length(which(!is.na(powerByGeneAlbert[,2])))
# beta (Albert & Noori): 1974
length(which(powerByGeneAlbert[,2] > 0.8 & !is.na(powerByGeneAlbert[,2])) )

# beta-binomial Albert:
# power > 0.5: 138
# power > 0.8: 36

length(which(powerByGeneNoori[,2] > 0.5 & !is.na(powerByGeneNoori[,2])) )
# beta-binomial Noori:
# power > 0.5: 541
# power > 0.8: 132
# higher power in Noori, probably because of higher counts there?

# and how many are powered in both?
length(which(powerByGeneAlbert[,2] > 0.8 & !is.na(powerByGeneAlbert[,2]) & powerByGeneNoori[,2] > 0.8 & !is.na(powerByGeneNoori[,2])) )
# 0.5: 131
# 0.8: 34
# so the well powered Albert genes are mostly contained in the powered Noori genes


# for comparing to real data, make eQTL and ASE lists the same
jointASEDataForeQTL <- jointASEData[localeQTLForASE$gene,]

# how many genes DO we see significant at Bonferroni?
bonfCut <- 0.05/nrow(jointASEData)
length(which(jointASEDataForeQTL[,"NooriB"] < bonfCut))
# 406
length(which(jointASEDataForeQTL[,"AlbertB"] < bonfCut))
# 145
length(which(jointASEDataForeQTL[,"AlbertB"] < bonfCut & jointASEDataForeQTL[,"NooriB"] < bonfCut))
# 100


directionByGeneAlbert <- apply(localeQTLForASE, 1, function(x){
    thisFC <- 2^(as.numeric(x["naturalLMCoeff"])*2)
    thisCov <- rowSums(AlbertCounts)[x["gene"]]
    if (is.na(thisCov) | is.na(thisFC)){return(NA)}
    theseDirs <- sapply(1:nSimByGene, function(i){
        simCounts <- simulateCountsBetaBinom(thisFC, thisCov, rho=0.01284717)
        sign(log2(simCounts[1]/simCounts[2])) == sign(log2(thisFC))
    })
    length(which(theseDirs))/length(theseDirs)
})

directionByGeneNoori <- apply(localeQTLForASE, 1, function(x){
    thisFC <- 2^(as.numeric(x["naturalLMCoeff"])*2)
    thisCov <- rowSums(NooriCountByGeneDownSampled)[x["gene"]]
    if (is.na(thisCov) | is.na(thisFC)){return(NA)}
    theseDirs <- sapply(1:nSimByGene, function(i){
        simCounts <- simulateCountsBetaBinom(thisFC, thisCov, rho=0.009470975)
        sign(log2(simCounts[1]/simCounts[2])) == sign(log2(thisFC))
    })
    length(which(theseDirs))/length(theseDirs)
})


# how many are majority same direction?
length(which(!is.na(directionByGeneAlbert) & directionByGeneAlbert > 0.8))
# beta binomial Albert:
# 0.8: 233

length(which(!is.na(directionByGeneNoori) & directionByGeneNoori > 0.8))
# 0.8: 392


# OK, the above is for numbers of genes.
# what about gene identities: how many of the individual powered genes DO we see as significant?
# attach power to local eQTL

localeQTLForASEWithPower <- cbind(localeQTLForASE, powerByGeneNoori, powerByGeneAlbert, directionByGeneNoori, directionByGeneAlbert, jointASEDataForeQTL[,1:4], NooriCountByGeneDownSampled[localeQTLForASE$gene,], AlbertCounts[localeQTLForASE$gene,])
colnames(localeQTLForASEWithPower)[21:24] <- c("NooriPowerNominal", "NooriPowerBonferroni", "AlbertPowerNominal", "AlbertPowerBonferroni")
colnames(localeQTLForASEWithPower)[31:34] <- c("NooriBYCount", "NooriRMCount", "AlbertBYCount", "AlbertRMCount")

nrow(localeQTLForASEWithPower[localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 & localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8,])
# 35 for >= 0.8, 34 for
# how many are NOT detected?
nrow(localeQTLForASEWithPower[localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 & localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8 & localeQTLForASEWithPower$sigDeciderBonferroni == 0 ,])
# 5:
#"YJR060W": CBF1
#"YKL210W": UBA1
#"YKR059W": TIF1
#"YLR261C" VPS63
#"YOR273C" TPO4

# how many found in only one dataset?
nrow(localeQTLForASEWithPower[localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 & localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8 & localeQTLForASEWithPower$sigDeciderBonferroni == 1 ,])
# 4
# how many found in two dataset?
nrow(localeQTLForASEWithPower[localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 & localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8 & localeQTLForASEWithPower$sigDeciderBonferroni == 2 ,])
# 26

# make plots where ASE power is highlighted:

pdf("ASE_vs_eQTLFC_withPower.pdf", width=12, height=12)
plotCol <- rep("#00000022", nrow(localeQTLForASEWithPower))
plotCol[localeQTLForASEWithPower$sigDeciderBonferroni == 1] <- "#FF000022"
plotCol[localeQTLForASEWithPower$sigDeciderBonferroni == 2] <- "red"
names(plotCol) <- localeQTLForASEWithPower$gene

circleCol <- rep("#00000000", nrow(localeQTLForASEWithPower))
circleCol[localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8] <- "lightblue"
circleCol[localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 & localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8] <- "blue"
names(circleCol) <- localeQTLForASEWithPower$gene

circleLwd <- rep(1, nrow(localeQTLForASEWithPower))
circleLwd[localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 & localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8] <- 2
names(circleLwd) <- localeQTLForASEWithPower$gene


par(mfrow=c(2,2))

plot(localeQTLForASEWithPower$ASE, localeQTLForASEWithPower$"naturalLMCoeff" * 2, pch=19, col=plotCol, ylab="local eQTL fold change", xlab="log2(ASE signal)", main="significant eQTL")
points(localeQTLForASEWithPower$ASE, localeQTLForASEWithPower$"naturalLMCoeff" * 2, col=circleCol, lwd=circleLwd)
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")

legend("bottomright", box.lty=0, col=c("#00000022", "#FF000022", "red","lightblue", "blue"), pch=c(19, 19, 19, 1, 1), legend=c("no significant ASE", "one significant ASE dataset", "two significant ASE datasets", "ASE power > 0.8 in one dataset", "ASE power > 0.8 in two datasets"))

textDecider <- ((localeQTLForASEWithPower$"sigDeciderBonferroni" == 0) & (localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 & localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8)) | (sign(localeQTLForASEWithPower$"naturalLMCoeff") != sign(localeQTLForASEWithPower$"ASE") & localeQTLForASEWithPower$"sigDeciderBonferroni" == 2 & sign(localeQTLForASEWithPower$"NooriFC") == sign(localeQTLForASEWithPower$"AlbertFC"))

text(localeQTLForASEWithPower$ASE[textDecider], localeQTLForASEWithPower$"naturalLMCoeff"[textDecider] * 2, labels = allNames[localeQTLForASEWithPower$gene[textDecider]], pos=4)


plot(localeQTLForASEWithPower$ASE, localeQTLForASEWithPower$"naturalLMCoeff" * 2, pch=19, col=plotCol, ylab="local eQTL fold change", xlab="log2(ASE signal)", main="significant eQTL", xlim=c(-2, 2), ylim=c(-2, 2))
points(localeQTLForASEWithPower$ASE, localeQTLForASEWithPower$"naturalLMCoeff" * 2, col=circleCol, lwd=circleLwd)
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")

plot(localeQTLForASEWithPower$ASE[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], localeQTLForASEWithPower$"naturalLMCoeff"[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8] * 2, pch=19, col=plotCol[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], ylab="local eQTL fold change", xlab="log2(ASE signal)", main="significant eQTL with significant ASE or high ASE power")
for(i in localeQTLForASEWithPower$gene[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8]){
        lines(-localeQTLForASEWithPower[i, c("NooriFC", "AlbertFC")], rep(localeQTLForASEWithPower[i,]$"naturalLMCoeff" * 2, 2), col=plotCol[i])
}
points(localeQTLForASEWithPower$ASE[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], localeQTLForASEWithPower$"naturalLMCoeff"[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8] * 2, col=circleCol[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], lwd=circleLwd[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8])
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")

plot(localeQTLForASEWithPower$ASE[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], localeQTLForASEWithPower$"naturalLMCoeff"[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8] * 2, pch=19, col=plotCol[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], ylab="local eQTL fold change", xlab="log2(ASE signal)", main="significant eQTL with significant ASE or high ASE power", xlim=c(-2, 2), ylim=c(-2, 2))
for(i in localeQTLForASEWithPower$gene[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8]){
    lines(-localeQTLForASEWithPower[i, c("NooriFC", "AlbertFC")], rep(localeQTLForASEWithPower[i,]$"naturalLMCoeff" * 2, 2), col=plotCol[i])
}
points(localeQTLForASEWithPower$ASE[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], localeQTLForASEWithPower$"naturalLMCoeff"[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8] * 2, col=circleCol[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8], lwd=circleLwd[localeQTLForASEWithPower$"sigDeciderBonferroni" > 0 | localeQTLForASEWithPower$NooriPowerBonferroni >= 0.8 | localeQTLForASEWithPower$AlbertPowerBonferroni >= 0.8])
abline(h=0, lwd=2, lty=2, col="grey")
abline(v=0, lwd=2, lty=2, col="grey")
abline(0, 1, lwd=2, lty=2, col="grey")


dev.off()

# discordant directions:
localeQTLForASEWithPower[(sign(localeQTLForASEWithPower$"naturalLMCoeff") != sign(localeQTLForASEWithPower$"ASE") & localeQTLForASEWithPower$"sigDeciderBonferroni" == 2 & sign(localeQTLForASEWithPower$"NooriFC") == sign(localeQTLForASEWithPower$"AlbertFC")),]
# 3 genes:
# YGR192C (TDH3)
# YMR089C (YTA12)
# YOR046C (DBP5)

# write a table for the paper:
datForLocalEQTLTable <- localeQTLForASEWithPower
datForLocalEQTLTable$naturalLMCoeff <- 2*datForLocalEQTLTable$naturalLMCoeff
colnames(datForLocalEQTLTable)[colnames(datForLocalEQTLTable) == "naturalLMCoeff"] <- "log2eQTLFoldChange"
colnames(datForLocalEQTLTable)[colnames(datForLocalEQTLTable) == "ASE"] <- "log2ASEFoldChange"

write.table(datForLocalEQTLTable, file="ASE_localeQTL_table_170323.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")





