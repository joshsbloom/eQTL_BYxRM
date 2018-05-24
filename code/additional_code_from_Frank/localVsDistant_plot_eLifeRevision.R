library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

serverPrefix = ""

load(paste0(serverPrefix, "R_allPeaksODPad_161213.RData", sep=""))
eQTLGRangesPad <- GRanges(seqnames = allPeaksODPad[,"chr"], ranges = IRanges(start=sapply(allPeaksODPad[,"CI.l"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])}), end=sapply(allPeaksODPad[,"CI.r"], function(x){as.integer(strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])})))
eQTLGRangesPad$gene <- allPeaksODPad$gene
eQTLGRangesPad$cis <- allPeaksODPad$cis

# all genes
sumAfterThis <- 5
varExpSummed <- t(sapply(unique(allPeaksODPad$gene), function(x){
    #    print(x)
    ret <- rep(0, sumAfterThis+3)
    transPeakVars <- sort(allPeaksODPad$"var.exp"[(!allPeaksODPad$cis) & allPeaksODPad$gene == x], decreasing=TRUE)
    ret[1] <- sum(allPeaksODPad$"var.exp"[allPeaksODPad$cis & allPeaksODPad$gene == x])
    ret[2] <- sum(transPeakVars)
    if(length(transPeakVars) > 0){
    for (i in 1:(min(c(sumAfterThis, length(transPeakVars))))){
        ret[i+2] <- transPeakVars[i]
    }
    if(length(transPeakVars) > sumAfterThis){
        ret[sumAfterThis+3] <- sum(transPeakVars[(sumAfterThis+1):length(transPeakVars)])
    }
    }
    ret
}))
colnames(varExpSummed) <- c("local", "distant", "distant1", "distant2", "distant3", "distant4", "distant5", "distant6Higher")
#colnames(varExpSummed) <- c("local", "distant", "distant1", "distant2", "distant3", "distant4Higher")

varExpSummedWithBoth <- varExpSummed[varExpSummed[,"local"] > 0 & varExpSummed[,"distant"] > 0,]

save(varExpSummed, file="varExpSummed.RData")

geneSets <- list(varExpSummedWithBoth, varExpSummed)
mainTitles <- c("at least one local and one distant", "all genes")
pdfNames <- c("figure1D_eLife_oneLocalOneDistant.pdf", "figure1D_eLife_allGenes.pdf")

for (i in 1:2){
# histogram-style plots
    pdf(pdfNames[i], width=7, height=4)
    plotCols=c("red", "blue")
    # compute the fraction distant/total
    #fractionDistantWithBoth <- varExpSummedWithBoth[,"distant"] / rowSums(varExpSummedWithBoth)
    fractionDistantWithBoth <- geneSets[[i]][,"distant"] / rowSums(geneSets[[i]][,1:2])
    fiftyFiftyPoint <- which(abs(sort(fractionDistantWithBoth) - 0.5) == min(abs(sort(fractionDistantWithBoth) - 0.5))) / length(fractionDistantWithBoth)

    bP <- barplot(t(geneSets[[i]][,1:2][order(fractionDistantWithBoth, geneSets[[i]][,"distant"]),]), beside=FALSE, border=NA, col=plotCols, xlab="Genes", ylab="Summed eQTL variance", xaxt="n", main="Sorted by distant fraction", sub=mainTitles[i])
    legend("topright", fill=plotCols, legend=c("local", "distant"), box.lty=0)
    abline(v=bP[fiftyFiftyPoint * length(fractionDistantWithBoth)], lty=2, col="black")
    axis(1, at = bP[fiftyFiftyPoint * length(fractionDistantWithBoth)], labels=round(fiftyFiftyPoint, 2))


    bP <- barplot(t((geneSets[[i]][,1:2]/rowSums(geneSets[[i]][,1:2]))[order(fractionDistantWithBoth),]), beside=FALSE, border=NA, col=plotCols, xlab="Genes", ylab="Fraction of eQTL variance", xaxt="n", main="Sorted by distant fraction", sub=mainTitles[i])
    legend("topright", fill=plotCols, legend=c("local", "distant"), box.lty=0)
    abline(h=0.5, col="white")
    abline(v=bP[fiftyFiftyPoint * length(fractionDistantWithBoth)], lty=2, col="white")
    axis(1, at = bP[fiftyFiftyPoint * length(fractionDistantWithBoth)], labels=round(fiftyFiftyPoint, 2))


    barplot(t(geneSets[[i]][,1:2][order(rowSums(geneSets[[i]][,1:2])),]), beside=FALSE, border=NA, col=plotCols, xlab="Genes", ylab="Summed eQTL variance", xaxt="n", main="Sorted by total", sub=mainTitles[i])
    legend("topleft", fill=plotCols, legend=c("local", "distant"), box.lty=0)

    plot((1:length(fractionDistantWithBoth) / length(fractionDistantWithBoth)), sort(fractionDistantWithBoth), type="l", ylab="Fraction of eQTL variance from trans eQTL", xlab="Fraction of genes", sub=mainTitles[1])
    abline(h=0.5, col="red")
    abline(v=fiftyFiftyPoint, lty=2)
    axis(1, at = fiftyFiftyPoint, labels=round(fiftyFiftyPoint, 2))

    dev.off()
}

# with the trans peaks split
#pdfNames <- c("figure1D_eLife_oneLocalOneDistant_withDistantSplit.pdf", "figure1D_eLife_allGenes_withDistantSplit.pdf")
#pdfNames <- c("figure1D_eLife_oneLocalOneDistant_withDistantSplit_invertedTransSort.pdf", "figure1D_eLife_allGenes_withDistantSplit_invertedTransSort.pdf")
pdfNames <- c("figure1D_eLife_oneLocalOneDistant_withDistantSplit_colorTest.pdf", "figure1D_eLife_allGenes_withDistantSplit_colorTest.pdf")
for (i in 1:2){
    # histogram-style plots
    pdf(pdfNames[i], width=7, height=4)
    plotCols=c("red", rev(brewer.pal(sumAfterThis + 1, "Blues")))
    
    plotCols1 <- plotCols
    plotCols1[4:7] <- plotCols1[5]
    # compute the fraction distant/total
    #fractionDistantWithBoth <- varExpSummedWithBoth[,"distant"] / rowSums(varExpSummedWithBoth)
    fractionDistant <- geneSets[[i]][,"distant"] / rowSums(geneSets[[i]][,1:2])
    fiftyFiftyPoint <- which(abs(sort(fractionDistant) - 0.5) == min(abs(sort(fractionDistant) - 0.5))) / length(fractionDistant)
    fractionOfDistantInFirst <- geneSets[[i]][,"distant1"] / rowSums(geneSets[[i]][,3:8])
    
    bP <- barplot(t(geneSets[[i]][,c(1, 3:8)][order(-geneSets[[i]][,1], -geneSets[[i]][,3], geneSets[[i]][,4], geneSets[[i]][,5], geneSets[[i]][,6]),]), beside=FALSE, border=NA, col=plotCols1, xlab="Genes", ylab="Summed eQTL variance", xaxt="n", main="Sorted by variance per local and distant eQTL", sub=mainTitles[i], ylim=c(0, 1))
    legend("topright", fill=plotCols[1:4], legend=c("local", "1st distant", "2nd distant", "sum of distant 3rd and higher"), box.lty=0)
    #abline(v=bP[fiftyFiftyPoint * length(fractionDistant)], lty=2, col="black")
    #axis(1, at = bP[fiftyFiftyPoint * length(fractionDistant)], labels=round(fiftyFiftyPoint, 2))
    axis(1, at = bP[pretty(1:length(fractionDistant))+1], pretty(1:length(fractionDistant)))
    
    # this plot looks better with many colors
    plotCols=c("red", rev(brewer.pal(sumAfterThis + 1, "Blues")))
    plotCols[5:7] <- plotCols[6]
    bP <- barplot(t((geneSets[[i]][,c(1, 3:8)]/rowSums(geneSets[[i]][,c(1, 2)]))[order(fractionDistant, 1/fractionOfDistantInFirst),]), beside=FALSE, border=NA, col=plotCols, xlab="Genes", ylab="Fraction of eQTL variance", xaxt="n", main="Sorted by fraction of variance in distant eQTL", sub=mainTitles[i])
    legend("topright", fill=plotCols[1:5], legend=c("local", "1st distant", "2nd distant", "3rd distant", "sum of distant 4th and higher"), box.lty=0)
    abline(h=0.5, col="white")
    abline(v=bP[fiftyFiftyPoint * length(fractionDistant)], lty=2, col="white")
    axis(1, at = bP[fiftyFiftyPoint * length(fractionDistant)], labels=round(fiftyFiftyPoint, 2))
    axis(1, at = bP[pretty(1:length(fractionDistant))+1], pretty(1:length(fractionDistant)))

    barplot(t(geneSets[[i]][,c(1, 3:8)][order(-rowSums(geneSets[[i]][,c(1, 3:8)])),]), beside=FALSE, border=NA, col=plotCols, xlab="Genes", ylab="Summed eQTL variance", xaxt="n", main="Sorted by total eQTL variance", sub=mainTitles[i], ylim=c(0, 1))
    legend("topright", fill=plotCols[1:5], legend=c("local", "1st distant", "2nd distant", "3rd distant", "sum of distant 4th and higher"), box.lty=0)
    axis(1, at = bP[pretty(1:length(fractionDistant))+1], pretty(1:length(fractionDistant)))

    plot((1:length(fractionDistant) / length(fractionDistant)), sort(fractionDistant), type="l", ylab="Fraction of eQTL variance from trans eQTL", xlab="Fraction of genes", sub=mainTitles[1])
    abline(h=0.5, col="red")
    abline(v=fiftyFiftyPoint, lty=2)
    axis(1, at = fiftyFiftyPoint, labels=round(fiftyFiftyPoint, 2))

    dev.off()
}






