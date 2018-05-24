load("R_allPeaksODPad_161213.RData", sep="")

getGcoords = function ( chr , pos, spacing=0, sgd.table="sacCer3ChromLenghts.txt" ) {
    offind = as.vector(cumsum(read.table(sgd.table, header=FALSE, sep="\t")[,2] + spacing))
    offind=    offind[-length(offind)]
    offind= c(0, offind)
    names(offind) = as.character(read.table(sgd.table, header=FALSE, sep="\t")[,1])
    chr.off=as.numeric(sapply(chr, function(x) {offind[[x]]}))
    return(chr.off+pos)
}
sepBetweenChr <- 1e5
chromosomeDividers <- c(0, 230218, 1043402, 1360022, 2891955, 3468829, 3738990, 4829930, 5392573, 5832461, 6578212, 7245028, 8323205, 9247636, 10031969, 11123260, 12071326)

pdf("2D_alleQTL.pdf", width=8, height=8)
nCols=6

# scale by size
cutPal <- rep("#00000033", nCols)
cutBreaks <- seq(0, 1, length.out=nCols)
cutVectorCex <- cut(abs(allPeaksODPad$"var.exp"), breaks=cutBreaks, labels=FALSE)
cutVectorCol <- cutVectorCex
scaleFac=2

plot(allPeaksODPad$"marker.gcoord", allPeaksODPad$"gene.gcoord", pch=19, ylab="Gene position (chromosome)", xlab="eQTL position (chromosome)", xaxt="n", yaxt="n", xaxs="i", yaxs="i", col=cutPal[cutVectorCol], cex=(cutBreaks+0.1)[cutVectorCex]*scaleFac)

for(i in 1:length(chromosomeDividers)){
    abline(v=chromosomeDividers[i], col="lightblue", lty=2)
    abline(h=chromosomeDividers[i], col="lightblue", lty=2)
}
axis(1, at=sapply(1:16, function(i){chromosomeDividers[i] + ((chromosomeDividers[i+1] - chromosomeDividers[i])/2)}), labels=as.roman(1:16), tick=FALSE)
axis(2, at=sapply(1:16, function(i){chromosomeDividers[i] + ((chromosomeDividers[i+1] - chromosomeDividers[i])/2)}), labels=as.roman(1:16), tick=FALSE)

legend("topright", box.lty=1, title="variance fraction", legend=paste(cutBreaks[1:(nCols-1)], cutBreaks[2:nCols], sep=" - "), pch=19, col="black", pt.cex=(cutBreaks+0.1)*scaleFac, bg="white")

dev.off()



