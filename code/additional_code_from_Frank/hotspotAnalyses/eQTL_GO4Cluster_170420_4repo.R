args=(commandArgs(TRUE))
# we passed thisHotSpot directly from the shell
#eval(parse(text=args[[i]]))
# or, it might be safer (less cross-script messiness) to pass just the value and assign it in here

########################
#
# this code waas written to run on MSI at the University of Minnesota
# note hardcoded file names
#
########################


thisHotSpot = as.integer(args[1])

print(thisHotSpot)

library(topGO)
library(Rgraphviz)

setwd("RESULTFOLDER")
objectPath = "OBJECTFOLDER"
load(paste0(objectPath, "B.Forward.RData"))

peaksForHeatLasso = t(B)
colnames(peaksForHeatLasso) <- gsub(":", "_", colnames(peaksForHeatLasso))
colnames(peaksForHeatLasso) <- gsub("/", "_", colnames(peaksForHeatLasso))
#colnames(peaksForHeatLasso) <- sapply(colnames(peaksForHeatLasso), function(x){paste(strsplit(x, "_")[[1]][1], collapse="_")})
# careful, hardcoded:
#peaksForHeatBeta = t(ols.coefficient.matrix[15:125,])
#colnames(peaksForHeatBeta) <- sapply(colnames(peaksForHeatBeta), function(x){paste(strsplit(x, "_")[[1]][c(2:5)], collapse="_")})


myGene2GO = read.table("PATH_TO_ANNOTATION/gene_association_MOD_fromGO_160229.txt", stringsAsFactors = FALSE, sep="\t", header=FALSE)
# make that gene to GO list format they want:
gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
names(gene2GOList) = unique(myGene2GO[,1])

maxTopGenes = 100

getUpGenes <- function(thisScore, limitNumberReturned = TRUE, maxReturnNumber = maxTopGenes) {
    ret <- thisScore > 0
    if (limitNumberReturned & length(ret[ret]) > maxReturnNumber){
        ret <- thisScore >= sort(thisScore, decreasing=TRUE)[maxReturnNumber]
    }
    return(ret)
}
getDownGenes <- function(thisScore, limitNumberReturned = TRUE, maxReturnNumber = maxTopGenes) {
    ret <- thisScore < 0
    if (limitNumberReturned & length(ret[ret]) > maxReturnNumber){
        ret <- thisScore <= sort(thisScore, decreasing=FALSE)[maxReturnNumber]
    }
    return(ret)
}
getNonZeroGenes <- function(score, limitNumberReturned = TRUE, maxReturnNumber = maxTopGenes) {
    ret = getUpGenes(score, limitNumberReturned, maxReturnNumber) | getDownGenes(score, limitNumberReturned, maxReturnNumber)
    return(ret)
}

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

outGOFolder = paste0("topGO/topGenes", maxTopGenes, "/", sep="")
thisSigThres = 0.0005

geneSelFunctions = c("getUpGenes", "getDownGenes", "getNonZeroGenes")
scoreObjects = c(rep("peaksForHeatLasso", 3))
scoreDirections = c(1, 1, 1)
statisticsForGO = c("fisher", "fisher", "fisher")
outFileNames = c("FisherUp", "FisherDown", "FisherBoth")


for (GOcat in c("BP", "CC", "MF")){
#for (GOcat in c("BP")){
    
    hotspotGOResult = vector("list", ncol(peaksForHeatLasso))
    hotspotGOData = vector("list", ncol(peaksForHeatLasso))
    hotspotGOTable = vector("list", ncol(peaksForHeatLasso))
    names(hotspotGOResult) = colnames(peaksForHeatLasso)
    names(hotspotGOData) = colnames(peaksForHeatLasso)
    names(hotspotGOTable) = colnames(peaksForHeatLasso)
    
    #thisHotSpot gets read in from the shell script
    #for (thisHotSpot){
        
        hotspotGOResult[[thisHotSpot]] = vector("list", length(geneSelFunctions))
        hotspotGOData[[thisHotSpot]] = vector("list", length(geneSelFunctions))
        hotspotGOTable[[thisHotSpot]] = vector("list", length(geneSelFunctions))
        
        print(c(thisHotSpot, GOcat))
        
        for (i in 1:length(geneSelFunctions)){
            theseAllGenes = scoreDirections[i] * get(scoreObjects[i])[,thisHotSpot]
            # necesary because with the betas, sometimes a single gene is NA, which blows up the test
            theseAllGenes = theseAllGenes[!is.na(theseAllGenes)]
            GOData <- new("topGOdata", ontology=GOcat,  geneSel = get(geneSelFunctions[i]), allGenes = theseAllGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
            hotspotGOData[[thisHotSpot]][[i]] = GOData
            # write the "significant" genes
            write.table(t(t(sigGenes(GOData))), file=paste0(outGOFolder, GOcat, "/", "GOGraph_hotspot_", colnames(peaksForHeatLasso)[thisHotSpot], "_", outFileNames[i],"_sigGenes.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
            if(length(sigGenes(GOData)) > 0){
                GOresult <- runTest(GOData, algorithm="classic", statistic=statisticsForGO[i])
                #GOresult <- runTest(GOData, algorithm="weight01", statistic=statisticsForGO[i])
                hotspotGOResult[[thisHotSpot]][[i]] = GOresult
                hotspotGOTable[[thisHotSpot]][[i]] = GenTable(GOData, pValue = GOresult, topNodes = length(score(GOresult)), numChar=140)
                write.table(hotspotGOTable[[thisHotSpot]][[i]], file= paste0(outGOFolder, GOcat, "/", "GOGraph_hotspot_", colnames(peaksForHeatLasso)[thisHotSpot], "_", outFileNames[i], "_GOResultTable.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
            }
        }
        
        # now plot the results
        # and write the tables
        pdf(paste0(outGOFolder, GOcat, "/", "GOGraph_hotspot_", colnames(peaksForHeatLasso)[thisHotSpot], ".pdf"), width=15, height=15)
        for (i in 1:length(hotspotGOResult[[thisHotSpot]])){
            tryTester <- try({plotGOToTree(hotspotGOData[[thisHotSpot]][[i]], hotspotGOResult[[thisHotSpot]][[i]], sigThres = thisSigThres)})
            if(class(tryTester) == "try-error"){
                plot(1,1)
                legend("topright", legend="something went wrong when plotting, perhaps because there were no significant genes to be tested", cex=0.4, box.lty=0)
            }
        }
        dev.off()
}
