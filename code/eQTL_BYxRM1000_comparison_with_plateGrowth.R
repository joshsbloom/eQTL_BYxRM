library(topGO)
library(org.Sc.sgd.db)
library(fdrtool)
library(sisus)

# setup for GO enrichment analysis ------------------------------------------------------------------------
# slightly modified by JB, original code from Frank Albert
myGene2GO.full.table=read.delim('/data/eQTL/reference/go_annotations_sgd.txt',  sep='\t', header=F, stringsAsFactors=F)

# data frame of gene then GO
myGene2GO=cbind( sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,5])
myGene2GO=na.omit(myGene2GO)

SYS2ORF=unique(cbind(sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,3]))
SYS2ORF.key=list()
SYS2ORF.key[SYS2ORF[,1]]=SYS2ORF[,2]
gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
names(gene2GOList) = unique(myGene2GO[,1])


plotGOToTree <- function(GOdat, GOres, sigThres = 0.0005){
     # only plot if there are any significant GO terms (SEE ABOVE for 
    #"significance"; I am somewhat lenient here):
    #     # we need these extra lines because very small p-values are 
    #reported as a text string "< X", rather than a numeric
     toTest <- as.numeric(GenTable(GOdat, pVal = GOres)[1,6])
     if(is.na(toTest)){toTest <- 0.000000000000000000000000000001}
     if (toTest < sigThres){
         showSigOfNodes(GOdat, score(GOres), firstSigNodes = 
        length(score(GOres)[score(GOres) < sigThres]), useInfo = "all")
     }else{
         plot(1,1)
         legend("topright", legend="no significant GO categories", 
        cex=0.4, box.lty=0)
     }
}
#--------------------------------------------------------------------------------------------------------------


# Load latest phenotype data 
load('/data/rr/Phenotyping/NORMpheno.RData')
BYxRM_plate.phenos=lapply(NORMpheno, function(x) {x[grep('^A', names(x))] })
plate.pheno.avg=sapply(BYxRM_plate.phenos, function(x) sapply(x, mean, na.rm=T) )
# load eQTL expression dataa 
load('/data/eQTL/RData/log2_t.tpm.matrix.RData')
eQTL.seg.names=sapply(strsplit(rownames(t.tpm.matrix), '-'), function(x) x[1]) 

plate.oset=which((rownames(plate.pheno.avg) %in% eQTL.seg.names))
eQTL.oset=which(  eQTL.seg.names %in% rownames(plate.pheno.avg) )

ppa=plate.pheno.avg[plate.oset,]
eoa=t.tpm.matrix[eQTL.oset,]
# load genotype data 
load('/data/eQTL/RData/gdata.RData')

# load marker and transcript annotation data
load('/data/eQTL/RData/markerAnnotation.RData')
load('/data/eQTL/RData/geneAnnotation.RData')

goa=gdata[eQTL.oset,]
colnames(ppa)
#
#R> (organisms[18,])
#           Name      orgPackage      keggPrefix 
#        "Yeast" "org.Sc.sgd.db"           "sce" 
#library(FGNet)
load('/data/eQTL/RData/covariates.OD.RData')

for(i in 1:41) { 
    print(i)
    cg=cor(ppa[,i],goa, use='pairwise.complete.obs')
    drexp= scale(residuals(lm(eoa~covariates.OD[eQTL.oset,])))
      
    c1=cor(ppa[,i], drexp, use='pairwise.complete.obs')
    #lasso expresion
    #    lasso.test=cv.EBglmnet(drexp, as.vector(ppa[,i]))
    # for fluconazole
    #                                    Estimate Std. Error t value Pr(>|t|)
    #(Intercept)                           0.0135     0.2479    0.05     0.96
    #drexp[, lasso.test$fit[, 1]]YAL056W   1.4779     0.2487    5.94  4.0e-09
    #drexp[, lasso.test$fit[, 1]]YAR028W  -1.9066     0.2482   -7.68  4.0e-14
    #drexp[, lasso.test$fit[, 1]]YHL008C  -1.6796     0.2492   -6.74  2.8e-11
    #drexp[, lasso.test$fit[, 1]]YMR202W   4.7299     0.2497   18.94  < 2e-16

    png(filename=paste0('/data/eQTL/RData/GO/', filename.clean(colnames(ppa)[i]), '.png'), width=1920, height=1080 )
    par(mfrow=c(2,1))
    plot(marker.GR$gcoord, cg[1,]^2, ylab='r^2', main=colnames(ppa)[i])
    plot(gene.GR$gcoord, c1[1,]^2, ylab='r^2')
    #identify(gene.GR$gcoord, c1[1,]^2, colnames(eoa))
    dev.off()

    fc1=fdrtool(c1[1,], statistic='correlation')
    if(sum(fc1$qval<.4)>10 ){
    geneList= names(which(fc1$qval<.40)) 
    # is important this is coded as numbers and not true or false !! why??!!!??, this is awful
    testGenes= factor(0+(colnames(t.tpm.matrix) %in% geneList)) # colnames(t.tpm.matrix)[(which(vcA[,1]>.7))]) #totest )
    names(testGenes)=colnames(t.tpm.matrix)

    testGeneSet=names(testGenes[testGenes==1])
    for( thisOntology in c("BP", "MF", "CC") ) {
        GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
        GOresult = runTest(GOData, algorithm="classic", statistic="fisher")
        pdf(file=paste0('/data/eQTL/RData/GO/', filename.clean(colnames(ppa)[i]),'_', thisOntology, '.pdf'), width=25, height=25)
        plotGOToTree(GOData, GOresult, sigThres = 0.00005)
        dev.off()
        gt=GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult)))
        
        genesinterms=genesInTerm(GOData, gt[,1])
        genes.enriched.list=lapply(genesinterms, function(x) x[x%in%names(testGenes[testGenes==1])])
        genes.enriched.list.simple=lapply(genes.enriched.list, function(x) as.character(SYS2ORF.key[x]))
        gt$Genes=as.vector(sapply(genes.enriched.list.simple, paste, collapse=','))
        gt$GenesSystematic=   as.vector(sapply(genes.enriched.list, paste, collapse=','))
        write.table(gt, file=paste0('/data/eQTL/RData/GO/', filename.clean(colnames(ppa)[i]),'_', thisOntology, '.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
    }
    
    #file=paste0("topGO/hotspotGenes/hotspotGOresult_", thisOntology, 
    #"_bootPeakGenes.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, 
    #FGNet report
    #uallGenes=as.character(na.omit(unlist(as.list(org.Sc.sgdGENENAME)[colnames(c1)])))
    #geneList= names(which(fc1$qval<.50)) #which(c1[1,]^2>.025))
    #geneLabels = as.character(na.omit(unlist(as.list(org.Sc.sgdGENENAME)[geneList])))
    #annotations=c("GO_BP", "GO_MF", "GO_CC")
    #gout=fea_topGO(geneLabels, geneIdType="GENENAME", organism="Sc", genesUniverse=uallGenes, jobName=paste0('/data/eQTL/RData/GO/', filename.clean(colnames(ppa)[i])))
    #rm(gout)
    }
    #FGNet_report(gout)
}

# playing around with expression as a structured covariance term 
A.subset=tcrossprod(gdata[eQTL.oset,])/ncol(gdata)
rgeA=tcrossprod(rge)/ncol(rge)
AA = A.subset*A.subset
regress(ppa[,1]~1, ~A.subset, verbose=T)

rtvc=list()
for( i in 1:41) {
    # contrast 31 tunicamycin vs 12 fluconazole
    #rtvc[[colnames(ppa)[i]]]$r1=regress(ppa[,i]~1, ~A.subset, verbose=T, pos=c(T,T,T))
    #rtvc[[colnames(ppa)[i]]]$r2=regress(ppa[,i]~1, ~A.subset+rgeA, verbose=T, pos=c(T,T,T))
    #rtvc[[colnames(ppa)[i]]]$r3=regress(ppa[,i]~1, ~A.subset+AA+rgeA, verbose=T, pos=c(T,T,T,T))
    #rtvc[[colnames(ppa)[i]]]$r4=regress(ppa[,i]~1, ~A.subset+AA, verbose=T, pos=c(T,T,T,T))
}

#A + AA model
rta=sapply(rtvc, function(x) x$r4$sigma)
# A +AA + expression model
rtt=sapply(rtvc, function(x) x$r3$sigma)
vrtt= t(t(rtt)/rowSums(t(rtt)))
vrta= t(t(rta)/rowSums(t(rta)))

# some of the most interesting traits
# fluconazole 
# cobalt
#ph3
# formamide
#galactose
par(mfrow=c(2,1), oma=c(10,1,1,1))
barplot(vrtt[1:3,], las=2, ylim=c(0,1), main='A+AA+expression' )
#par(oma=c(10,1,1,1))
barplot(vrta[1:2,], las=2, ylim=c(0,1), main='A+AA')

