library(GenomicRanges)

serverPrefix = ""

allVariants <- read.table(paste0(serverPrefix, "gdata_42k_VEP_wExtraMarkersAttached.txt", sep=""), sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
allVariants <- cbind(allVariants, sapply(allVariants[,1], function(x){strsplit(x, ":")[[1]][1]}))
colnames(allVariants)[ncol(allVariants)] <- "chr"
allVariants <- cbind(allVariants, as.integer(sapply(allVariants[,1], function(x){strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1]})))
colnames(allVariants)[ncol(allVariants)] <- "pos"


# TFs:
TFs <- read.table(paste0(serverPrefix, "transcription_factor_activity_sequencespecific_DNA_binding_annotations_withChildTerms_170403_MOD.txt", sep=""), stringsAsFactors=FALSE, head=TRUE, sep="\t")

length(unique(TFs$Gene))
# 209

# how many variants, missense, etc per TF?
variantsInTFs <- t(sapply(unique(TFs$Gene.Systematic.Name), function(x){
    ret=c(0, 0, 0)
    ret[1] <- length(which(allVariants[,9] == x))
    ret[2] <- length(which(allVariants[,9] == x & allVariants[,5] == "HIGH"))
    ret[3] <- length(which(allVariants[,9] == x & allVariants[,5] == "MODERATE"))
    ret
}))
colnames(variantsInTFs) <- c("all", "HIGH", "MODERATE")

length(which(variantsInTFs[,"all"] > 0))
# 165
length(which(variantsInTFs[,"all"] == 0))
# 44


length(which(variantsInTFs[,"HIGH"] > 0))
# 8
allVariants[allVariants[,9] %in% rownames(variantsInTFs[which(variantsInTFs[,"HIGH"] > 0),]) & allVariants[,5] == "HIGH", c(1, 4, 6, 13)]

#4484                                            chrII:541484_CA/C
#4485          chrII:541486_TCATCATCATCATCATCATCATCATCATCATCATCA/T
#6146  chrIII:148614_T/TTGTTGGAATAAAAATCAACTATCATCTACTAACTAGTATTTA
#14010                                     chrVI:96024_CAT/CATATAT
#14012                                         chrVI:96038_G/GATAT
#29421                                           chrXI:411481_T/TA
#32209                                          chrXII:510208_GT/G
#36854                                         chrXIII:312171_G/GA
#43294                                            chrXV:389923_G/A
#45063                                           chrXV:656377_CG/C
#V4   V6               V13
#4484                          frameshift_variant TBS1               D/X <- this IS a hotspot at chrII:541139_G/A; "unknown function, so no TF motifs, no targets in SGD)
#4485                          frameshift_variant TBS1  NDDDDDDDDDDDD/KX
#6146      "stop_gained,protein_altering_variant" SRD1 E/VNTS**MIVDFYSNK <- close to chrIII:143912; no annotated targets; "involved in the processing of pre-rRNA to mature rRNA; contains a C2/C2 zinc finger motif"
#14010                         frameshift_variant GAT1           AY/AYIX <- clear hotspot chrVI:95849_T/C; TF enrichment is for GLN3, which binds the same consensus sequence according to SGD
#14012                         frameshift_variant GAT1             D/DIX
#29421 "frameshift_variant,stop_retained_variant" PUT3                 * <- clearly a hotspot chrXI:411756_T/C; no matching TFBS enrichment, no obvious reason why
#32209                         frameshift_variant RFX1               K/X <- clearly chrXII:510559_C/T; no TFBS enrichment but major targets are RNR genes, and these have huge effects here; interesting regulation OF this gene too
#36854                         frameshift_variant STB4              E/EX <- somewhat far away-ish from chrXIII:333449_G/A
#43294                                stop_gained HMS1               Q/* <- in CI for chrXV:377605_C/T; no matching TFBS enrichment, no good reason why
#45063                         frameshift_variant YRM1               A/X <- not in a hotspot region
                                                                      

length(which(variantsInTFs[,"MODERATE"] > 0))
# 141

