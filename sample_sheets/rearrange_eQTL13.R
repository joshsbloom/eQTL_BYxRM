sample.sheet=read.delim('/data/eQTL/sample_sheets/sample_sheet_eQTL13.csv', sep=',', header=T, stringsAsFactors=F)
rotate = function(mat) t(mat[nrow(mat):1,,drop=FALSE])
nmat=t(matrix(paste0(rep(toupper(letters[1:8]), each=12), seq(1:12)), 12,8))
rnmat=rotate(rotate(nmat))

#match(sample.annot.df$growth.well, rnmat)
on=cbind(as.vector(nmat), as.vector(rnmat))
mm=match(sample.sheet$Sample_Well, on[,1])
m2=on[mm,2]
m3=match(sample.sheet$Sample_Well, m2)


ndf=data.frame(sample.sheet[96:1, c(1,2,3)], sample.sheet[,c(4,5,6,7,8)])
ndf$Sample_ID=seq(1:96)
write.table(ndf, sep=',', file='/data/eQTL/sample_sheets/sample_sheet_eQTL13_flipped.csv', quote=F, row.names=F)
