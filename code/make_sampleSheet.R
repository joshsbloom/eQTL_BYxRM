
#iterated 1-13 manually ----------------------------------------------------------------------------------------------------------------
rearranged.csv='/data/eQTL/rearrange/13.csv'
plate.name='BYxRM_eQTL_13'
sample.sheet='/data/eQTL/sample_sheets/sample_sheet_eQTL_13.csv'

x=read.delim(rearranged.csv, header=T, sep=',')

D501='TATAGCCT'
D502='ATAGAGGC'
D503='CCTATCCT'
D504='GGCTCTGA'
D505='AGGCGAAG'
D506='TAATCTTA'
D507='CAGGACGT'
D508='GTACTGAC'
Dcol=c(D501,       D502,       D503,      D504,       D505,       D506,       D507,       D508)
Dcoln=c('D501',       'D502',       'D503',       'D504',       'D505',       'D506',       'D507',       'D508')


 D701='ATTACTCG'
 D702='TCCGGAGA'
 D703='CGCTCATT'
 D704='GAGATTCC'
 D705='ATTCAGAA'
 D706='GAATTCGT'
 D707='CTGAAGCT'
 D708='TAATGCGC'
 D709='CGGCTATG'
 D710='TCCGCGAA'
 D711='TCTCGCGC'
 D712='AGCGATAG'
Drow=c(D701,     D702,     D703,     D704,     D705,     D706,     D707,     D708,     D709,     D710,     D711,     D712)
Drown=c('D701',     'D702',     'D703',     'D704',     'D705',     'D706',    'D707',     'D708',     'D709',     'D710',     'D711',     'D712')

#nums=do.call('rbind', strsplit(names(ym.by.plate[[1]]), '_'))[,2]
rows=rep(toupper(letters)[1:8], each=12)
cols=rep(1:12, 8)
newwell=t(matrix(paste(rows, cols, sep=''), 12,8))

idplate1=cbind(Dcoln,Dcoln,Dcoln,Dcoln,Dcoln,Dcoln,Dcoln,Dcoln,Dcoln,Dcoln,Dcoln,Dcoln)
idplate2=rbind(Drown,Drown,Drown,Drown,Drown,Drown,Drown,Drown)

iplate1=cbind(Dcol,Dcol,Dcol,Dcol,Dcol,Dcol,Dcol,Dcol,Dcol,Dcol,Dcol,Dcol)
iplate2=rbind(Drow,Drow,Drow,Drow,Drow,Drow,Drow,Drow)

x2=x[match(as.vector(newwell), as.character(x$new.well)),]

dfo=data.frame(
Sample_ID=seq(1,96),
Sample_Name=paste(x2$name, x2$old.plate, x2$old.well, sep="-"),
Sample_Plate=plate.name,
Sample_Well=as.vector(newwell),
I7_Index_ID=as.vector(idplate2),
index=as.vector(iplate2),
I5_Index_ID=as.vector(idplate1),
index2=as.vector(iplate1))

write.table(dfo, file=sample.sheet, row.names=F, sep=',', quote=F)


#---------------------------------------------------------------------------------------------------------------------



load('/media/kserver/kruglyak/raid1/home/jbloom/1000BYxRM/phenotyping/AllResultsProcessed.bin') 
yn=grep('YNB', names(Results.Processed))
RR=Results.Processed[yn]
RR=RR[grep('TRUE', names(RR))]
rr2=do.call('rbind', RR)
rr3=rr2[order(rownames(rr2)),]
srr3=split(rr3[,'g.effr.norm'], rownames(rr3))
sapply(srr3, mean, na.rm=T)->ym.n
srr3=split(rr3[,'g.effr'], rownames(rr3))
sapply(srr3, mean, na.rm=T)->ym
srr3=split(rr3[,'m.pxs'], rownames(rr3))
sapply(srr3, mean, na.rm=T)->ym.p

#yn=grep('YPD', names(Results.Processed))
#RR=Results.Processed[yn]
#RR=RR[grep('TRUE', names(RR))]
#rr2=do.call('rbind', RR)
#rr3=rr2[order(rownames(rr2)),]
#srr3=split(rr3[,'g.effr.norm'], rownames(rr3))
#sapply(srr3, mean, na.rm=T)->m.n
#srr3=split(rr3[,'g.effr'], rownames(rr3))
#sapply(srr3, mean, na.rm=T)->m

ym=ym.n

ym= ym[1:1056]
ym[is.na(ym)]=mean(ym, na.rm=T)

ym.by.plate=split(ym, do.call('rbind',strsplit(names(ym), '_'))[,1])
#convert 1-96 to well position
plate.name=do.call('rbind', strsplit(names(ym), '_'))[,1]

nums=do.call('rbind', strsplit(names(ym.by.plate[[1]]), '_'))[,2]
rows=rep(toupper(letters)[1:8], each=12)
cols=rep(1:12, 8)

nums=rep(nums, 11)
rows=rep(rows,11)
cols=rep(cols,11)

newwell=paste(rows, cols, sep='')

p.key=data.frame(name=names(ym), ynb.growth=ym, ynb.rank=rank(ym), old.plate=plate.name, old.well=newwell)

sp.key=p.key[order(p.key$ynb.growth),]


nums=do.call('rbind', strsplit(names(ym.by.plate[[1]]), '_'))[,2]
rows=rep(toupper(letters)[1:8], each=12)
cols=rep(1:12, 8)
newwell=paste(rows, cols, sep='')

spr.key = split(sp.key, rep(1:11, each=96))
for(i in 1:11) {
    spr.key[[i]]=data.frame(spr.key[[i]], new.plate=paste('BYxRM_eQTL_',sprintf("%02d",i), sep=''), new.well=newwell[sample(1:96)],vol=5) 
}

save(spr.key, file = '/data/eQTL/rearrange/rearrange_by_YNB.RData')

for(i in 1:11) {
    write.table(spr.key[[i]], file=paste('/data/eQTL/rearrange/', i, '.csv', sep=''),
                row.names=F, sep=',',quote=F)
}
