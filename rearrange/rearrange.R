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

ym->ynb.plate.growth
save(ynb.plate.growth, file='/data/eQTL/rearrange/ynb_plate_growth_1000BYxRM.RData')

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
