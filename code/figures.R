# eQTL per gene 
eQTL.per.gene.table=table(sapply(peaks.per.gene, nrow))
eQTL.per.gene.table=c( 5720-length(peaks.per.gene), eQTL.per.gene.table)
names(eQTL.per.gene.table)[1]='0'
pdf(file='/data/eQTL/FiguresForPublication/Figure1.pdf', width=8, height=5,useDingbats=F)
barplot(eQTL.per.gene.table, ylab ='Number of genes', xlab='Number of eQTL per gene')
dev.off()




pdf(file='/data/eQTL/FiguresForPublication/Figure2.pdf', width=8, height=8,useDingbats=F )
par(xaxs='i', yaxs='i')
plot(h2As,qCVs, xlim=c(0,1), ylim=c(0,1), col="#00000022", pch=19,
     xlab=expression(paste('Narrow-sense heritability', " (h"^"2", ")")), ylab='Phenotypic variance explained by QTL')
abline(0,1)
dev.off()




pdf(file='/data/eQTL/FiguresForPublication/Figure3.pdf', width=11, height=6,useDingbats=F )
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(dtrans, xlab='Fraction of phenotypic variance', ylab='Density', xlim=c(0,0.2),  type='n', main='')
polygon(dtrans, col=rgb(0, 0, 1,0.3), border=NA)
polygon(dcis, col=rgb(1, 0, 0,0.3), border=NA)
legend('topright', 
       legend=c('Distant eQTL effects', 'Local eQTL effects'), fill=c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)),
       box.lty=0)
sf= max(c(dtrans$y, dcis$y))
points((pvar), n1000*sf, type='l',col='black',lwd=3, lty=1, xlab=NA, ylab=NA)
#points((pvar), n1000.int*sf, type='l',col='red',lwd=3, lty=1,  xlab=NA, ylab=NA)
axis(side=1, at =c(.01,.02,.03,.04), labels=c(.01,.02,.03,.04))
axis(side=4, at=seq(0,1, .2)*sf, labels= seq(0,1, .2))
mtext(side = 4, line = 3, "Power")
dev.off()




pdf(file='/data/eQTL/FiguresForPublication/Figure5.pdf', width=12, height=12,useDingbats=F )
eQTL_bigPlot(all.peaks.OD, gcoord.key,marker.GR)
dev.off()


B.count.plot=(apply(B, 1, function(x) sum(x!=0)))
xBcp=marker.GR$gcoord[match(names(B.count.plot), marker.GR$mname)]
pdf(file='/data/eQTL/FiguresForPublication/Figure6.pdf', width=12, height=6,useDingbats=F )
par(yaxs='i')
plot(xBcp, B.count.plot, type='h', lwd=3, xaxt='n', xlab='Genome position of hotspots (chromosome name)', ylab='Number of genes' ,bty='n')
abline(v=gcoord.key, col='grey' , lty=2)
midpoints=(c(gcoord.key[-1])+gcoord.key[-17])/2
axis(1, at=midpoints, labels=as.roman(1:16), tick=F)
dev.off()




pdf(file='/data/eQTL/FiguresForPublication/Figure8.pdf', width=6, height=6,useDingbats=F )
B.count.plot=(apply(B, 1, function(x) sum(x!=0)))
#cor.test(B.count.plot, sapply(Blist, function(x) median(abs(x))))
#cor.test(B.count.plot, sapply(Blist, function(x) mean(abs(x))))
#plot(B.count.plot, sapply(Blist, function(x) median(abs(x))))
B.count.plot.max=apply(B, 1, function(x) {
                        y=x[x!=0];
                        #max(y)
                        median(sort(abs(y), decreasing=T)[1:50])
                        }
                       )
cor.test(B.count.plot, B.count.plot.max, method='spearman')
plot(B.count.plot, B.count.plot.max, xlab='Number of transcripts linking to a hotspot',
     ylab='Median effect size of top 50 linkages to a hotspot (in sd units) ')
dev.off()


pdf(file='/data/eQTL/FiguresForPublication/Figure9.pdf', width=6, height=6,useDingbats=F )
B.count.ploth=(apply(B, 2, function(x) sum(x!=0)))




pdf(file='/data/eQTL/FiguresForPublication/SFigure1.pdf', width=6, height=6,useDingbats=F )
plot(reads, c(mean(h2A.OD), apply(h2A.downsample,2, mean)), ylim=c(0,.5), xlab='reads per individual', ylab=expression(h^2), xlim=c(0,6e6), type='p' ,cex=1)
yin=as.numeric(c(mean(h2A.OD), apply(h2A.downsample,2, mean)))
xin=reads
nfit=nls(yin~d*xin/(b+xin))
x2=seq(1e4, 20e6,length.out=100)
predicty=((coef(nfit)['d']*x2))/(coef(nfit)['b']+x2)
points(x2, predicty, col='black', type='l')
dev.off()


pdf(file='/data/eQTL/FiguresForPublication/SFigure2.pdf', width=6, height=6,useDingbats=F )
plot(colMeans(t.tpm.matrix), h2A.OD, xlab='Average tpm per gene', ylab=expression(h^2), col='#00000044' , ylim=c(0,1), pch=20)    
dev.off()

pdf(file='/data/eQTL/FiguresForPublication/SFigure3.pdf', width=6, height=6,useDingbats=F )
ppg.count =sapply(peaks.per.gene, nrow)
dppg=cbind( ppg.count, h2A.OD[names(ppg.count)])
dppg=rbind(dppg, cbind(0, h2A.OD[names(h2A.OD)[!(names(h2A.OD) %in% names(ppg.count))]]))
plot(dppg[,1], dppg[,2],  xlab= 'number of eQTL per gene', ylab=expression(h^2), ylim=c(0,1), col='#00000022' ,pch=20 )
dev.off()

pdf(file='/data/eQTL/FiguresForPublication/SFigure4.pdf', width=12, height=6,useDingbats=F )
par(mfrow=c(1,2))
max.exp.q=sapply(peakModel, function(x) max(x$var.exp))
plot(max.exp.q, h2A.OD[names(max.exp.q)],  xlab= 'variance explained by strongest eQTL', ylab=expression(h^2), ylim=c(0,1), col='#00000022' ,pch=20)
abline(0,1, lty=2, col='grey')

max.exp.q=sapply(peakModel, function(x) max(x$var.exp))
plot(max.exp.q/(h2A.OD[names(max.exp.q)]) , h2A.OD[names(max.exp.q)],  xlab= 'fraction of additive heritability explained by strongest eQTL', ylab=expression(h^2), xlim=c(0,1.2), ylim=c(0,1), col='#00000022' ,pch=20)
abline(v=1, lty=2, col='grey')
#abline(0,1, lty=2, col='grey')
dev.off()

# cis no trans histogram
pdf(file='/data/eQTL/FiguresForPublication/histogram_figure.pdf', width=6, height=6,useDingbats=F )

atcto=data.frame(do.call('rbind', tctos))
stripchart(atcto, vertical=T, pch=19, cex=.5, ylim=c(0,3000), method='jitter', 
           jitter=.25, at=c(0:6), group.names=c(0:6), xlim=c(-.25,6.25), col='#00000022',
           xlab='(trans eQTL peaks not in hotspots) per local eQTL bin', ylab='number of local eQTL bins'
           )
tcto=table(cisTransOver)
tcto[7]=sum(tcto[7:length(tcto)])
tcto=tcto[1:7]
points(c(0:6), tcto, pch='-', cex=3, col='red')
dev.off()


vc_mc_int$AA/vc_mc_int$A
vcr=vcmc.norm[,2]/vcmc.norm[,1]
#mean(vcmc.norm[,2])/mean(vcmc.norm[,1])

vcr=vcr[-5322]
vcr=vcr[-1591]

vcmc.norm=vc_mc_int/rowSums(vc_mc_int)
# interaction variance 
vc_mc_int=data.frame(t(sapply(vc_multiple_components, function(x) x[['A_AA']]$sigma)))

pdf(file='/data/eQTL/FiguresForPublication/interaction_variance_figure.pdf', width=6, height=8,useDingbats=F )

stripchart(vc_mc_int[,-3], vertical=T, method='jitter', pch=19, col='#00000010', ylim=c(0,1), 
           ylab='fraction of phenotypic variance', group.names=c('A','AA'), at=c(0,.5), xlim=c(-.1, .75) )
boxplot(t(sapply(vc_multiple_components, function(x) x[['A_AA']]$sigma[c(1,2)])), add=T, border='darkred',
        outline=FALSE, boxwex=.10, xaxt='n', at=c(0.175,.675))
dev.off()


# 2D connected points 
pdf(file='/data/eQTL/FiguresForPublication/Figure_connected_interactions.pdf', width=8, height=8,useDingbats=F )

lm.2D.df.long=do.call('rbind', lm.2D.df)
cint=cbind(lm.2D.df.long$m2.gcoord,lm.2D.df.long$m1.gcoord)
cint=t(apply(cint,1, function(x) if(x[1]<x[2]) {return(x)} else { return(rev(x)) }))
xlim.ind=1:length(marker.GR)
par(xaxs='i', yaxs='i')
plot(cint[,1], lm.2D.df.long$t.gcoord,    pch=20, xlab='QTL positions', ylab='transcript position',
      cex=lm.2D.df.long$vexp*1,
     xaxt='n', yaxt='n',xlim=range( marker.GR$gcoord[xlim.ind]), ylim=range( marker.GR$gcoord[xlim.ind]),
     type='n'
     )
#points(cint[,2], lm.2D.df.long$t.gcoord,   pch=20,cex=lm.2D.df.long$vexp*1)
segments(cint[,1], lm.2D.df.long$t.gcoord, cint[,2], lm.2D.df.long$t.gcoord, lwd=lm.2D.df.long$vexp*25, col='#00000088')
         #col=ifelse(lm.2D.df.long$coef2D>0, 'orange', 'purple'))
         #ifelse(lm.2D.df.long$partner.cis, '#00000055', '#ff000044'))
abline(v=gcoord.key, lty=2, col='lightblue')
abline(h=gcoord.key, lty=2, col='lightblue')
gk=(gcoord.key+gcoord.key[-1])/2
gk=gk[1:16]
axis(1, at=gk, labels=as.roman(1:16), tick=F)
axis(2, at=gk, labels=as.roman(1:16), tick=F)
dev.off()


# heritability tables
heritability_A=vcA.OD
rownames(heritability_A)=names(h2A.OD)
colnames(heritability_A)=c('A', 'E')

heritability_A_AA=(t(sapply(vc_multiple_components, function(x) x[['A_AA']]$sigma)))
colnames(heritability_A_AA)[3]='E'

heritability_tables=list(heritability_A=data.frame(heritability_A), 
                         heritability_A_AA=data.frame(heritability_A_AA))
WriteXLS(heritability_tables,'/data/eQTL/RData/Supplementary_Data_4.xlsx',row.names=T) 


lm.2D.df.long
marginal2d.table
full2d.table

marginal2d.scan=data.frame(marginal2d.table)
marginal2d.scan=marginal2d.scan[,c(1,4,6,7)]

full2d.scan=data.frame(full2d.table)
full2d.scan=full2d.scan[,c(1,4,6,7)]

between_additive_loci.scan=data.frame(lm.2D.df.long[,c(9,1,2,3,4,8)])
colnames(between_additive_loci.scan)[1]='trait'

twoD_tables=list(between_additive_loci.scan=between_additive_loci.scan,
                 marginal2d.scan=marginal2d.scan,
                 full2d.scan=full2d.scan)
WriteXLS(twoD_tables, '/data/eQTL/RData/Supplementary_Data_15.xlsx')
