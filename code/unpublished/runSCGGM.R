#
load('/data/eQTL/RData/hotspots_031016.RData')

hotspot.vector=sort(match(as.character(unlist( sapply(hotspots, function(x) sapply(x, function(y) as.character(y$peak))))), colnames(gdata)))

g.h=gdata.scaled[,hotspot.vector]
#local QTL mediation analysis
1:6288
t.tpm.matrix
i=1
p=scale(residuals.tpm.gbatch.OD)/

p.batch.cis=p[,1:6288]
for(i in 1:6288) {
    print(i)
    p.batch.cis[,i]=scale(residuals(lm(t.tpm.matrix[,i]~gbatch.fact+gdata[,closest.marker.to.transcript[i]])))
}

p.batch.cis.OD=p[,1:6288]
for(i in 1:6288) {
    print(i)
    p.batch.cis.OD[,i]=scale(residuals(lm(t.tpm.matrix[,i]~gbatch.fact+gdata[,closest.marker.to.transcript[i]]+OD.cov)))
}

vcA.bc=calcA(p.batch.cis, A)
h2.bc=(vcA.bc[,1]/(vcA.bc[,1]+vcA.bc[,2]))

#Ae=tcrossprod(p.batch.cis)/ncol(p.batch.cis)

library(pls)
p1=pcr(p.batch.cis[,1:6288]~p.batch.cis-1,ncomp=6)
pbc_mpc=p1$residuals[,,6]
pbc_mpc.s=scale(pbc_mpc)

vcA.p1=calcA(pbc_mpc.s, A)
h2.p1=(vcA.p1[,1]/(vcA.p1[,1]+vcA.p1[,2]))


vcA.p3=calcA(p.batch.cis.OD,A)
h2.p3=(vcA.p3[,1]/(vcA.p3[,1]+vcA.p3[,2]))

x11()
plot(h2.bc, h2.p3)

x11()
plot(h2.bc, h2.p1)





p2=pcr(p.batch.cis.OD[,1:6288]~p.batch.cis.OD-1,ncomp=6)
pbc_mpc2=p2$residuals[,,6]
pbc_mpc.s2=scale(pbc_mpc2)
vcA.p2=calcA(pbc_mpc.s2, A)
h2.p2=(vcA.p2[,1]/(vcA.p2[,1]+vcA.p2[,2]))



vcA.e=cbind(rep(NA, ncol(p.batch.cis)), rep(NA, ncol(p.batch.cis)))
rownames(vcA.e)=colnames(p.batch.cis)
eigA.e=doEigenA_forMM(pbc_mpc.s2,A)
    # calculate mixed model, one term for additive variance  -------------------------------------------
for(i in 1:ncol(p.batch.cis)){
       if(is.na(sd(p.batch.cis[,i]))) {
        next;
        }
        vcA.e[i,]=m.S(pbc_mpc.s2[,i], K=A,  theta=eigA.e$theta, Q=eigA.e$Q)
     print(i)
    }
h2.w=(vcA.e[,1]/(vcA.e[,1]+vcA.e[,2]))
plot(h2.w, h2.bc, xlab='h2 with first 6 PCs removed', ylab='h2')


peaks.per.gene
OD.mediation=list()
for(gs in names(peaks.per.gene) ){
        print(gs)
    #
    #gs=names(peaks.per.gene)[50]
    #gs=names(peaks.per.gene)[51]
    #gs=names(peaks.per.gene)[53]

     ppgl=peaks.per.gene[[gs]]
     gind=match(gs, colnames(t.tpm.matrix))
   
     #testing if growth (OD)  mediates QTL effects on expression levels  -----------------------
     sy=scale(t.tpm.matrix[,gind])
     sOD=scale(OD.cov)
    
     M1= lm(sy  ~ gbatch.fact+gdata.scaled[,ppgl$pmarker]-1)
     M2= lm(sOD ~ gbatch.fact+gdata.scaled[,ppgl$pmarker]-1)
     M3= lm(sy  ~ gbatch.fact+sOD+gdata.scaled[,ppgl$pmarker]-1)
     sM1= summary(M1)$coefficients
     sM2= summary(M2)$coefficients
     sM3= summary(M3)$coefficients
     
     #https://en.wikipedia.org/wiki/Sobel_test
     tau=sM1[which(grepl('gdata', rownames(sM1))),'Estimate']
     tau.se=sM1[which(grepl('gdata', rownames(sM1))),'Std. Error']

     a=sM2[which(grepl('gdata', rownames(sM2))),'Estimate']
     a.se=sM2[which(grepl('gdata', rownames(sM2))),'Std. Error']

     b=sM3[match('sOD', rownames(sM3)),'Estimate']
     b.se=sM3[match('sOD', rownames(sM3)),'Std. Error']    
    
     tau.prime=sM3[which(grepl('gdata', rownames(sM3))),'Estimate']
     tau.prime.se=sM3[which(grepl('gdata', rownames(sM3))),'Std. Error']    

     pooled.se = a^2 * b.se^2 + b^2*a.se^2 + a.se^2+b.se^2
     #Zscore=(tau-tau.prime)/pooled.se
     Zscore=(a*b)/pooled.se
     Zscore[abs(tau.prime)>abs(tau)]=NA
     p.med=pnorm(abs(Zscore), lower.tail=F)*2
     dff=data.frame(tau, tau.se, tau.prime, tau.prime.se, a,a.se, b,b.se,Zscore, p.med)
     rownames(dff)=ppgl$pmarker
     #----------------------------------------------------------------------------------------------
     OD.mediation[[gs]]=data.frame(ppgl, dff)
}
OD.med.pmarkers= do.call('c', sapply(OD.mediation, function(x) x$pmarker))
OD.med.pmed= do.call('c', sapply(OD.mediation, function(x) x$p.med))
OD.med.g= do.call('c', sapply(OD.mediation, function(x) x$gene))
OD.med.gind=match(OD.med.g, colnames(t.tpm.matrix))
OD.med.p=match(OD.med.pmarkers, colnames(gdata))

length(which(OD.med.pmed < .05/36455))
par(xaxs='i', yaxs='i')
plot(OD.med.p[which(OD.med.pmed < .05/36455)], OD.med.gind[which(OD.med.pmed < .05/36455)])
abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])$lengths), col='grey')
abline(h=cumsum(rle(gene.annot.df$chr)$lengths), col='grey')


sOD=scale(residuals(lm(OD.cov~gbatch.fact)))
ODmat=cbind(scale(OD.cov), sOD)

LODS_od=fasterLOD(nrow(pheno.scaled),ODmat,gdata.scaled)

plot(LODS_od[2,], type='l', ylab='LOD')
points(LODS_od[1,], type='l', col='blue')
abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])$lengths), col='grey')


#----------------------------------------------------------------------------------------------



     #testing if gene expression level differences mediate QTL effects on growth --------------------
     sy=scale(t.tpm.matrix[,gind])
     sOD=scale(OD.cov)
         
     M1= lm(sy  ~ gbatch.fact+gdata.scaled[,ppgl$pmarker]-1)
     M2= lm(sOD ~ gbatch.fact+gdata.scaled[,ppgl$pmarker]-1)
     M3= lm(sy  ~ gbatch.fact+sOD+gdata.scaled[,ppgl$pmarker]-1)
     sM1= summary(M1)$coefficients
     sM2= summary(M2)$coefficients
     sM3= summary(M3)$coefficients
     
     #https://en.wikipedia.org/wiki/Sobel_test
     tau=sM1[which(grepl('gdata', rownames(sM1))),'Estimate']
     tau.se=sM1[which(grepl('gdata', rownames(sM1))),'Std. Error']

     a=sM2[which(grepl('gdata', rownames(sM2))),'Estimate']
     a.se=sM2[which(grepl('gdata', rownames(sM2))),'Std. Error']

     b=sM3[match('sOD', rownames(sM3)),'Estimate']
     b.se=sM3[match('sOD', rownames(sM3)),'Std. Error']    
    
     tau.prime=sM3[which(grepl('gdata', rownames(sM3))),'Estimate']
     tau.prime.se=sM3[which(grepl('gdata', rownames(sM3))),'Std. Error']    

        #Arorian's second order exact solution (most conservative)
     pooled.se = a^2 * b.se^2 + b^2*a.se^2 + a.se^2+b.se^2
     Zscore=(tau-tau.prime)/pooled.se
     p.med=pnorm(Zscore, lower.tail=T)
     dff=data.frame(tau, tau.se, tau.prime, tau.prime.se, a,a.se,b,b.se,Zscore, p.med)
     rownames(dff)=ppgl$pmarker
     #----------------------------------------------------------------------------------------------
     OD.mediation[[gs]]=data.frame(ppgl, dff)
}






     dmatrix=data.frame(sOD, gbatch.fact, gdata.scaled[,ppgl$pmarker])
     dymatrix=data.frame(sy, OD.cov, gbatch.fact, gdata.scaled[,ppgl$pmarker])


     summary(Zm)$coefficients
    
     #Zm=lm(sOD~ gbatch.fact + chrI.62654_C.T  + chrIV.749852_G.A + 
     #      chrXI.225168_G.A + chrXII.646705_A.G + chrXII.135031_G.A+ 
     #      chrXIV.376313_C.T + chrXIV.456295_T.C + chrXV.141631_G.A -1, data=dymatrix)
     #Zy=lm(sy~ sOD + gbatch.fact + chrI.62654_C.T  + chrIV.749852_G.A + 
     #      chrXI.225168_G.A + chrXII.646705_A.G + chrXII.135031_G.A+ 
     #      chrXIV.376313_C.T + chrXIV.456295_T.C + chrXV.141631_G.A -1, data=dymatrix)
    
     #Zm=lm(sOD~ gbatch.fact + chrI.62654_C.T  + chrIV.749852_G.A + 
     #      chrXI.225168_G.A + chrXII.646705_A.G + chrXII.135031_G.A+ 
     #      chrXIV.376313_C.T + chrXIV.456295_T.C + chrXV.141631_G.A -1, data=dymatrix)
     #Zy=lm(sy~ sOD + gbatch.fact + chrI.62654_C.T  + chrIV.749852_G.A + 
     #      chrXI.225168_G.A + chrXII.646705_A.G + chrXII.135031_G.A+ 
     #s      chrXIV.376313_C.T + chrXIV.456295_T.C + chrXV.141631_G.A -1, data=dymatrix)
     
     mm=mediate(Zm, Zy, treat='chrXV.141631_G.A', mediator='sOD', data=dymatrix)
     
     
     
     


summary(Zm)$coefficients

sy~gbatch.fact+gdata.scaled[,ppgl$pmarker]-1  
     



     Z0=lm(scale(t.tpm.matrix[,gind])~gbatch.fact+gdata.scaled[,ppgl$pmarker]-1)
     Z2=lm(scale(t.tpm.matrix[,gind])~gbatch.fact+gdata.scaled[,ppgl$pmarker]+OD.cov)
     

       Z=

     M0=lm(scale(t.tpm.matrix[,gind])~gbatch.fact+gdata.scaled[,ppgl$pmarker]-1)
     M2=lm(scale(t.tpm.matrix[,gind])~gbatch.fact+gdata.scaled[,ppgl$pmarker]+OD.cov)
     
     M0.c=coef(M0)[14:length(coef(M0))] 
     M0.c.se=sqrt(diag(vcov(M0)))[14:length(coef(M0))] 
     M2.c=coef(M2)[-length(coef(M2))][14:length(coef(M0))]
     M2.c.se=sqrt(diag(vcov(M2)))[14:length(coef(M0))] 

     M2.minus.M0=M2.c-M0.c 
     names(M2.minus.M0)=colnames(gdata.scaled[,ppgl$pmarker])

     L0=lm(scale(OD.cov)~gbatch.fact+gdata.scaled[,ppgl$pmarker]-1)
     L2=lm(scale(OD.cov)~gbatch.fact+gdata.scaled[,ppgl$pmarker]+scale(t.tpm.matrix[,gind])-1)
     L2.minus.L0=coef(L2)[-length(coef(L2))][14:length(coef(L0))] - coef(M0)[14:length(coef(L0))] 
     names(L2.minus.L0)=colnames(gdata.scaled[,ppgl$pmarker])

      
    
     M1.minus.M0=coef(M1)[-length(coef(M2))][14:length(coef(M0))] - coef(M0)[14:length(coef(M0))] 

    mtest=mediate(M1,Q1)
    OD.cov
}






spbc=svd(p.batch.cis)
pigc=rep(NA,1012)
for(i in 1:1012) { pigc[i]=cor.test(spbc$u[,i], OD.cov)$p.value } 

spbco=svd(p.batch.cis.OD)
x11()
barplot(spbco$d/sum(spbco$d))
pigco=rep(NA,1012)
for(i in 1:1012) { pigco[i]=cor.test(spbco$u[,i], OD.cov)$p.value } 




most.h2.bc=names(which(h2.bc>.19))
p.batch.cis.h2thresh=p.batch.cis[,match(most.h2.bc, colnames(p.batch.cis))]

lambda1_seq=c(0.08, 0.04, 0.02, 0.01)
lambda2_seq=c(0.32, 0.16, 0.08, 0.04, 0.02)

#lambda1_seq = rev(seq(.001,.3, .1))
#lambda2_seq = rev(seq(.001,.3, .1))
tpm.ds=p.batch.cis.h2thresh[,sort(sample(1:3166, 200))]

opt_cv = scggm_cv(g.h, p.batch.cis.h2thresh, 3, lambda1_seq, lambda2_seq)

# values around .05 and .05 gave reasonable results
ot_05_05=scggm(g.h, p.batch.cis.h2thresh, 0.05,  0.05)




# 100 last round  0.051, optimal lambda2 = 0.251"
opt_cv = scggm_cv(g.h, tpm.ds, 3, lambda1_seq, lambda2_seq)






'YCR039C'
'YCR040W'
'YHR084W'

pm=match('YHR084W', colnames(p))
pm=match('YHR005C', colnames(p))

#
hh=12


load(
coef(lm(p[,mediator]~g.h[,hh]))
summary((lm(p[,mediator]~g.h[,hh]+p[,mediated])))


load('/home/jbloom/Dropbox/Public/eQTL/hotspot_coefficients_lasso.RData')
#lasso.matrix
 

mat.hotspot.transcripts=names(which(abs(lasso.matrix[13,])>0))
tmat=p[,match(mat.hotspot.transcripts, colnames(p))]
X1=model.matrix(tmat[,1]~g.h[,12])
#these are mediated
lm1=lm.fit(X1,tmat)

lm2=list()
for( tr in 1:ncol(tmat) ){
    print(tr)
    X2=model.matrix(tmat[,1]~g.h[,12]+p[,tr])
    lm2[[colnames(tmat)[tr]]]=coef(lm.fit(X2, tmat,singular=FALSE))
}

lm2.mat=do.call('rbind', lapply(lm2, function(x) x[2,]))
diag(lm2.mat)=NA
# rows are mediator
# columns are mediated

#lm1.c=coef(lm1)
#med.mat=apply(lm2.mat, 2, function(x) x-lm1.c[1,] ) 
#hist(as.vector(coef(lm1)))

#lm2.mat[1:10,1:10]
#as.vector(coef(lm1))[1:10]
#lm3=

#hist(lm2.mat[,1401], breaks=100)


colnames(tmat)[1401]

hist((lm2.mat[1401,]^2) ,breaks=100)
abline(v=(lm1.c[1,1401]^2))

plot( lm2.mat[1401,], lm1.c[1,])

plot(lm1.c[1,], lm2.mat[,1], xlab='m1', ylab='m2' )


for(i in 2:1400) 


    transp=rgb(0,0,0,.1)
plot(lm2.mat[1,], lm1.c[1,], xlab='m2', ylab='m1' ,xlim=c(-.25, .25), ylim=c(-.25, .25))
for(i in 2:500) {
    print(i)
points(lm2.mat[i,], lm1.c[1,], xlab='m2', ylab='m1', col=transp )
}

abline(0,1, col='red')
plot(lm2.mat[690,], lm1.c[1,], xlab='m2', ylab='m1' )
colnames(tmat)[690]
 
tr=690
X2=model.matrix(tmat[,1]~g.h[,12]+p[,tr]-1)
lcrz=lm.fit(X2, tmat,singular.ok=FALSE)
lm(


#mediated ste2 or 3 or mfa1
#mediator matalpha2 or matalpha1 
#matalpha
mediator=match('YCR040W', colnames(p))
#mfa1
#mediated=match('YDR461W', colnames(p))

#ste2
mediated=match('YFL026W', colnames(p))

coef(lm(p[,mediated]~g.h[,hh]))
summary((lm(p[,mediated]~g.h[,hh]+p[,mediator])))

coef(lm(p[,mediator]~g.h[,hh]))
summary((lm(p[,mediator]~g.h[,hh]+p[,mediated])))


coefs=list()
for( i in c(1:6288)[-pm] ) {
    coefs[[colnames(p)[i]]]=coef(lm(p[,pm]~g.h[,hh]+p[,i]))
}
coefs=do.call('rbind', coefs)

plot(coefs[,2], coefs[,3], xlab='B_31(QTL)', ylab='B_32(Mediator)')
abline(v=(coef(lm(p[,pm]~g.h[,hh]))[2]))
identify(coefs[,2], coefs[,3],labels=rownames(coefs))

rownames(coefs)[2984]
peaks.per.chr=split(all.peaks, all.peaks$chr)
#library(EBglmnet)
library(caret)

source('~/Local/SCGGM_code/scggm.R')

for(cc in unique.chrs[-17][1:16]) { 

   cc=unique.chrs[1]
   ppc=peaks.per.chr[[cc]]
   print(cc)
   keep.transcripts=colnames(t.tpm.matrix)
   # remove local QTL
   #keep.transcripts=keep.transcripts[gene.annot.df$chr!=cc]  
   ppc2=ppc[ppc$gene %in% keep.transcripts,]
   keep.transcripts=ppc2$gene
   keep.transcripts=unique(keep.transcripts)
  
   tmm=matrix(0, nrow(t.tpm.matrix), length(keep.transcripts))
   mid=match(keep.transcripts, colnames(t.tpm.matrix))
   colnames(tmm)=colnames(t.tpm.matrix)[mid]
   rownames(tmm)=rownames(t.tpm.matrix)
   #Correct transcripts for background QTL effects 
   q.cnt=rep(0, length(keep.transcripts))
   q.pos=list()
   names(q.cnt)=colnames(tmm)
   for(i in 1:length(keep.transcripts)) { 
         #print(i)
         mid.i= mid[i]
         gene.i=colnames(tmm)[i]
         ppgi.a=peaks.per.gene[[gene.i]]
         q.cnt[i]=sum(ppgi.a$chr==cc)
         q.pos[[gene.i]]=ppgi.a$pcind[ppgi.a$chr==cc]
         ppgi=ppgi.a[ppgi.a$chr!=cc,]
         if(is.null(ppgi)) {
             tmm[,i]=scale(residuals(lm(t.tpm.matrix[,mid.i]~gbatch.fact))) 
         }else {
             if(nrow(ppgi)==0) {
                 tmm[,i]=scale(residuals(lm(t.tpm.matrix[,mid.i]~gbatch.fact)))
             next;
             } else {
                 bQTL=gdata[,ppgi$gcind[ppgi$chr!=cc]]
                 tmm[,i]=scale(residuals(lm(t.tpm.matrix[,mid.i]~gbatch.fact+bQTL) ) ) 
             }
         }
   }

   presid=tmm
   # avoid unecessary scaling
   # presid=scale(preal)
   Aloco=A.mat(do.call('cbind', (gdata.by.chr[-match(cc, names(gdata.by.chr))])))/2
             
    # Eskin trick to speed up variance component calculation
    eigA=doEigenA_forMM(presid,Aloco)
    # calculate mixed model, one term for additive variance  -------------------------------------------
    pb=txtProgressBar(min=1, max=ncol(presid), style=3)
    for(tp in 1:ncol(presid)){
         setTxtProgressBar(pb,tp)
         #print(tp)
         rr=m.S(presid[,tp], K=Aloco,  theta=eigA$theta, Q=eigA$Q)
         W=solve(rr[1]*Aloco+rr[2]*diag(nrow(presid)))
         if(rr[1]>0) {
             blups=calc.BLUPS(rr[1]*Aloco,diag(nrow(presid)),W,presid[,tp],matrix(1,nrow(presid),1),0 )[,1]
             #presid[,tp]=as.vector(scale(presid[,tp]-blups)[,1])
             presid[,tp]=as.vector(presid[,tp] - blups)

         }
     }
     rm(W)
     close(pb)
     #
     #preal=bresids[[cc]]
     #presid=scale(preal)
     tpm_gQA=scale(presid)
    
     fL=fasterLOD(nrow(tpm_gQA),  tpm_gQA,  gdata.scaled) 
    # abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), ':'))[,1])$lengths), col='grey')

 #   abline(v=cumsum(rle(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])$lengths), col='grey')

     g=gdata.s.by.chr[[cc]]
     XXr=crossprod(g)/(ncol(tpm.matrix)-1)
     #XYr=crossprod(g,tmm)/(ncol(tpm.matrix)-1)
     #scanone.all=(-1012*log(1-XYr^2))/(2*log(10))
     fCX=sort(findCorrelation(XXr, cutoff=.99, exact=F))
     g.s=g[,-fCX]

# lambda1 is XY and lambda2 ix YY
#       = 0.007, optimal lambda2 = 0.004"
# search 1 
     #.001 to .01

# search 2      
# "lambda1 = 0.001, lambda2 = 0.001, cross validation error = 0.951074312198161"
#[1] "training sCGGM with optimal regularization parameters:"
#[1] "optimal lambda1 = 0.021, optimal lambda2 = 0.041"


#lambda1_seq = rev(seq(.001,.1, .02))
#lambda2_seq = rev(seq(.001,.1, .02))

lambda1_seq = rev(seq(.001,.1, .01))
lambda2_seq = rev(seq(.001,.5, .05))
tpm.ds=tpm_gQA[,sort(sample(1:899, 300))]

# 100 last round  0.051, optimal lambda2 = 0.251"
opt_cv = scggm_cv(g.s,tpm.ds, 3, lambda1_seq, lambda2_seq)

which(abs(opt_cv$Theta$xy)>0)
which((abs(opt_cv$Theta$yy)>0))


opt_cv = scggm_cv(g.s,tpm_gQA, 3, lambda1_seq, lambda2_seq)

save(opt_cv, file='/data/eQTL/RData/chrI_scggm.RData')

scggm_out_chrIII=scggm(g.s, tpm_gQA, 0.021,  0.041)
scggm_out_chrIII.iso=scggm_indirect_SNP_overall(scggm_out_chrIII$Theta)
scggm_out_chrIII.cov_decomp=scggm_cov_decompose(scggm_out_chrIII$Theta, g.s, tpm_gQA)
isd=list()
for(i in 1 :ncol(tpm_gQA)) {
    isd[[i]]=scggm_indirect_SNP_decompose(scggm_out_chrIII$Theta,i)
}

save(file='/data/eQTL/RData/chrIII_scggm.RData'
     , list=c('tpm_gQA', 'g.s', 'scggm_out_chrIII', 'scggm_out_chrIII.iso', 'scggm_out_chrIII.cov_decomp', 'isd'))


# with a subset of 100
scggm_out_chrI=scggm(g.s, tpm_gQA, 0.041,  0.081)

scggm_out_chrI.iso=scggm_indirect_SNP_overall(scggm_out_chrI$Theta)
scggm_out_chrI.cov_decomp=scggm_cov_decompose(scggm_out_chrI$Theta, g.s, tpm_gQA)
isd=list()
for(i in 1 :ncol(tpm_gQA)) {
    isd[[i]]=scggm_indirect_SNP_decompose(scggm_out_chrI$Theta,i)
}
#note these will clobber variables in your workspace
save(file='/data/eQTL/RData/chrI_scggm.RData'
     , list=c('tpm_gQA', 'g.s', 'scggm_out_chrI', 'scggm_out_chrI.iso', 'scggm_out_chrI.cov_decomp', 'isd'))


#fit lasso for each trait 
tpm_gQA
tpm.lasso=list()
for(i in 77:899)  {
tpm.lasso[[i]]=cv.EBglmnet(g.s,tpm_gQA[,i])
}

lasso.matrix=matrix(0, ncol(g.s), ncol(tpm_gQA) ) 
for(i in 1:899) {
    print(i)
    lh=tpm.lasso[[i]]
    lasso.matrix[lh$fit[,1],i]=lh$fit[,3]
}
B=lasso.matrix
 xy =B %*% cov(tpm_gQA)






lasso.hotspots=list()
for(i in 391:ncol(bc.resid)) {
    print(i)
    lasso.hotspots[[colnames(bc.resid)[i]]]=cv.EBglmnet(gdata.s[,hotspot.vector], bc.resid[,i])
}








(scggm_indirect_SNP_decompose(scggm_out_chrIII$Theta,1))


<- function(Theta, k){
    iThetayy = solve(Theta$yy)
    Beta_k = -Theta$xy[,k] %*% t(iThetayy[,k])
    return(Beta_k)
}

opt_cv 
system.time({
        test.e.scggm <- scggm(g.s, tpm_gQA, 0.01, 0.01)
    })


     


      
      
       #sXXr=crossprod(g.s)/(ncol(tpm.matrix)-1)
       #sXYr=crossprod(g.s,tmm)/(ncol(tpm.matrix)-1)
       #scanone.red=(-1012*log(1-sXYr^2))/(2*log(10))
       #cor(apply(scanone.all,2,max), apply(scanone.red, 2, max))
       #plot(apply(scanone.all,2,max), apply(scanone.red, 2, max), xlim=c(0,40), ylim=c(0,40) )
       for( i in 1:ncol(tmm)) {
          # print(i)
           # with preset hyperparameters
          # optimize this 
           bl[[cc]][[colnames(tmm)[i]]]=EBglmnet(g.s, tmm[,i], hyperparameters=c(.1,.1),verbose=0)
           # choose hyperparameters
           #bl[[cc]][[colnames(tmm)[i]]]=cv.EBglmnet(g.s, tmm[,i],verbose=-1)
       }

       gind=match(names(bl[[cc]]), gene.annot.df[,1])
       q.bl.c=lapply(bl[[cc]], function(x) nrow(x$fit))
       q.bl.c=do.call('c', q.bl.c)
       greps=rep(gind,as.vector(q.bl.c))
       bl.loc=sapply(bl[[cc]], function(x) x$fit[,1] )
       bl.loc=sapply(bl[[cc]], function(x) x$fit[,1] )
       bl.beta=sapply(bl[[cc]], function(x) x$fit[,'beta'] )
       bl.var=sapply(bl[[cc]], function(x) x$fit[,'posterior variance'] )
       bl.t=sapply(bl[[cc]], function(x) x$fit[,'t-value'] )
       bl.p=sapply(bl[[cc]], function(x) x$fit[,'p-value'] )
       #plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.beta)))
       #plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.var)))
       #plot(unlist(as.vector(bl.loc)), unlist(as.vector(bl.t)))
       #plot(unlist(as.vector(bl.loc)), -log10(unlist(as.vector(bl.p))))
 }


   tmp= fasterLOD(nrow(tmm),tmm,gdata.s.by.chr[[cc]])
   max.obsLOD=apply(tmp,1,max)   
   n.perm=100
   permLOD=matrix(0, length(max.obsLOD), n.perm)
   print('doing permutations')
   pb =txtProgressBar(min = 1, max =n.perm, style = 3)
   for(i in 1:n.perm) {
             setTxtProgressBar(pb, i)
             permLOD[,i]=apply(fasterLOD(nrow(tmm),tmm[sample(1:nrow(tmm)),], gdata.s.by.chr[[cc]]),1,max)
   }
   close(pb)
       
   obsPcnt = sapply(seq(1.5, 6, .01), function(thresh) { sum(max.obsLOD>thresh) }   )
   names(obsPcnt) = seq(1.5, 6, .01)
       
   if(sum(obsPcnt)<5) {break;}
               
   # expected number of QTL peaks with LOD greater than threshold
   expPcnt = sapply(seq(1.5, 6, .01),  
                         function(thresh) { 
                                #print(thresh); 
                                mean(apply(permLOD, 2, function(ll) {sum(ll>thresh) }) )
                            } )
   names(expPcnt) = seq(1.5, 6, .01)
   pFDR = expPcnt/obsPcnt
   pFDR = rev(cummax(rev(pFDR)))
   if(sum(is.na(pFDR))>2) {break;}
   fdrFX=approxfun(pFDR, seq(1.5,6,.01))
   thresh=fdrFX(.05)
   print(paste('FDR Thresh 5%', thresh)) 
  
   #thresh
   keep.transcripts=names(which(max.obsLOD>thresh))
   print(length(keep.transcripts))
   
   if(length(keep.transcripts)<10) {break;}

   LOD.kt=tmp[keep.transcripts,]
   LOD.kt.wmax=apply(tmp[keep.transcripts,], 1, which.max)
   LOD.kt.max=apply(tmp[keep.transcripts,], 1, max)
   #plot(apply(LOD.kt,2,max))
   #plot(apply(LOD.kt,2,sum))

   bb=seq(0,ncol(tmp),30)
   if(max(ncol(tmp))>max(bb)) { bb=c(bb, max(bb)+30) }
   l.in.bins=cut(LOD.kt.wmax, bb)
   r.l.in.bins=rle(sort(as.vector(l.in.bins)))
   max.bin=r.l.in.bins$values[which.max(r.l.in.bins$lengths)]

   mb=gsub('\\(', '', max.bin)
   mb=gsub('\\]', '', mb)
   mb=as.numeric(strsplit(mb, ',')[[1]])

   g.s.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g.s),':'), function(x)x[2])), '_'), function(x) x[1]))
   g.f.pos=  as.numeric(sapply(strsplit((sapply(strsplit(colnames(g),':'), function(x)x[2])), '_'), function(x) x[1]))

   to.scan=names(LOD.kt.max)[which(l.in.bins %in% max.bin)]
   max.to.scan=LOD.kt.max[to.scan]
   max.to.scan=(sort(max.to.scan, decreasing=T))
   #max.to.scan=max.to.scan[max.to.scan>3.5]
   to.scan=names(max.to.scan)
   #nn=min(800, length(keep.transcripts))
   #if(cc=='chrXIV' & nn>799) {nn=500}
   #if(cc=='chrXIV' & nn>799) {nn=500}
   #nn=50
   if(length(to.scan)>500) { to.scan=to.scan[1:500] }
   if(length(to.scan)<3)  {break;}
   
   nn=length(to.scan)
   #sub.t=tmm[,names(sort(max.obsLOD,decreasing=T)[1:nn])]
   sub.t=tmm[,to.scan]
   #regress out effects of background QTL on chromosome

   if(mb[2]>ncol(g)) {mb[2]=ncol(g) }
   if(mb[1]>1) {mb[1]=1 }

   for(ttt in colnames(sub.t)) {
       #print(ttt)
        bcl=bl[[cc]][[ttt]]$fit[,'locus1']
        bq=bcl[g.s.pos[bcl]<g.f.pos[mb[1]] | g.s.pos[bcl]>g.f.pos[mb[2]]]

        if(length(bq)>0) {
            print(length(bq))
        sub.t[,ttt]=scale(residuals(lm(sub.t[,ttt]~g.s[,bq])))[,1] 
        }
   }
  
   #  save(gdata.s.by.chr, file='/media/juno/bootstrap_hotspot/gdata.s.by.chr.RData' )
    #what information do we need to do the hotspot bootstrapping? 
   #1) sub.t (de-regressed phenotype values)
   #2) gdata.s.by.chr  and the target chromosome 
   #3) that is it ... sample from the number of segregants with replacment, recalculate across chromosome, record max, and repeat 
    

   #library(Brobdingnag)
   sss=rep(NA, ncol(gdata.s.by.chr[[cc]]))
   pb =txtProgressBar(min = 1, max =length(sss), style = 3)
   
   #optimize ... SST doesn't need to be recalculated during this scan!!!!! at least 2X speedup 
   SST=crossprod((sub.t))
   t5=svd(SST)  

   pc.to.retain=which(cumsum(t5$d/sum(t5$d))<.999)
   print(length(pc.to.retain))
   #  http://arxiv.org/abs/1305.5870
   for(i in 1:length(sss)){
       setTxtProgressBar(pb, i)
       test=lm(sub.t~gdata.s.by.chr[[cc]][,i])
       t2=residuals(test)
       SSE=crossprod(t2)
       t4=svd(SSE)
       X2=-(1012-1-((nn+2)/2))*(sum(log(t4$d[pc.to.retain]))-sum(log(t5$d[pc.to.retain])))
       sss[i]=X2
   }
   close(pb)
   #plot(-log10(pchisq(sss,nn,lower.tail=F)))
   #sss=-log10(pchisq(sss,nn,lower.tail=F))
   # plot(sss)
   max(sss)
   if(!is.finite(max(sss)))  {break;}
   if(min(sss)< (-1)) {break;}
   if((max(sss)/4.605)<5)  {break;}
   #mvlod[[cc]][[hs.n]]=sss
   names(sss)=colnames(gdata.s.by.chr[[cc]])
   #rl.w=rle(which(sss>(max(sss)-1.5 )))
   hs.int.range=range(which((sss/4.605)>(max(sss/4.605)-1.5 ))) 
   names(sss)[hs.int.range]
   
   plot(sss/4.605, col='red')
   pmarker=c(pmarker, names(which.max(sss)))
   l=names(sss)[hs.int.range][1]
   r=names(sss)[hs.int.range][2]

   hotspots[[cc]][[hs.n]]=data.frame(maxneglog10p=max(sss), n.trans=length(keep.transcripts), peak=names(which.max(sss)), cI.l=l, cI.r=r)
   total.hotspot.n=total.hotspot.n+1
   attr(sub.t, 'chr')=cc
   attr(sub.t, 'hs.n')=total.hotspot.n
   save(sub.t, file=paste0('/media/juno/bootstrap_hotspot/', total.hotspot.n, '.RData' ))

   print(cc)
   print(hs.n)
   print(hotspots[[cc]][[hs.n]])   
   hs.n=hs.n+1
  }
}

