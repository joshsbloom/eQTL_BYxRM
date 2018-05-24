library(foreach)
library(doMC)
registerDoMC(cores=6)

# BUILD KALLISTO INDEX ------------------------------------------------------------
# run once only within reference/
#kallisto index -i transcripts.idx Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
#-----------------------------------------------------------------------------------


root_dir='/data/eQTL/'
transcript.idx=paste0(root_dir,  'reference/transcripts.idx')
fastq_dir=paste0(root_dir, 'fastq_trimmed/')
kallisto_dir=paste0(root_dir, 'kallisto_out_merged_stranded_boot/')
#kallisto_dir=paste0(root_dir, 'kallisto_out_merged_stranded_l150_s20/')
dir.create(kallisto_dir)
fastq.files=list.files(fastq_dir)
f.size=file.size(paste0(fastq_dir, fastq.files))

fmat=do.call('rbind', strsplit(fastq.files,'-'))
colnames(fmat)=c('sample', 'cross', 'well', 'growth', 'growth.well', 'lane')
fmat=data.frame(fmat, f.size,seq(1:nrow(fmat)), stringsAsFactors=F)
colnames(fmat)[8]='index'
f.split=split(fmat, fmat$sample)
#t2=test[sapply(test,nrow)>2]
t3=lapply(f.split, function(x) split (x, x$growth))
best.batch=sapply(t3, function(x) names(which.max(sapply(x, function(y) sum(y$f.size)))))

foreach(samp=names(best.batch[-c(1,2)])) %dopar% {
    imat=f.split[[samp]]
    imat=imat[imat$growth==best.batch[samp],]
    out.name=paste(imat$sample[1], imat$cross[1], imat$well[1], imat$growth[1], imat$growth.well[1], sep='-')
    input.indices=imat$index
    in.files=paste0(fastq_dir, fastq.files[input.indices], collapse= ' ' )
    print(samp)
    print(out.name)
    print(in.files)
    out.file=paste0(kallisto_dir, out.name)
    # original v42.1
    #system(paste('/home/jbloom/Local/bin/kallisto quant --single -l 150 -i', transcript.idx, '-o', out.file, in.files))
    # merged stranded boot 
    system(paste('/home/jbloom/Local/kallisto_linux-v0.43.0/kallisto quant -t 10 --rf-stranded --single -l 150 -s 8 -b 100 -i', transcript.idx, '-o', out.file, in.files))
    # not stranded, but new version
    # merged no boot not stranded
    #system(paste('/home/jbloom/Local/kallisto_linux-v0.43.0/kallisto quant -t 10 --single -l 150 -s 8 -i', transcript.idx, '-o', out.file, in.files))
    # merged different expected fragment length
    #system(paste('/home/jbloom/Local/kallisto_linux-v0.43.0/kallisto quant -t 10 --rf-stranded --single -l 100 -s 10 -i', transcript.idx, '-o', out.file, in.files))
    #system(paste('/home/jbloom/Local/kallisto_linux-v0.43.0/kallisto quant -t 10 --rf-stranded --single -l 150 -s 20 -i', transcript.idx, '-o', out.file, in.files))
}
#-------------------------------------------------------------------------------------------------------------------------------------


# Run kallisto on merged parental values----------------------------------------------------------------------------
library(foreach)
library(doMC)
registerDoMC(cores=12)

root_dir='/data/eQTL/'
transcript.idx=paste0(root_dir,  'reference/transcripts.idx')
fastq_dir=paste0(root_dir, 'fastq_trimmed/')
kallisto_dir=paste0(root_dir, 'kallisto_out_merged_parents_stranded_boot/')
dir.create(kallisto_dir)

fastq.files=list.files(fastq_dir)
f.size=file.size(paste0(fastq_dir, fastq.files))

fmat=do.call('rbind', strsplit(fastq.files,'-'))
colnames(fmat)=c('sample', 'cross', 'well', 'growth', 'growth.well', 'lane')
fmat=data.frame(fmat, f.size,seq(1:nrow(fmat)), stringsAsFactors=F)
colnames(fmat)[8]='index'

fmat=fmat[grep('1879|1950', fmat$sample),]

f.split=split(fmat, paste(fmat$sample, fmat$well, sep='_') )
#t2=test[sapply(test,nrow)>2]
t3=lapply(f.split, function(x) split (x, x$growth))
best.batch=sapply(t3, function(x) names(which.max(sapply(x, function(y) sum(y$f.size)))))

foreach(samp=names(best.batch)) %dopar% {
    imat=f.split[[samp]]
    imat=imat[imat$growth==best.batch[samp],]
    out.name=paste(imat$sample[1], imat$cross[1], imat$well[1], imat$growth[1], imat$growth.well[1], sep='-')
    input.indices=imat$index
    in.files=paste0(fastq_dir, fastq.files[input.indices], collapse= ' ' )
    print(samp)
    print(out.name)
    print(in.files)
    out.file=paste0(kallisto_dir, out.name)
#    system(paste('/home/jbloom/Local/bin/kallisto quant --single -l 150 -i', transcript.idx, '-o', out.file, in.files))
    system(paste('/home/jbloom/Local/kallisto_linux-v0.43.0/kallisto quant -t 10 --rf-stranded --single -l 150 -s 8 -b 100 -i', transcript.idx, '-o', out.file, in.files))
}

#---------------------------------------------------------------------------------------------------------------------------------------------





# Running alignments with STAR ---------------------------------------------------------------------------------------------------------------
## For segregants run STAR alignments
#root_dir='/data/eQTL/'
#transcript.idx=paste0(root_dir,  'reference/transcripts.idx')
#fastq_dir=paste0(root_dir, 'fastq_trimmed/')
#star_dir_p=paste0(root_dir, 'star_out_parents/')
#star_dir=paste0(root_dir, 'star_out/')
#
#dir.create(star_dir_p)
#
#fastq.files=list.files(fastq_dir)
#f.size=file.size(paste0(fastq_dir, fastq.files))
#
#fmat=do.call('rbind', strsplit(fastq.files,'-'))
#colnames(fmat)=c('sample', 'cross', 'well', 'growth', 'growth.well', 'lane')
#fmat=data.frame(fmat, f.size,seq(1:nrow(fmat)), stringsAsFactors=F)
#colnames(fmat)[8]='index'
#
#fmat.p=fmat[grep('1879|1950', fmat$sample),]
#f.split.p=split(fmat.p, paste(fmat.p$sample, fmat.p$well, sep='_') )
##t2=test[sapply(test,nrow)>2]
#t3.p=lapply(f.split.p, function(x) split (x, x$growth))
#best.batch.p=sapply(t3.p, function(x) names(which.max(sapply(x, function(y) sum(y$f.size)))))
#
#foreach(samp=names(best.batch.p)) %do% {
#    imat=f.split.p[[samp]]
#    imat=imat[imat$growth==best.batch.p[samp],]
#    out.name=paste(imat$sample[1], imat$cross[1], imat$well[1], imat$growth[1], imat$growth.well[1], sep='-')
#    input.indices=imat$index
#    in.files=paste0(fastq_dir, fastq.files[input.indices], collapse= ' ' )
#    print(samp)
#    print(out.name)
#    print(in.files)
#    fq.merged.file=tempfile()
#    system(paste('cat', in.files, '>', fq.merged.file))
#    out.file=paste0(star_dir_p, out.name, '/')
#    dir.create(out.file)
#    system(paste('STAR --genomeDir /data/eQTL/reference/STAR --runThreadN 10 --readFilesIn', fq.merged.file, '--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ', out.file, '--outReadsUnmapped Fastx --outWigType bedGraph --outWigNorm None' ))
#    unlink(fq.merged.file)
#     #majority of reads are accounted for here
#    #system(paste('express --r-stranded /data/eQTL/reference/R64-1.1.cdna.ncrna.fa ',  out.file, 'Aligned.toTranscriptome.out.bam -o  ', out.file, 'F/', sep='') )
#    #system(paste('express --f-stranded /data/eQTL/reference/R64-1.1.cdna.ncrna.fa ',  out.file, 'Aligned.toTranscriptome.out.bam -o  ', out.file, 'R/', sep='') )
#}
#
#f.split=split(fmat, fmat$sample)
##t2=test[sapply(test,nrow)>2]
#t3=lapply(f.split, function(x) split (x, x$growth))
#best.batch=sapply(t3, function(x) names(which.max(sapply(x, function(y) sum(y$f.size)))))
#
#bb=sapply(t3, function(x) {sapply(x, function(y) sum(y$f.size))})
#bb.reps=bb[sapply(bb,length)>1]
#bb.reps=bb.reps[sapply(bb.reps, function(x) x[1]>1e7 & x[2]>1e7 )][-c(1,2)]
##t(do.call('cbind', bb[sapply(bb, length)>1]))
#
#for( samp in names(best.batch)[419:1058]) { 
#    #[each(samp=names(best.batch[-c(1,2)])) %do% {
#    imat=f.split[[samp]]
#    imat=imat[imat$growth==best.batch[samp],]
#    out.name=paste(imat$sample[1], imat$cross[1], imat$well[1], imat$growth[1], imat$growth.well[1], sep='-')
#    input.indices=imat$index
#    in.files=paste0(fastq_dir, fastq.files[input.indices], collapse= ' ' )
#    print(samp)
#    print(out.name)
#    print(in.files)
#    fq.merged.file=tempfile()-
#    system(paste('cat', in.files, '>', fq.merged.file))
#    out.file=paste0(star_dir, out.name, '/')
#    dir.create(out.file)
#    system(paste('STAR --genomeDir /data/eQTL/reference/STAR --runThreadN 10 --readFilesIn', fq.merged.file, '--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ', out.file, '--outReadsUnmapped Fastx --outWigType bedGraph --outWigNorm None' ))
#    unlink(fq.merged.file)
#    #system(paste('express --r-stranded /data/eQTL/reference/R64-1.1.cdna.ncrna.fa ',  out.file, 'Aligne-d.toTranscriptome.out.bam -o  ', out.file, 'F/', sep='') )
#    #system(paste('express --f-stranded /data/eQTL/reference/R64-1.1.cdna.ncrna.fa ',  out.file, 'Aligned.toTranscriptome.out.bam -o  ', out.file, 'R/', sep='') )
#
#}
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
