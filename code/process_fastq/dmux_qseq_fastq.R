# Code to demultiplex gz fastq files given sample sheets
# Run for each batch and lane.identifier

#mkdir -p /media/ramdisk
#mount -t tmpfs -o size=32048M tmpfs /media/ramdisk/

#find . -name '*eQTL_13-*' -delete
#find . -name '*eQTL_13-*'  -type d -exec rm -r {} \;

library(ShortRead)
library(foreach)
library(doMC)
n.cores=6
#n.cores=1
registerDoMC(cores=n.cores)
cl <- makeCluster(getOption("cl.cores", n.cores))

#batch='eQTL1'; lane.identifier='SxaQSEQsVB100L6'  #KKHwuMHm 
#batch='eQTL2'; lane.identifier='SxaQSEQsVB100L7'  #JcnasWvQ 
#batch='eQTL4';  lane.identifier='SxaQSEQsVB101L1' #:dlxbzKRH 
#batch='eQTL5';  lane.identifier='SxaQSEQsVB101L2' #:NESeWaFM 
#batch='eQTL6';  lane.identifier='SxaQSEQsVB101L3' #:lXVwGmRA
#batch='eQTL7';  lane.identifier='SxaQSEQsVB101L4' #:DUJTqxg7 
#batch='eQTL8';  lane.identifier='SxaQSEQsVB101L5' #:YXD6pNhA 

#batch='eQTL9';  lane.identifier='SxaQSEQsVB101L6' #:lGAg6hBx
#batch='eQTL10'; lane.identifier='SxaQSEQsVB101L7' #:7hYU2zdn 
#batch='eQTL11'; lane.identifier='SxaQSEQsVB101L8' #:E2HMbB5W 

#batch='eQTL3'; lane.identifier='SxaQSEQsWB117L2'  #:aeavGPns 
#batch='eQTL12'; lane.identifier='SxaQSEQsWB117L3' #::ltU8LhKN 
#batch='eQTL13'; lane.identifier='SxaQSEQsWB117L4' #:Nhq44d3W 

#batch='eQTL1';   lane.identifier='SxaQSEQsWB117L5' #:nD3ZGWNn 
#batch='eQTL2';   lane.identifier='SxaQSEQsWB117L6' #:e4qPexws 
#batch='eQTL4';   lane.identifier='SxaQSEQsWB117L7' #:P3FjBEmr 

#batch='eQTL5';   lane.identifier='SxaQSEQsWB117L8' #:uMyEf9td 
#batch='eQTL6';  lane.identifier='SxaQSEQsVA093L1' #:dqh4ryxZ 
#batch='eQTL7';  lane.identifier='SxaQSEQsVA093L2' #:84fBkvjP 
#batch='eQTL8';  lane.identifier='SxaQSEQsVA093L3' #:ysBWtfgG
#batch='eQTL9';  lane.identifier='SxaQSEQsVA093L4' #:Uq3Velkz 

#batch='eQTL10'; lane.identifier='SxaQSEQsVA093L5' #:8tt3rU7C 
#batch='eQTL11'; lane.identifier='SxaQSEQsVA093L6' #:euuGeBeK 
#batch='eQTL12'; lane.identifier='SxaQSEQsVB103L1' #:vxXntz7B
#batch='eQTL13'; lane.identifier='SxaQSEQsVB103L2' #:Ntt9uC3r
batch='eQTL3' ; lane.identifier='SxaQSEQsXa091L1'

dcr=function(x) {do.call('rbind', x) } 
dcc=function(x) {do.call('cbind', x) } 

readformat= 'gz.fastq'
#readformat = 'qseq'
#base.dir='/media/gen/eQTL/'
#base.dir='/media/ramdisk/'
base.dir='/data/eQTL/'
out.base.dir = base.dir
#out.base.dir='/data/eQTL/'

sample.sheet.file = paste(base.dir, 'sample_sheets/sample_sheet_', batch, '.csv', sep='')
base.dir.qseq=(paste(base.dir, 'qseq/', batch, '/', sep=''))

#read in a miseq formatted sample sheet
ss=read.delim(sample.sheet.file, header=T, sep=',', stringsAsFactors=F) #skip=20 # if miseq sample sheet
n700.inds=unique(ss$index)
n500.inds=unique(ss$index2)

ind384=paste(match(ss$index, n700.inds), match(ss$index2, n500.inds), sep='_')

# Version starting with gzipped fastq files  --------------------------------------------------
i1=paste(base.dir.qseq, lane.identifier, '_index1.fastq.gz', sep='')
i2=paste(base.dir.qseq, lane.identifier, '_index2.fastq.gz', sep='')
r1=paste(base.dir.qseq, lane.identifier, '_read1.fastq.gz', sep='')

nbuffer=2e6

#PRELOADED to here
fi1 =FastqStreamer(i1, nbuffer, readerBlockSize=1e7,verbose=F)
fi2 =FastqStreamer(i2, nbuffer, readerBlockSize=1e7)
fr1 =FastqStreamer(r1, nbuffer, readerBlockSize=1e7)
#fr2 =FastqStreamer(r2, nbuffer, readerBlockSize=1e7)
repeat {
    fq.i1=yield(fi1)   
    fq.i2=yield(fi2)   
    fq.r1=yield(fr1)   
    #fq.r2=yield(fr2)   

    if(length(fq.i1) ==0 ) {break }
    print('calculating distance to index1')
    # calculate distance matrix of reads indices (read1 is N700 read)
    index1.dist.calc=srdistance(sread(fq.i1), n700.inds	)
    d1mat=do.call('cbind', index1.dist.calc)
    # find closest match for each read
    d1min.ind=parRapply(cl, d1mat, which.min)
    # find edit distance of closest match
    d1min=d1mat[cbind(1:nrow(d1mat), d1min.ind)]

    # lop off the last base
    # s7=DNAStringSet(sread(fq.i2),start=1, end=8)
    print('calculating distance to index2')
    # calculate distance matrix of reads indices (read2 is N500 read)
    index2.dist.calc=srdistance(sread(fq.i2), n500.inds)
    d2mat=do.call('cbind', index2.dist.calc) 
    # find closest match for each read
    d2min.ind=parRapply(cl, d2mat, which.min)
   # find edit distance of closest match
    d2min=d2mat[cbind(1:nrow(d2mat), d2min.ind)]

    keep.reads=which(d1min<2 & d2min<2)
    print(paste( (length(keep.reads)/nbuffer)*100 ,  ' % of reads retained'))
    d1.split.ind=d1min.ind[keep.reads]
    d2.split.ind=d2min.ind[keep.reads]
    matchtoinput=paste(d1.split.ind, d2.split.ind, sep='_')
    #readID=ss$Sample_Name[match(matchtoinput, ind384)]
    readID=paste(ss$Sample_Name, ss$Sample_Plate, ss$Sample_Well, sep='-')[match(matchtoinput, ind384)]
    readID.split=split(keep.reads, readID)
    foreach(i=1:length(readID.split)) %dopar% {
        file1=paste(out.base.dir, 'fastq/', names(readID.split)[i] ,'-', lane.identifier, '.fq.gz', sep='')
        #file1=paste('/data/eQTL/fastq/', names(readID.split)[i], '_R1.fq.gz', sep='')
        #file2=paste('/data/eQTL/fastq/', names(readID.split)[i], '_R2.fq.gz', sep='')
        file3=paste(out.base.dir, 'fastq/Indices/', names(readID.split)[i], ':I7.fq.gz', sep='')
        file4=paste(out.base.dir, 'fastq/Indices/', names(readID.split)[i], ':I5.fq.gz', sep='')

        writeFastq(fq.r1[readID.split[[i]]], file=file1, mode='a', full=FALSE, compress=TRUE)
        #writeFastq(qr2[readID.split[[i]]], file=file2, mode='a', full=FALSE, compress=TRUE)
        writeFastq(fq.i1[readID.split[[i]]], file=file3, mode='a', full=FALSE, compress=TRUE)
        writeFastq(fq.i2[readID.split[[i]]], file=file4, mode='a', full=FALSE, compress=TRUE)
    }
}
#-----------------------------------------------------------------------------------------------------------------------


# Version for qseq files
# For qseq files ... figure out what is happening to read names????
#qseq.files=list.files(paste(base.dir, lane.identifier, sep=''))
#qseq.sub=dcr(strsplit(qseq.files, '_'))
#qseq.grps=unique(sort(qseq.sub[,4]))
#lane.numbers=unique(sort(qseq.sub[,2]))
#
#for(l in lane.numbers ) {
#    for(g in qseq.grps) {
#            qr1=readQseq(paste(base.dir, lane.identifier, '/s_', l, '_1_', g, '_qseq.txt.gz', sep=''))
#            qi1=readQseq(paste(base.dir, lane.identifier, '/s_', l, '_2_', g, '_qseq.txt.gz', sep=''))
#            qi2=readQseq(paste(base.dir, lane.identifier, '/s_', l, '_3_', g, '_qseq.txt.gz', sep=''))
#            #qr2=readQseq(paste(base.dir, 'L', l, '/s_', l, '_4_', g, '_qseq.txt.gz', sep=''))
#
#            print('calculating distance to index1')
#            index1.dist.calc=srdistance(sread(qi1), n700.inds	)
#            d1mat=do.call('cbind', index1.dist.calc)
#            # find closest match for each read
#            d1min.ind=parRapply(cl, d1mat, which.min)
#            # find edit distance of closest match
#            d1min=d1mat[cbind(1:nrow(d1mat), d1min.ind)]
#
#            #lop off the last base
#            s7=DNAStringSet(sread(qi2),start=1, end=7)
#            print('calculating distance to index2')
#            # calculate distance matrix of reads indices (read2 is N500 read)
#            index2.dist.calc=srdistance(s7, n500.inds)
#            d2mat=do.call('cbind', index2.dist.calc) 
#            # find closest match for each read
#            d2min.ind=parRapply(cl, d2mat, which.min)
#            # find edit distance of closest match
#            d2min=d2mat[cbind(1:nrow(d2mat), d2min.ind)]
#
#            keep.reads=which(d1min<2 & d2min<2)
#            #print(paste( (length(keep.reads)/nbuffer)*100 ,  ' % of reads retained'))
#            d1.split.ind=d1min.ind[keep.reads]
#            d2.split.ind=d2min.ind[keep.reads]
#            matchtoinput=paste(d1.split.ind, d2.split.ind, sep='_')
#            #readID=ss$Sample_Name[match(matchtoinput, ind384)]
#            readID=paste(ss$Sample_Name, ss$Sample_Plate, ss$Sample_Well, sep='-')[match(matchtoinput, ind384)]
#            readID.split=split(keep.reads, readID)
#            foreach(i=1:length(readID.split)) %dopar% {
#                file1=paste('/data/eQTL/fastq/', names(readID.split)[i] ,'-', lane.identifier, '.fq.gz', sep='')
#                #file1=paste('/data/eQTL/fastq/', names(readID.split)[i], '_R1.fq.gz', sep='')
#                #file2=paste('/data/eQTL/fastq/', names(readID.split)[i], '_R2.fq.gz', sep='')
#                file3=paste('/data/eQTL/fastq/Indices/', names(readID.split)[i], '_I7.fq.gz', sep='')
#                file4=paste('/data/eQTL/fastq/Indices/', names(readID.split)[i], '_I5.fq.gz', sep='')
#
#                writeFastq(qr1[readID.split[[i]]], file=file1, mode='a', full=FALSE, compress=TRUE)
#                #writeFastq(qr2[readID.split[[i]]], file=file2, mode='a', full=FALSE, compress=TRUE)
#                writeFastq(qi1[readID.split[[i]]], file=file3, mode='a', full=FALSE, compress=TRUE)
#                writeFastq(qi2[readID.split[[i]]], file=file4, mode='a', full=FALSE, compress=TRUE)
#            }
#    }
#}



