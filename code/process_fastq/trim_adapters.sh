#!/bin/bash
# This script only does adapter trimming and conversion from Illumin 1.3+ to Phred33

genome='/data/eQTL/reference/sacCer3.fasta'
cputhreads=10
root_dir='/data/eQTL/'

# sample prefixes (exclude file extension)
#snames=(`cat "${root_dir}all.sample.names" `)
#snames=(`cat "${root_dir}sample.names_eQTL1" `)
#snames=(`cat "${root_dir}sample.names_eQTL2" `)
#snames=(`cat "${root_dir}sample.names_eQTL4-10" `)
#snames=(`cat "${root_dir}sample.names_eQTL_3_12_13" `)
#snames=(`cat "${root_dir}sample.names_eQTL_repeat1" `)
#snames=(`cat "${root_dir}sample.names_eQTL_13" `)
#snames=(`cat "${root_dir}sample.names_eQTL_repeat12_13" `)
#snames=(`cat "${root_dir}sample.names_eQTL_repeat_3" `)
snames=(`cat "/data/eQTL/reference/all.samples.101215" `)

#convert to phred33, trim, and align with bwa mem
for sample in ${snames[@]}
do
    #adapter trimming
    java -jar /home/jbloom/Local/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads ${cputhreads} \
    ${root_dir}/fastq/${sample}.fq.gz  ${root_dir}/fastq_trimmed/${sample}.fq.gz \
    ILLUMINACLIP:${root_dir}reference/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 TOPHRED33
done


#for sample in ${snames[@]}
#do
 # run kallisto quant with default settings on untrimmed data
    #kallisto quant --single -l 150 -i ${root_dir}reference/transcripts.idx -o ${root_dir}kallisto_out/$sample ${root_dir}fastq/${sample}.fq.gz
   # kallisto quant --single -l 150 -i ${root_dir}reference/transcripts.idx -o ${root_dir}kallisto_out_trimmed/$sample ${root_dir}fastq_trimmed/${sample}.fq.gz
   # alignment with BWA
   # bwa mem -t $cputhreads ${genome} ${root_dir}fastq_trimmed/${sample}.fq.gz | samtools view -huS - | samtools sort - ${root_dir}bwa_bam/${sample}
   # samtools index ${root_dir}bwa_bam/${sample}.bam
#done 

#install emboss 
#revseq Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa -reverse -complement -outseq  Saccharomyces_cerevisiae.R64-1-1.cdna.all.rev.fa
#kallisto index -i rev_transcripts.idx Saccharomyces_cerevisiae.R64-1-1.cdna.all.rev.fa.gz
#cat Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa Saccharomyces_cerevisiae.R64-1-1.cdna.all.rev.fa > plus_minus_cat.R64.fa
#kallisto index -i plus_minus_transcripts.idx plus_minus_cat.R64.fa.gz
#kallisto index -i test.idx test.fa.gz 
# ---------------------------------------------------------------------------------

#for sample in ${snames[@]}
#do
#  kallisto quant --single -l 150 -i ${root_dir}reference/rev_transcripts.idx -o ${root_dir}kallisto_out_rev/$sample ${root_dir}fastq/${sample}.fq.gz
#done
#for sample in ${snames[@]}
#do
#  kallisto quant --single -l 150 -i ${root_dir}reference/plus_minus_transcripts.idx -o ${root_dir}kallisto_out_plus_minus/$sample ${root_dir}fastq/${sample}.fq.gz
#done


##run kallisto quant with default settings on trimmed data
#for sample in ${snames[@]}
#do
#  kallisto quant --single -l 150 -i ${root_dir}reference/transcripts.idx -o ${root_dir}kallisto_out_trimmed/$sample ${root_dir}fastq_trimmed/${sample}.fq.gz
#done
#
##ehttps://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
#STAR --genomeDir /data/eQTL/reference \
#     --genomeFastaFiles /data/eQTL/reference/sacCer3.fasta \
#     --runMode genomeGenerate \
#     --runThreadN 11 --sjdbGTFfile  /data/eQTL/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.gtf
#
## combine this with run merged kallisto script and input multiple files
#for sample in ${snames[@]}
#do
#mkdir ${root_dir}star_out/${sample}
#STAR --genomeDir /data/eQTL/reference/STAR \
#     --runThreadN ${cputhreads} \
#     --readFilesIn ${root_dir}fastq_trimmed/${sample}.fq.gz \
#    --readFilesCommand zcat \
#    --outSAMtype BAM SortedByCoordinate \
#    --quantMode TranscriptomeSAM GeneCounts \
#    --outFileNamePrefix ${root_dir}star_out/${sample}/ \
#    --outReadsUnmapped Fastx \
#    --outWigType bedGraph \
#    --outWigNorm None
#done
#
## directory names reflect correct forward or reverse strand mappings
#for sample in ${snames[@]}
#    do
#    #majority of reads are accounted for here
#    express --r-stranded /data/eQTL/reference/R64-1.1.cdna.ncrna.fa  ${root_dir}star_out/${sample}/Aligned.toTranscriptome.out.bam -o  ${root_dir}star_out/${sample}/F/
#    #the rest 
#    express --f-stranded /data/eQTL/reference/R64-1.1.cdna.ncrna.fa  ${root_dir}star_out/${sample}/Aligned.toTranscriptome.out.bam -o   ${root_dir}star_out/${sample}/R/
#done
#
#
#
#
#
########TESTING CUFFLINKS, TOPHAT, other ...#######################################################################################################3
#
#
#
#gffread /data/eQTL/reference/Saccharomyces_cerevisiae.R64-1-1.81.gff3 -g /data/eQTL/reference/sacCer3.fasta -w transcripts_from_gff3.fa
## consider workflow with bowtie or STAR
## and express
##http://homer.salk.edu/homer/basicTutorial/mapping.html
##http://bio.math.berkeley.edu/eXpress/manual.html
## kallisto
##http://arxiv.org/pdf/1505.02710v2.pdf
#
## express workflow
##http://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html
##http://piquant.readthedocs.org/en/latest/quantifiers.html#express
#
## BWA alignment (standalone)
##for sample in  ${snames[@]}
##    do 
##    echo $sample
##    # note, this is for PE reads ... adjust for SE reads
##    bwa mem -t $cputhreads ${genome} ${root_dir}fastq_trimmed/${sample}.fq.gz | samtools view -huS - | samtools sort - ${root_dir}bam/${sample}
##    samtools index ${root_dir}bam/${sample}.bam
##done
#
## testing Cufflinks
## note I modified gtf chromosome names to match reference to circumvent an annoying error
##note -p for threading
#tophat2 -G /data/eQTL/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.gtf \
#    /data/eQTL/reference/Saccharomyces_cerevisiae_UCSC_sacCer3/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome \
#    A11_91-A11-H7-BYxRM_eQTL_08-H9-SxaQSEQsVB101L5.fq.gz \
#    -o tophat_out
#cufflinks -o A11_91 accepted_hits.bam
#
#
## improved
#tophat2 -o /data/eQTL/tophat_out1/ \
#    -G /data/eQTL/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.gtf \
#    /data/eQTL/reference/Saccharomyces_cerevisiae_UCSC_sacCer3/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome \
#    /data/eQTL/A10_48-A10-D12-BYxRM_eQTL_01-G9-SxaQSEQsVB100L6.fq.gz 
#cufflinks -o /data/eQTL/tophat_out1/A10_48 /data/eQTL/tophat_out1/accepted_hits.bam
#
#cuffmerge -g /data/eQTL/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.gtf \
#    -s /data/eQTL/reference/Saccharomyces_cerevisiae_UCSC_sacCer3/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome.fa \
#    /data/eQTL/cuffassemblies.txt
#
##cuffdiff -o /data/eQTL/diff_out \
##    -b /data/eQTL/reference/Saccharomyces_cerevisiae_UCSC_sacCer3/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome.fa \
##    -u /data/eQTL/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.gtf \
##     /data/eQTL/tophat_out/accepted_hits.bam \
##     /data/eQTL/tophat_out1/accepted_hits.bam
