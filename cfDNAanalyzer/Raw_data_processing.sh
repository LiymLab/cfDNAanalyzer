fastqc  -o ./6 ./6/*trimmed_cut.fq.gz

trim_galore --quality 20 --length 20 -o ./9 ./9/ENCFF937AFR.fastq.gz 

# mapping
#bwa
bwa index ref.fa
bwa mem -t 10 -M ../hg38.fa ./SRR3819940.sralite.1_1.fastq.gz ./SRR3819940.sralite.1_2.fastq.gz | samtools view -S -bF 4 - > ./SRR3819940_bwa_hg38.bam

#bowtie2
bowtie2 -p 15 -x ../bowtie2_index/bowtie2_index -1 SRR3819940.sralite.1_1.fastq.gz -2 SRR3819940.sralite.1_2.fastq.gz | samtools sort -O bam -@ 15 -o - > ./SRR3819940_bowtie2_hg19.bam

# filter bam
samtools view -h -@ 6 -b -q 10 -F 1796 ENCFF903GQY.bam > ENCFF903GQY_filtered.bam 