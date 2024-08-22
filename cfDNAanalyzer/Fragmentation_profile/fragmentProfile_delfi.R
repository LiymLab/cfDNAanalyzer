#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
usage=function(){
  cat('Description:This script was revised from https://github.com/cancer-genomics/delfi_scripts. It calculate GC-corrected fragmentation profile in 100kb bins of hg19 genome.\n',file=stderr())
  cat('Usage: fragmentProfile_delfi.R -i=input.bam -o=output_bins_100kb.tsv\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe input file in bam format (The corresponding index file should in the same directory).\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput tsv file\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

if(length(args)==0 || args[1]=="-h"){
  usage()
}else{
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the input file -i')
      }else{
        inFile=arg.split[2]
      }
    }else if(grepl('^-o=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the input file -i')
      }else{
        outFile=arg.split[2]
      }
    }
  }
}
# 2024.6.3 Zhou Junpeng modified
bedpath <- args[3]
reference <- args[4]
script_dir <- args[5]
# End
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(biovizBase)
library(httr)
library(tidyverse)
library(RCurl)

###### Functions ######
gc.correct <- function(coverage, bias) {
  i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend <- loess(coverage ~ bias)
  coverage.model <- loess(predict(coverage.trend, i) ~ i)
  coverage.pred <- predict(coverage.model, bias)
  coverage.corrected <- coverage - coverage.pred + median(coverage)
}
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# 2024.6.10 Zhou Junpeng modified
if (reference == "hg19") {
  library(BSgenome.Hsapiens.UCSC.hg19)
  filters_hg19_path = paste0(script_dir,"/Fragmentation_profile/hg19/filters.hg19.rda")
  gaps_hg19_path = paste0(script_dir,"/Fragmentation_profile/hg19/gaps.hg19.rda")
  load(filters_hg19_path)
  load(gaps_hg19_path)
} else if (reference == "hg38") {
  library(BSgenome.Hsapiens.UCSC.hg38)
  filters_hg38_path = paste0(script_dir,"/Fragmentation_profile/hg38/filters.hg19.rda")
  gaps_hg38_path = paste0(script_dir,"/Fragmentation_profile/hg38/gaps.hg19.rda")
  load(filters_hg38_path)
  load(gaps_hg38_path)
}
# end

###### Get 100-kb bins based on A/B compartments data ######
### A/B compartments of of lymphoblastoid cells in 100b bins from Hi-C data https://github.com/Jfortin1/HiC_AB_Compartments
# ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
# 2024.6.3 Zhou Junpeng modified
AB <- read.table(bedpath, header=TRUE)
# End
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)
chromosomes <- GRanges(paste0("chr", 1:22), IRanges(0, seqlengths(Hsapiens)[1:22]))
tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]
arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
### Chromosome 13 (as well as chromosomes 14, 15, 21 and 22) is an acrocentric chromosome. Short arms of acrocentric chromosomes do not contain any genes. 
### All genes are located in the long arm.
# 2024.6.10 Zhou Junpeng modified
if (reference == "hg19") {
arms <- arms[-c(25,27,29,41,43)]
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
} else if (reference == "hg38") {
armlevels <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
}
#end
arms$arm <- armlevels
AB <- AB[-queryHits(findOverlaps(AB, gaps.hg19))]
AB <- AB[queryHits(findOverlaps(AB, arms))]
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]
seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)

###### Read fragments from bam ###### only apply to small file due to high memory demand
param <- ScanBamParam(flag = scanBamFlag(isSecondaryAlignment=FALSE, isUnmappedQuery=FALSE), mapqFilter=30)
galp <- readGAlignmentPairs(inFile, param = param)
frags <- granges(keepSeqlevels(galp, paste0("chr", c(1:22,"X","Y")), pruning.mode="coarse"), on.discordant.seqnames="drop")
w.all <- width(frags)
q.all <- quantile(w.all, c(0.001, 0.999))
frags <- frags[which(w.all > q.all[1] & w.all < q.all[2])]
frags <- trim(frags)
seqlengths(frags) <- seqlengths(Hsapiens)[names(seqlengths(frags))]
seqlengths(AB) <- seqlengths(Hsapiens)[names(seqlengths(AB))]
gcs <- GCcontent(Hsapiens, unstrand(frags))
frags$gc <- gcs

###### Count fragments in each 100kb bin and do GC-correction ######
# frags <- frags[-queryHits(findOverlaps(frags, filters.hg19))]
# 2024.6.10 Zhou Junpeng modified
frags <- frags[-queryHits(findOverlaps(frags, filters.hg19))]
#end
w.all <- width(frags)
frags <- frags[which(w.all >= 100 & w.all <= 220)]
w <- width(frags)
frag.list <- split(frags, w)
counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
### Remove those no coverage for short or long fragments
uniqWidth <- as.numeric(colnames(counts))
short <- rowSums(counts[,which(uniqWidth<=150)]) ### 100-150 bp
long <- rowSums(counts[,which(uniqWidth>150)]) ### 151-220
noCov = unique(c(which(short==0),which(long==0)))
counts=counts[-noCov,]
AB=AB[-noCov]
short <- rowSums(counts[,which(uniqWidth<=150)])
long <- rowSums(counts[,which(uniqWidth>150)])
ratio <- short/long

olaps <- findOverlaps(frags, AB)
bin.list <- split(frags[queryHits(olaps)], subjectHits(olaps))
bingc <- rep(NA, nrow(counts))
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))

short.corrected=gc.correct(short, bingc)
long.corrected=gc.correct(long, bingc)
nfrags.corrected=gc.correct(short+long, bingc)
ratio.corrected=gc.correct(ratio, bingc)

##### Output results ######
AB$short <- short
AB$long <- long
AB$ratio <- short/long
AB$nfrags <- short+long
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected

AB$mode <- Mode(w)
AB$mean <- round(mean(w), 2)
AB$median <- median(w)
AB$quantile.25 <- quantile(w, 0.25)
AB$quantile.75 <- quantile(w, 0.75)
AB$frag.gc <- bingc
write.table(data.frame(AB), file=outFile, quote=F, sep="\t", row.names=F, col.names=T)
