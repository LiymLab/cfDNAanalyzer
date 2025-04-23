#!/bin/bash

##
## Author : Ahfyth
## Contact: sunkun@cuhk.edu.hk
## Date   : Jan 18, 2019
##
## This program is designed to analyze the plasma DNA fragmentation pattern 

## Note: this program only works for Paired-end sequencing data
## OCF: Orientation-aware CfDNA Fragmentation
##
## Please cite the following paper if you use this program in your work:
## Sun et al. Orientation-aware plasma cell-free DNA fragmentation analysis in open chromatin
## regions informs tissue of origin. Genome Res 2019.
##
# Revise by Junpeng Zhou @ 11/24/2024 for any regions for one bed input

set -o nounset
set -o errexit

if [ $1 == "-h" ];then
	echo "Please set four parameters: inputRegion.bed4 data prefixForOutput outDir"
	exit 2
fi

prg=`dirname $0`
inputRegion=$1
data=$2
code=$3
outDir=$4

echo Processing: $data

filename=$(basename "$data" .bam)

$prg/bedtools intersect -a $inputRegion -b $data -wo -sorted > $outDir/$code.ol.$filename

#perl $prg/sync.pl     $code.ol.$filename $filename.vs.$code &
perl $prg/sync.end_custom.pl $outDir/$code.ol.$filename $outDir/$filename.vs.$code 
#R --slave --args $filename $filename <$prg/plot.sync.end.R

echo -n > $outDir/${code}.${filename}.OCF
find $outDir -maxdepth 1 -name "$filename.vs.${code}*.sync.end" | while read -r file; do
    perl $prg/quant.pl "$file" >> $outDir/${code}.${filename}.OCF
done


for file in $outDir/$filename.vs.${code}*.sync.end; do
    rm "$file"
done



