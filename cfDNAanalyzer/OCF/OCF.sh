#!/bin/bash

##
## Author : Ahfyth
## Contact: sunkun@cuhk.edu.hk
## Date   : Jan 18, 2019
##
## This program is designed to analyze the plasma DNA fragmentation pattern in tissue-specific
## open chromatin regions for infering the tissue origin of cfDNA, which information could be
## valuable in predicting the tumor origin after a positive cancer testing
##
## Note: this program only works for Paired-end sequencing data
## OCF: Orientation-aware CfDNA Fragmentation
##
## Please cite the following paper if you use this program in your work:
## Sun et al. Orientation-aware plasma cell-free DNA fragmentation analysis in open chromatin
## regions informs tissue of origin. Genome Res 2019.
##
## Revise by Yumei Li @ 11/15/2022 for any regions

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

$prg/bedtools intersect -a $inputRegion -b $data -wo > $outDir/$code.ol.$filename

#perl $prg/sync.pl     $code.ol.$filename $filename.vs.$code &
perl $prg/sync.end_custom.pl $outDir/$code.ol.$filename $outDir/$filename.vs.$code 
#R --slave --args $filename $filename <$prg/plot.sync.end.R

> $outDir/${code}.${filename}.OCF
ls $outDir/$filename.vs.${code}*.sync.end | while read file;do
	perl $prg/quant.pl $file >> $outDir/${code}.${filename}.OCF
done;
rm $outDir/$filename.vs.${code}*.sync.end



