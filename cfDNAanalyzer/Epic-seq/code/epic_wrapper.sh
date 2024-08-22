#!/bin/bash
# This is the wrapper for the EPIC-Seq Technique (by Alizadeh and Diehn laboratories at Stanford University)
# The aim of this wrapper is to combine multiple steps which lead to full featurization of either Whole-Genome or Targeted cell-free DNA samples
# This wrapper calls a collection of R scripts and one single Python script. One major dependency of the wrapper is 'samtools' for fragmentation analysis
bamdir=$1
tssfile=$2 # This is the design file-- the TSSs of interest and other columns
outdir=$3 
selector=$4 # This is the actual coordinates in the designed selector  
mode=$5 #[either wgs or epicseq]
groupsize=$6 # If we want to group **GENES** which is used for training the model for gene expression prediction
groupmode=$7 # wholeblood / pbmc / etc.
mapq=$8 # mapping quality for the samtools view

skiphist=$9 # If 1 the histograms for the TSS genes will not be made
skiphistocf=${10} # If 1 histograms for the OCF will not be made
skipentropy=${11} # If 1 the entropy is not calculated (it's assumed that the entropy file does already exist)
skipdepth=${12} # If 1 the depth is not caluclated (it's assumed that the depth file does already exist)
skipocf=${13} # If 1 the OCF is not calculated (it's assumed that the ocf file does already exist)

# 20240711 zjp modified
epic_dir=${14}
# end 

ocfout="$outdir/ocfoutputs/"
mkdir -p $ocfout
scriptpath=$epic_dir

fragpath="$scriptpath/epicseq_histextract.R"
ocfpath="$scriptpath/ocf_parallel.sh"
convertpath="$scriptpath/epicseq_conv2r.R"
matrixgenpath="$scriptpath/epicseq_pathgen.R"
featurization="$scriptpath/epicseq_featurization.R"
calcoverlap="$scriptpath/epicseq_overlapsel.R"

tv=1
# NDR module:
#Rscript $calcoverlap --tssinfo $tssfile --selector $selector --outdir $outdir --path2dep $scriptpath
if [ "$skiphist" -ne 1 ]
then
echo "working on selector overlap calculations..."
# 20240711 zjp modified
# Rscript $calcoverlap --tsspath $tssfile --selector $selector --outdir $outdir --path2dep $scriptpath
# Rscript $fragpath --bamdir $bamdir --outdir $outdir  --tsspath $tssfile --mapq $mapq --path2dep $scriptpath
Rscript $calcoverlap --tsspath $tssfile --selector $selector --outdir $outdir --path2dep $scriptpath
Rscript $fragpath --bamdir $bamdir --outdir $outdir  --tsspath $tssfile --mapq $mapq --path2dep $scriptpath
# end
fi
# OCF module:
echo "************"
echo "$skiphistocf"
if [ "$skiphistocf" -ne "$tv" ]
then
$ocfpath $bamdir $tssfile $ocfout $mapq $scriptpath
fi

# 20240711 zjp modified
# Rscript $convertpath --tssinfo $tssfile --outdir $outdir --mode $mode --path2dep $scriptpath
Rscript $convertpath --tssinfo $tssfile --outdir $outdir --mode $mode --path2dep $scriptpath
# end 
if [ "$skiphistocf" -ne "$tv" ] || [ "$skiphist" -ne 1 ]
then
# 20240711 zjp modified
# Rscript $convertpath --tssinfo $tssfile --outdir $outdir --mode $mode --path2dep $scriptpath  # Here we convert all the above analysis into histogram formats + OCF format
Rscript $convertpath --tssinfo $tssfile --outdir $outdir --mode $mode --path2dep $scriptpath  # Here we convert all the above analysis into histogram formats + OCF format
# end
fi
########### Summarize data to create features (featurization)

newdirocf="$outdir/ocffiles/"
mkdir $newdirocf
cp $ocfout/*txt $newdirocf

### Second, summarize the entropy/depth/OCFmeasures
echo "Working on generating matrices"
echo $skipocf
# 20240711 zjp modified
# Rscript $matrixgenpath --mode $mode --newdirocf $newdirocf --outdir $outdir --groupsize $groupsize --tssinfopath $tssfile --sortmode $groupmode --skipdepth $skipdepth --skipentropy $skipentropy --skipocf $skipocf --path2dep $scriptpath
Rscript $matrixgenpath --mode $mode --newdirocf $newdirocf --outdir $outdir --groupsize $groupsize --tssinfopath $tssfile --sortmode $groupmode --skipdepth $skipdepth --skipentropy $skipentropy --skipocf $skipocf --path2dep $scriptpath
# end 
### Third, featurization step:
# 20240711 zjp modified
# Rscript $featurization --outdir $outdir --tssinfo $tssfile --groupsize $groupsize
Rscript $featurization --outdir $outdir --tssinfo $tssfile --groupsize $groupsize
#  end 