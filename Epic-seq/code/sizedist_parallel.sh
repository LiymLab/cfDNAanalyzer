# This uses the Python function to have the bin counts
bamdir=$1
selector=$2
outdir=$3
quality=$4
path2dep=$5
rm $outdir/bamfiles.txt
for i in $bamdir/*.bam; do
echo $i >> $outdir/bamfiles.txt
done
echo $selector
parallel --results logfiles --jobs 50 --xapply bash $path2dep/sizedist_subroutine.sh \
:::: $outdir/bamfiles.txt \
::: $selector \
::: $outdir \
::: $quality \
::: $path2dep
