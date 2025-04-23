# This uses the Python function to have the bin counts
bamdir=$1
selector=$2
outdir=$3
mapq=$4
path2dep=$5
cat $selector | awk '{if (NR>1) print}' > $selector.noheader
rm $outdir/bamfiles.txt
for i in $bamdir/*.bam; do
echo $i >> $outdir/bamfiles.txt
done
echo $selector
parallel --results logfiles --jobs 20 --xapply bash $path2dep/ocf_subroutine.sh \
:::: $outdir/bamfiles.txt \
::: $selector \
::: $outdir \
::: $mapq \
::: $path2dep
