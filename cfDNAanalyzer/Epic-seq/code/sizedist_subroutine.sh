# This uses the Python function to have the bin counts
file=$1
selector=$2
outdir=$3
quality=$4 
path2dep=$5
name=$(basename $file)
name=${name%.bam}
outputcurr="${outdir}/${name}"
mkdir $outputcurr
out="${outputcurr}/bin_out_$name.txt"
echo -n > $out
while read chrin start end genename category tss
do
out2="${outputcurr}/fracs.size"
len=$(($end-$start))
samtools view $file ${chrin}:${start}-${end}  -F3084 -q $quality | awk '{if (($2==99 || $2==97 || $2==147 || $2==145 || $2==163 || $2==161 || $2==83 || $2==81) && $7=="="){if ($9>0){print $1, $9}else{print $1, -$9}}}'| sort | uniq -c | awk '{if ($1==2) print $3}'  > "${outputcurr}/temppP.txt"
echo $chrin
chrnew=${chrin#chr}
samtools view $file ${chrnew}:${start}-${end}  -F3084 -q $quality | awk '{if (($2==99 || $2==97 || $2==147 || $2==145 || $2==163 || $2==161 || $2==83 || $2==81) && $7=="="){if ($9>0){print $1, $9}else{print $1, -$9}}}'| sort| uniq -c | awk '{if ($1==2) print $3}'  >> "${outputcurr}/temppP.txt"
printf "$chrin\t$start\t$end\t$genename\t" >> $out
echo -e "0" >> "${outputcurr}/temppP.txt"
rt=$(python $path2dep/len2dist.py "${outputcurr}/temppP.txt")
echo -e $rt >> $out
done < $selector
