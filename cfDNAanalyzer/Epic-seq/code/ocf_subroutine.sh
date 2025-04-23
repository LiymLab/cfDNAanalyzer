#!/bin/bash
# this function is written to extract the information regarding the lengths of the fragments for the "OCF" analysis 
# This uses the Python function to go after TSS coordinates U & D Signals
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
while read chrin tss gene category strand tss_id # genename category
do
  # Zhou Junpeng 20241105 modified
  out2="${outputcurr}/fracs.size"
  start=$(($tss-500))
  end=$(($tss+500))
  samtools view $file ${chrin}:${start}-${end} -F3084 -q $quality | awk '{if ($2=="99"||$2=="163"||$2=="97"||$2=="161" && $7=="=") print $2"\t"$4"\t"$9}' > "${outputcurr}/temppP.txt"
  printf "$chrin\t$tss\t$gene\t$strand\t" >> $out
  # out2="${outputcurr}/fracs.size"
  # formatted_tss=$(printf "%.0f" $tss)
  # start=$(($formatted_tss-500))
  # end=$(($formatted_tss+500))
  # formatted_start=$(printf "%.0f" $start)
  # formatted_end=$(printf "%.0f" $end)
  # samtools view $file ${chrin}:${formatted_start}-${formatted_end} -F3084 -q $quality | awk '{if ($2=="99"||$2=="163"||$2=="97"||$2=="161" && $7=="=") print $2"\t"$4"\t"$9}' > "${outputcurr}/temppP.txt"
  # printf "$chrin\t$formatted_tss\t$gene\t$strand\t" >> $out
  # end
  
  rt=$(python $path2dep/len2ocf.py "${outputcurr}/temppP.txt" $tss)
  echo -e $rt >> $out
done <$selector.noheader
