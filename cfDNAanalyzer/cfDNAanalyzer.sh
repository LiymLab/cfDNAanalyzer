#!/bin/bash

usage(){
  cat << EOF
Description:
  cfDNAanalyzer (cell-free DNA sequencing data analyzer) is a toolkit for cfDNA whole-genome sequencing data analysis which includes two main parts:
  (1) the extraction of genomic and fragmentatomic features at whole-genome or genomic-region levels;
  (2) the processing of extracted features and the building of machine learning models for disease detection and classification.

Usage：
   bash cfDNAanalyzer.sh -I <InputFile> -o <OutputDirectory> -F <Features> [Options]
  
Options: 
-- General options
  -I  FILE      A text file containing all input BAM files with one BAM file per line. 
                BAM files generated using both Bowtie2 and BWA are accepted.  
  -o  DIR       Output directory for all the results. Default: [./]
  -F  STR       Features to extract, including CNA, NOF, WPS, EM, EMR, FP, FPR, NP, OCF, PFE, and TSSC. 
                Features should be set as a string separated by comma, e.g., CNA,NOF. 
                Default: All available features will be extracted.
                Note: The following features are specifically designed for paired-end sequencing data: FP, FPR, NP, PFE, and OCF.
  -g  STR       Genome version of input BAM files (hg19/hg38). Default: [hg38] 
  -b  FILE      A BED3 file specifying the regions to extract features.
                The file should contain three TAB-delimited columns: chromosome start end.              

-- Options specific for Copy Number Alterations (CNA)
  -B  INT       Bin size in kilobases (10, 50, 500, or 1000). Default: [1000]
  --CNA  STR    Additional parameter setting for software ichorCNA. 
                The full parameter list is available by running Rscript cfDNAanalyzer/ichorCNA/ichorCNA/scripts/runIchorCNA.R --help. [optional]

-- Options specific for Nucleosome Occupancy and Fuzziness (NOF)
  --NOF  STR    Additional parameter setting for software DANPOS2. 
                The full parameter list is available by running python cfDNAanalyzer/DANPOS3/danpos.py dpos -h. [optional]

-- Options specific for Windowed Protection Score (WPS)
  -x  INT       Min fragment length used for long fragments WPS calculation. Default: [120]
  -X  INT       Max fragment length used for long fragments WPS calculation. Default: [180]
  -w  INT       Window size for long fragments WPS calculation. Default: [120]
  -m  INT       Min fragment length used for short fragments WPS calculation. Default: [35]
  -M  INT       Max fragment length used for short fragments WPS calculation. Default: [80]
  -W  INT       Window size for short fragments WPS calculation. Default: [16]

-- Options specific for End Motif frequency and diversity (EM)
  -f  FILE      Reference genome in FASTA format. 
                For example, hg38 reference genome can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz.

-- Options specific for Nucleosome Profile (NP)
  -l  DIR       Directory containing a list of files with each file for a set of sites.
                The file must have at least two columns with the following column names: 
                "Chrom" for the chromosome name and "position" for the site positions. 
                If not provided, the 377 TF binding site lists from the referenced Nucleosome Profile paper will be used . 

-- Options for Promoter Fragmentation Entropy (PFE)
  -T  FILE      A TAB-delimited TSS sites gene information file without any header.
                The file must have six columns: (1) chromosome, (2) 1-base TSS coordinate, (3) Hugo symbol of the gene corresponding to the TSS,
                (4) Category (e.g., WGS), (5) gene transcript strand (+1/-1) and (6) a column for TSS ID (e.g., for genes with multiple TSS, this should be geneID_1, geneID_2, etc.)
                If not provided, 20139 genes in TSS sites from the referenced Promoter Fragmentation Entropy paper will be used.
  --PFE  STR    Addtional parameter setting for PFE analysis. 
                The full parameter list is available by running Rscript cfDNAanalyzer/Epic-seq/code/runEPIC.R -h.[optional]

-- Options for TSS Coverage (TSSC)
  -u  INT                    Number of base pairs upstream of TSS used for calculating TSSC. Default: [1000]
  -d  INT                    Number of base pairs downstream of TSS used for calculating TSSC. Default: [1000]
  -S  FILE                   A BED6 file specifying the coordinates of TSSs used for calculating TSSC.  
  --bamCoverage  STR         Additional parameter seting for the "bamCoverage" command of deeptools. [optional] 
  --multiBigwigSummary  STR  Additional parameter seting for the "multiBigwigSummary" command of deeptools. [optional]
EOF
    exit 0
}

if [ $# -eq 1 ] && [ $1 == -h ]; then usage; fi

script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"

inputBamList=""
outputDir=$(pwd)
Feature=CNA,EM,FP,NOF,NP,WPS,OCF,EMR,FPR,PFE,TSSC
genomeVer=hg38
fasta=""
sequencing=pair
bedfiles=""
clonalCN=7
binsize=1000
upstream=1000
downstream=1000
TSSfile="$script_dir/TSScoverage/test_hg19.bed"
Minlength_long=120
Maxlength_long=180
window_long=120
Minlength_short=35
Maxlength_short=80
window_short=16
sitespath=$script_dir/Griffin/Ref/sites
tssinfo=$script_dir/Epic-seq/code/priordata/sample_hg19.txt
addCNA=""
addNOF=""
addbamCoverage=""
addmultiBigwigSummary=""
addEpicSeq=""

# Use getopt to handle long options
TEMP=$(getopt -o hI:o:F:g:f:s:b:c:B:u:d:S:x:X:w:m:M:W:l:T: \
               --long help,CNA:,NOF:,bamCoverage:,multiBigwigSummary:,PFE: \
               -- "$@")

# Check whether getopt succeeded
if [ $? != 0 ]; then
    echo "Error in command line parsing" >&2
    exit 1
fi

# Reset position parameters
eval set -- "$TEMP"

# Process options
while true; do
    case "$1" in
        -h|--help)                usage; shift;;
        -I)                       inputBamList=$2; shift 2;;
        -o)                       outputDir=$2; shift 2;;
        -F)                       Feature=$2; shift 2;;
        -g)                       genomeVer=$2; shift 2;;
        -f)                       fasta=$2; shift 2;;
        -s)                       sequencing=$2; shift 2;;
        -b)                       bedfiles=$2; shift 2;;
        -c)                       clonalCN=$2; shift 2;;
        -B)                       binsize=$2; shift 2;;
        -u)                       upstream=$2; shift 2;;
        -d)                       downstream=$2; shift 2;;
        -S)                       TSSfile=$2; shift 2;;
        -x)                       Minlength_long=$2; shift 2;;
        -X)                       Maxlength_long=$2; shift 2;;
        -w)                       window_long=$2; shift 2;;
        -m)                       Minlength_short=$2; shift 2;;
        -M)                       Maxlength_short=$2; shift 2;;
        -W)                       window_short=$2; shift 2;;
        -l)                       sitespath=$2; shift 2;;
        -T)                       tssinfo=$2; shift 2;;
        --CNA)                    addCNA=$2; shift 2;;
        --NOF)                    addNOF=$2; shift 2;;
        --bamCoverage)            addbamCoverage=$2; shift 2;;
        --multiBigwigSummary)     addmultiBigwigSummary=$2; shift 2;;
        --PFE)                    addEpicSeq=$2; shift 2;;
        --)                       shift; break;;
        *)                        echo "Unexpected option: $1"; exit 1;;
    esac
done

# Main function
run_analysis() {
  # Ensure input file list is provided
  if [ -z "$inputBamList" ]; then
    echo "Error: Input bam file list not provided."
    exit 1
  fi

  # Create output directory if it doesn't exist
  mkdir -p "$outputDir"

  # Validate reference type
  if [[ "$genomeVer" != "hg19" ]] && [[ "$genomeVer" != "hg38" ]]; then
    echo "Error: incorrect type of reference fasta file"
    exit 1
  fi
  
  if [[ "$Feature" == *"NOF"* ]] || [[ "$Feature" == *"NP"* ]] || [[ "$Feature" == *"WPS"* ]] || [[ "$Feature" == *"OCF"* ]] || [[ "$Feature" == *"EMR"* ]] || [[ "$Feature" == *"FPR"* ]] && [ -z "$bedfiles" ]; then
    echo "Error: Region bed file is not provided."
    exit 1
  fi
  
  # Convert relative paths to absolute paths
  inputBamList=$(realpath "$inputBamList")
  outputDir=$(realpath "$outputDir")
  fasta=$(realpath "$fasta")
  bedfiles=$(realpath "$bedfiles")
  TSSfile=$(realpath "$TSSfile")
  sitespath=$(realpath "$sitespath")
    
    
  # Process each bam file
  for inputBam in $(cat "$inputBamList"); do
      if [ ! -f "$inputBam" ]; then
        echo "Error: Input bam file not found: $inputBam"
        exit 1
      fi

      samtools index -@ 5 "$inputBam"

      local filename=$(basename "$inputBam" .bam)
      local bam_output_dir="$outputDir/$filename"
      mkdir -p "$bam_output_dir"
      
      # Features for whole genome
      if [[ "$Feature" == *"CNA"* ]]; then calculate_cna "$inputBam" "$bam_output_dir"; fi
      if [[ "$Feature" == *"EM"* ]]; then calculate_em "$inputBam" "$bam_output_dir"; fi
      if [[ "$Feature" == *"FP"* ]]; then calculate_fp "$inputBam" "$bam_output_dir"; fi
      
      # Features for specific regions
      if [[ "$Feature" == *"NOF"* ]]; then calculate_nof "$inputBam" "$bam_output_dir"; fi
      if [[ "$Feature" == *"NP"* ]]; then calculate_np "$inputBam" "$bam_output_dir"; fi
      if [[ "$Feature" == *"WPS"* ]]; then calculate_wps "$inputBam" "$bam_output_dir"; fi
      if [[ "$Feature" == *"OCF"* ]]; then calculate_ocf "$inputBam" "$bam_output_dir"; fi
      if [[ "$Feature" == *"EMR"* ]]; then calculate_emr "$inputBam" "$bam_output_dir"; fi      
      if [[ "$Feature" == *"FPR"* ]]; then calculate_fpr "$inputBam" "$bam_output_dir"; fi
      
      # Features for transcription start sites 
      if [[ "$Feature" == *"PFE"* ]]; then calculate_pfe "$inputBam" "$bam_output_dir"; fi
      if [[ "$Feature" == *"TSSC"* ]]; then calculate_tssc "$inputBam" "$bam_output_dir"; fi
      echo "Analysis for file: $inputBam completed."
  done

  wait
}

# Calculate Copy Number Alterations
calculate_cna() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating CNA for file: '$inputBam'*****'
  if [[ "$binsize" != "10" ]] && [[ "$binsize" != "50" ]] && [[ "$binsize" != "500" ]] && [[ "$binsize" != "1000" ]]; then
    echo "Error: incorrect bin size"
    exit 1
  fi
  mkdir -p $bam_output_dir/CNA
  $script_dir/ichorCNA/hmmcopy_utils/bin/readCounter --window $((binsize * 1000)) --quality 20 --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $inputBam > $bam_output_dir/CNA/CNA.wig
  if [[ "$genomeVer" == "hg19" ]];then
    Rscript $script_dir/ichorCNA/ichorCNA/scripts/runIchorCNA.R --id sample --WIG $bam_output_dir/CNA/CNA.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN $clonalCN --gcWig $script_dir/ichorCNA/ichorCNA/inst/extdata/gc_hg19_${binsize}kb.wig --mapWig $script_dir/ichorCNA/ichorCNA/inst/extdata/map_hg19_${binsize}kb.wig --centromere $script_dir/ichorCNA/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $bam_output_dir/CNA $addCNA
  elif [[ "$genomeVer" == "hg38" ]];then
    Rscript $script_dir/ichorCNA/ichorCNA/scripts/runIchorCNA.R --id sample --WIG $bam_output_dir/CNA/CNA.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN $clonalCN --gcWig $script_dir/ichorCNA/ichorCNA/inst/extdata/gc_hg38_${binsize}kb.wig --mapWig $script_dir/ichorCNA/ichorCNA/inst/extdata/map_hg38_${binsize}kb.wig --centromere $script_dir/ichorCNA/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir $bam_output_dir/CNA $addCNA
  fi
  cut -f1,2,3,4 $bam_output_dir/CNA/sample.cna.seg > $bam_output_dir/CNA/CNA.txt
  rm -r $bam_output_dir/CNA/sample
  rm $bam_output_dir/CNA/sample.cna.seg
  rm $bam_output_dir/CNA/sample.correctedDepth.txt
  rm $bam_output_dir/CNA/sample.params.txt
  rm $bam_output_dir/CNA/sample.RData
  rm $bam_output_dir/CNA/sample.seg
  rm $bam_output_dir/CNA/sample.seg.txt
  rm $bam_output_dir/CNA/CNA.wig
}

# Calculate End motif frequency and diversity
calculate_em() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating EM for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/EM
  bedtools bamtobed -i $inputBam > $bam_output_dir/EM/${filename}.bed
  python $script_dir/End_motif_frequency/endMotifFreq_WG.py -b $bam_output_dir/EM/${filename}.bed -f $fasta -o $bam_output_dir/EM
  rm $bam_output_dir/EM/${filename}.bed
}

# Calculate Fragmentation Profile
calculate_fp() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating FP for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/FP
  Rscript $script_dir/Fragmentation_profile/fragmentProfile_delfi.R -i=$inputBam -o=$bam_output_dir/FP/FragmentationProfile.txt $script_dir/Fragmentation_profile/hic_compartments_100kb_ebv_2014.txt $genomeVer $script_dir
  cut -f1,2,3,9,10,11 $bam_output_dir/FP/FragmentationProfile.txt > $bam_output_dir/FP/Fragmentation_Profile.txt
  rm $bam_output_dir/FP/FragmentationProfile.txt
}


# Calculate Nucleosome Occupancy and Fuzziness
calculate_nof() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating NOF for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/NOF
  python $script_dir/DANPOS3/danpos.py dpos $inputBam -a 10 -o $bam_output_dir/NOF $addNOF
  if [[ "$genomeVer" == "hg19" ]];then
    $script_dir/TSScoverage/wigToBigWig -clip $bam_output_dir/NOF/pooled/*.wig $script_dir/DANPOS3/hg19.chrom.sizes $bam_output_dir/NOF/NOF.bw
  elif [[ "$genomeVer" == "hg38" ]];then
    $script_dir/TSScoverage/wigToBigWig -clip $bam_output_dir/NOF/pooled/*.wig $script_dir/DANPOS3/hg38.chrom.sizes $bam_output_dir/NOF/NOF.bw
  fi
  awk 'BEGIN {OFS="\t"} {print $0, "region_" NR}' $bedfiles > $bam_output_dir/NOF/regions_with_ids.bed
  $script_dir/TSScoverage/bigWigAverageOverBed $bam_output_dir/NOF/NOF.bw  $bam_output_dir/NOF/regions_with_ids.bed $bam_output_dir/NOF/occupancy.tsv -bedOut=$bam_output_dir/NOF/occupancy.bed
  sed -i '1i\chr\tstart\tend\tID\toccupancy' $bam_output_dir/NOF/occupancy.bed
  awk -v OFS="\t" 'NR>1{print $1,$4-1,$4,$5,$6}' $bam_output_dir/NOF/pooled/*.xls |grep -vE "Lambda_NEB|pUC19"| bedtools intersect -a $bedfiles -b stdin -wo > $bam_output_dir/NOF/intersect.tsv 
  awk '{sum[$1$2$3]+=$8; count[$1$2$3]++; total[$1$2$3]=$1" "$2" "$3} END {for (key in sum) print total[key], sum[key]/count[key]}' $bam_output_dir/NOF/intersect.tsv | sort -k1,1V -k2,2n -k3,3n | awk 'BEGIN {print "chr\tstart\tend\tmeanfuziness"} {print}' > $bam_output_dir/NOF/meanfuziness.tsv
  rm $bam_output_dir/NOF/intersect.tsv
  rm $bam_output_dir/NOF/NOF.bw
  rm $bam_output_dir/NOF/regions_with_ids.bed
  rm $bam_output_dir/NOF/occupancy.tsv
}

# Calculate Nucleosome Profile
calculate_np() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating NP for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/NP
  # Defining path ariables
  griffin_path="$script_dir/Griffin"
  result_dir="$griffin_path/snakemakes"
  gc_correction_dir="$result_dir/griffin_GC_and_mappability_correction"
  nucleosome_profiling_dir="$result_dir/griffin_nucleosome_profiling"
        
  # Preparing for GC and mappability corrections
  sed -e "s!Griffin_path!$griffin_path!g" -e "s!reference_path!$fasta!g" \
      -e "s!result_name!$gc_correction_dir/result_${filename}!g" "$gc_correction_dir/config/config.yaml" \
      > "$gc_correction_dir/config/config_${filename}.yaml"
        
  sed "s!input!$inputBam!g" "$gc_correction_dir/config/samples.yaml" \
      > "$gc_correction_dir/config/samples_${filename}.yaml"
        
  sed -e "s!config_samples!$gc_correction_dir/config/samples_${filename}.yaml!g" \
      -e "s!config_config!$gc_correction_dir/config/config_${filename}.yaml!g" \
      -e "s!config_cluster!$gc_correction_dir/config/cluster_slurm.yaml!g" \
      "$gc_correction_dir/griffin_GC_and_mappability_correction.snakefile" \
      > "$gc_correction_dir/griffin_GC_and_mappability_correction_${filename}.snakefile"
        
  snakemake -s "$gc_correction_dir/griffin_GC_and_mappability_correction_${filename}.snakefile" --cores 1
        
  rm "$gc_correction_dir/config/config_${filename}.yaml" "$gc_correction_dir/config/samples_${filename}.yaml"
  rm "$gc_correction_dir/griffin_GC_and_mappability_correction_${filename}.snakefile"
  mv "$gc_correction_dir/result_${filename}/samples.GC.yaml" "$nucleosome_profiling_dir/config/samples.GC_${filename}.yaml"
        
  # Nucleosome profiling analysis
  sed -e "s!Griffin_path!$griffin_path!g" -e "s!reference_path!$fasta!g" \
      -e "s!result_name!$nucleosome_profiling_dir/result_${filename}!g" \
      -e "s!tmp_name!$nucleosome_profiling_dir/tmp_${filename}!g" \
      -e "s!config_cluster!$nucleosome_profiling_dir/config/cluster_slurm.yaml!g" \
      -e "s!config_sites!$nucleosome_profiling_dir/config/sites_${filename}.yaml!g" \
      -e "s!config_samples.GC!$nucleosome_profiling_dir/config/samples.GC_${filename}.yaml!g" \
      "$nucleosome_profiling_dir/config/config.yaml" \
      > "$nucleosome_profiling_dir/config/config_${filename}.yaml"
        
  # Generate site list for nucleosome analysis
  touch "$nucleosome_profiling_dir/config/sites_${filename}.yaml"
  echo "site_lists:" > "$nucleosome_profiling_dir/config/sites_${filename}.yaml"
  find "$sitespath" -type f | while read file; do
      echo " $(basename "$file"): $file" >> "$nucleosome_profiling_dir/config/sites_${filename}.yaml"
  done
        
  sed -e "s!config_cluster!$nucleosome_profiling_dir/config/cluster_slurm.yaml!g" \
      -e "s!config_sites!$nucleosome_profiling_dir/config/sites_${filename}.yaml!g" \
      -e "s!config_samples.GC!$nucleosome_profiling_dir/config/samples.GC_${filename}.yaml!g" \
      -e "s!config_config!$nucleosome_profiling_dir/config/config_${filename}.yaml!g" \
      "$nucleosome_profiling_dir/griffin_nucleosome_profiling.snakefile" \
      > "$nucleosome_profiling_dir/griffin_nucleosome_profiling_${filename}.snakefile"
        
  snakemake -s "$nucleosome_profiling_dir/griffin_nucleosome_profiling_${filename}.snakefile" --cores 1 --unlock
  snakemake -s "$nucleosome_profiling_dir/griffin_nucleosome_profiling_${filename}.snakefile" --cores 1
        
  # Extract columns required for nucleosome analysis results
  header=$(head -n 1 "$nucleosome_profiling_dir/result_${filename}/sample_name_1/sample_name_1.GC_corrected.coverage.tsv")
  mean_coverage_col=$(echo "$header" | tr '\t' '\n' | grep -n -m 1 "mean_coverage" | cut -d: -f1)
  central_coverage_col=$(echo "$header" | tr '\t' '\n' | grep -n -m 1 "central_coverage" | cut -d: -f1)
  amplitude_col=$(echo "$header" | tr '\t' '\n' | grep -n -m 1 "amplitude" | cut -d: -f1)
  site_name_col=$(echo "$header" | tr '\t' '\n' | grep -n -m 1 "site_name" | cut -d: -f1)
        
  awk -v mean_col="$mean_coverage_col" -v central_col="$central_coverage_col" \
      -v amp_col="$amplitude_col" -v sites_col="$site_name_col" \
      'BEGIN {FS=OFS="\t"} NR==1 {print $(sites_col), $(mean_col), $(central_col), $(amp_col); next} {print $(sites_col), $(mean_col), $(central_col), $(amp_col)}' \
      "$nucleosome_profiling_dir/result_${filename}/sample_name_1/sample_name_1.GC_corrected.coverage.tsv" \
      > "$nucleosome_profiling_dir/result_${filename}/NucleosomeProfile.txt"
        
  rm -r "$gc_correction_dir/result_${filename}"
  rm "$nucleosome_profiling_dir/config/sites_${filename}.yaml" "$nucleosome_profiling_dir/config/samples.GC_${filename}.yaml" "$nucleosome_profiling_dir/config/config_${filename}.yaml"
  rm "$nucleosome_profiling_dir/griffin_nucleosome_profiling_${filename}.snakefile"
  mv "$nucleosome_profiling_dir/result_${filename}" "$bam_output_dir/NP"
  rm -r "$outputDir/$filename/NucleosomeProfile/sample_name_1"
}

# Calculate Windowed Protection Score
calculate_wps() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating WPS for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/WPS
  > $bam_output_dir/WPS/WPS.txt
  # Extract the first three columns of each line of the file (chr，start，end)
  while IFS=$'\t' read -r col1 col2 col3; do
    if [ "$col2" -ge 500 ]; then
      ex_start=$((col2 - 500))
    else
      ex_start=$col2
    fi
    ex_end=$((col3 + 500))
    $script_dir/WPS/calculate_wps.sh $inputBam $bam_output_dir/WPS $col1:$col2-$col3 $col1:$ex_start-$ex_end $Minlength_long $Maxlength_long $window_long $Minlength_short $Maxlength_short $window_short $script_dir
    gzip -df $bam_output_dir/WPS/*.gz
    # Calculate the average value in the wig file (long/short) for each region and write it to the result file
    long_average=$(awk '{sum += $1} END {if (NR > 1) print sum / (NR - 1)}' "$bam_output_dir/WPS/WPS_long.wig")
    short_average=$(awk '{sum += $1} END {if (NR > 1) print sum / (NR - 1)}' "$bam_output_dir/WPS/WPS_short.wig")
    echo -e "$col1\t$col2\t$col3\t$long_average\t$short_average" >> $bam_output_dir/WPS/WPS.txt
  done < $bedfiles
  sed -i '1i\chr\tstart\tend\tlong_WPS\tshort_WPS' $bam_output_dir/WPS/WPS.txt
  rm $bam_output_dir/WPS/*.wig
}

# Calculate Orientation-aware CfDNA Fragmentation
calculate_ocf() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating OCF for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/OCF
  awk 'BEGIN {OFS="\t"} {print $0, "region_" NR}' $bedfiles > $bam_output_dir/OCF/regions_with_ids.bed
  $script_dir/OCF/OCF.sh $bam_output_dir/OCF/regions_with_ids.bed $inputBam sample $bam_output_dir/OCF
  awk -v OFS="\t" '{split($1,a,".");print a[4],$2}' $bam_output_dir/OCF/sample.${filename}.OCF > $bam_output_dir/OCF/unsorted.OCF
  awk '{print substr($1, 8), $0}' $bam_output_dir/OCF/unsorted.OCF | sort -n | cut -d' ' -f2- > $bam_output_dir/OCF/sorted.OCF
  join -1 4 $bam_output_dir/OCF/regions_with_ids.bed -2 1 $bam_output_dir/OCF/sorted.OCF > $bam_output_dir/OCF/OCF_id.txt
  cut -d ' ' -f2,3,4,5 $bam_output_dir/OCF/OCF_id.txt > $bam_output_dir/OCF/OCF.txt
  sed -i '1i\chr\tstart\tend\tOCF' $bam_output_dir/OCF/OCF.txt
  rm $bam_output_dir/OCF/OCF_id.txt
  rm $bam_output_dir/OCF/regions_with_ids.bed
  rm $bam_output_dir/OCF/sample.ol.${filename}
  rm $bam_output_dir/OCF/sample.${filename}.OCF
  rm $bam_output_dir/OCF/unsorted.OCF
  rm $bam_output_dir/OCF/sorted.OCF
}

# Calculate End Motif frequency and diversity for Regions
calculate_emr() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating EMR for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/EMR
  python $script_dir/End_motif_frequency/endMotifFreq_region.py -b $bedfiles -f $fasta -o $bam_output_dir/EMR
}

# Calculate Fragmentation Profile for Regions
calculate_fpr() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating FPR for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/FPR
  awk 'BEGIN {print "chr\tstart\tend"} {print}' $bedfiles > $bam_output_dir/FPR/regions_with_header.bed
  Rscript $script_dir/Fragmentation_profile/fragmentProfile_delfi.R -i=$inputBam -o=$bam_output_dir/FPR/FragmentationProfile.txt $bam_output_dir/FPR/regions_with_header.bed $genomeVer $script_dir
  cut -f1,2,3,7,8,9 $bam_output_dir/FPR/FragmentationProfile.txt > $bam_output_dir/FPR/Fragmentation_Profile_regions.txt
  rm $bam_output_dir/FPR/FragmentationProfile.txt
  rm $bam_output_dir/FPR/regions_with_header.bed
}

# Calculate Promoter Fragmentation Entropy
calculate_pfe() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating PFE for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/middle
  if [[ "$genomeVer" == "hg19" ]]; then
      control_file=$script_dir/Epic-seq/code/priordata/control_hg19.txt
  elif [[ "$genomeVer" == "hg38" ]]; then
      control_file=$script_dir/Epic-seq/code/priordata/control_hg38.txt
  fi

  if [[ "$tssinfo" == "$script_dir/Epic-seq/code/priordata/sample_hg19.txt" && "$genomeVer" == "hg38" ]]; then
    tssinfo=$script_dir/Epic-seq/code/priordata/sample_hg38.txt
  fi 
  
  cat "$tssinfo" "$control_file" > $bam_output_dir/middle/tssinfo_unsorted.txt
  sort -k1,1V -k2,2n $bam_output_dir/middle/tssinfo_unsorted.txt > $bam_output_dir/middle/tssinfo.txt
  awk '{print $1 "\t" $2-1000 "\t" $2+1000}' $bam_output_dir/middle/tssinfo.txt > $bam_output_dir/middle/panelbed.txt
  sed -i '1i\CHR\tTSS\tGene-Symbol\tCategory\tStrand\tTS_ID' $bam_output_dir/middle/tssinfo.txt
  
  bedtools intersect -a $inputBam -b $bam_output_dir/middle/panelbed.txt -wa > $bam_output_dir/middle/intersected.bam
  samtools index $bam_output_dir/middle/intersected.bam
  Rscript $script_dir/Epic-seq/code/runEPIC.R --epic_dir $script_dir/Epic-seq/code --bamdir $bam_output_dir/middle --tssinfo $bam_output_dir/middle/tssinfo.txt --panelbed $bam_output_dir/middle/panelbed.txt --targeted no --outdir $bam_output_dir/middle $addEpicSeq
  mkdir -p $bam_output_dir/PFE
  cut -f1,5 $bam_output_dir/middle/epicseq.features.merged.by.0.txt > $bam_output_dir/PFE/PFE.txt
  rm -r $bam_output_dir/middle
}

# Calculate TSS coverage
calculate_tssc() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating TSSC for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/TSSC

  output_file=$bam_output_dir/TSSC/TSS.bed

  # Create or clear the output file
  > "$output_file"

  # Process each line of the input file
  while IFS=$'\t' read -r chrom start end name score strand; do
      start=$(printf "%d" "$start")
      end=$(printf "%d" "$end")
      upstream_site=0
      downstream_site=0

      if [ "$strand" == "+" ]; then
          upstream_site=$((start - upstream))
          if [ "$upstream_site" -lt 0 ]; then
              upstream_site=0
          fi
          downstream_site=$((end + downstream))
      echo -e "${chrom}\t${upstream_site}\t${downstream_site}\t${name}\t${score}\t${strand}" >> "$output_file"
      elif [ "$strand" == "-" ]; then
          upstream_site=$((end + upstream))
          downstream_site=$((start - downstream))
          if [ "$downstream_site" -lt 0 ]; then
              downstream_site=0
          fi
      echo -e "${chrom}\t${downstream_site}\t${upstream_site}\t${name}\t${score}\t${strand}" >> "$output_file"
      else
          echo "Error: Strand must be either '+' or '-'"
          exit 1
      fi
  done < "$TSSfile"
  if [[ "$sequencing" == "pair" ]];then
    bamCoverage -b $inputBam -o $bam_output_dir/TSSC/coverage.bw --extendReads $addbamCoverage
  elif [[ "$sequencing" == "single" ]];then
    bamCoverage -b $inputBam -o $bam_output_dir/TSSC/coverage.bw --extendReads 167 $addbamCoverage
  fi
  multiBigwigSummary BED-file -b $bam_output_dir/TSSC/coverage.bw -o $bam_output_dir/TSSC/coverage.npz --BED $bam_output_dir/TSSC/TSS.bed --outRawCounts $bam_output_dir/TSSC/coverage.txt $addmultiBigwigSummary
  sed '1d' $bam_output_dir/TSSC/coverage.txt | sort -k1,1V -k2,2n -k3,3n | awk 'BEGIN {print "chr\tstart\tend\tcoverage"} {print}' > $bam_output_dir/TSSC/average_coverage.txt
  rm $bam_output_dir/TSSC/coverage.txt
  rm $bam_output_dir/TSSC/coverage.bw
  rm $bam_output_dir/TSSC/coverage.npz
  rm $bam_output_dir/TSSC/TSS.bed
}

# Run the analysis
run_analysis