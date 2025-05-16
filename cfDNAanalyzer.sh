#!/bin/bash

usage(){
  cat << EOF
Description:
cfDNAanalyzer (cell-free DNA sequencing data analyzer) is a toolkit for cfDNA genomic sequencing data analysis which includes three core modules:
(1) Feature Extraction, which extracts multiple genomic and fragmentomic features at whole-genome or specific genomic-region levels;
(2) Feature Processing and Selection, allowing for refinement and optimization of extracted features, suporting multiple processing and selection methods;
(3) Machine Learning Model Building, supporting the development of both single-modality and multiple-modality predictive models for disease detection and classification.

Usage:
   bash cfDNAanalyzer.sh -I <InputFile> -o <OutputDirectory> -F <Features> [Options]
  
-------------------- Options for feature extraction --------------------

-----General options -----
  -I  FILE      A text file containing all input BAM files with one BAM file per line. 
                BAM files generated using both Bowtie2 and BWA are accepted.  
  -o  DIR       Output directory for all the results. Default: [./]
  -F  STR       Features to extract, including CNA, NOF, WPS, EM, EMR, FP, FPR, NP, OCF, PFE, and TSSC.
                Features should be set as a string separated by comma, e.g., CNA,NOF. 
                Default: All available features will be extracted.
                The detailed description of each feature can be accessed at https://github.com/LiymLab/cfDNAanalyzer. 
                Note: The following features are specifically designed for paired-end sequencing data: FP, FPR, EM, EMR, NP, PFE, and OCF.
  -g  STR       Genome version of input BAM files (hg19/hg38). Default: [hg38] 
  -b  FILE      A BED3 file specifying the regions to extract features.
                The file should contain three TAB-delimited columns: chromosome start end.
  -s  STR       Sequencing method of input BAM files (single/pair). Default: [pair]
  -t  INT       Number of threads to use. Default: [1]              

----- Options specific for Copy Number Alterations (CNA) -----
  -B     INT    Bin size in kilobases (10, 50, 500, or 1000). Default: [1000]
  --CNA  STR    Additional parameter setting for CNA analysis. 
                The full parameter list is available by running: Rscript cfDNAanalyzer/ichorCNA/ichorCNA/scripts/runIchorCNA.R --help. [optional]

----- Options specific for Nucleosome Occupancy and Fuzziness (NOF) -----
  --NOF  STR    Additional parameter setting for NOF analysis. 
                The full parameter list is available by running: python cfDNAanalyzer/DANPOS3/danpos.py dpos -h. [optional]

----- Options specific for Windowed Protection Score (WPS) -----
  -x  INT       Min fragment length used for long fragments WPS calculation. Default: [120]
  -X  INT       Max fragment length used for long fragments WPS calculation. Default: [180]
  -w  INT       Window size for long fragments WPS calculation. Default: [120]
  -m  INT       Min fragment length used for short fragments WPS calculation. Default: [35]
  -M  INT       Max fragment length used for short fragments WPS calculation. Default: [80]
  -W  INT       Window size for short fragments WPS calculation. Default: [16]

----- Options specific for End Motif frequency and diversity (EM) -----
  -f  FILE      Reference genome in FASTA format. 
                For example, hg38 reference genome can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz.

----- Options specific for Nucleosome Profile (NP) -----
  -l  DIR       Directory containing a list of files with each file for a set of sites.
                The file must have at least two columns with the following column names: 
                "Chrom" for the chromosome name and "position" for the site positions. 
                If not provided, the 377 TF binding site lists (cfDNAanalyzer/Griffin/Ref_hg19/sites or cfDNAanalyzer/Griffin/Ref_hg38/sites) from the original Nucleosome Profile paper will be used . 

----- Options for Promoter Fragmentation Entropy (PFE) -----
  -T     FILE   A TAB-delimited TSS information file without any header.
                The file must have six columns: (1) chromosome, (2) 1-base TSS coordinate, (3) gene name, (4) strand (+1/-1) and (5) TSS ID (e.g., for genes with multiple TSS, this could be geneName_1, geneName_2, etc.)
                If not provided, TSSs (cfDNAanalyzer/Epic-seq/code/priordata/sample_hg19.txt or cfDNAanalyzer/Epic-seq/code/priordata/sample_hg38.txt) from the original Promoter Fragmentation Entropy paper will be used.
  --PFE  STR    Addtional parameter setting for PFE analysis. 
                The full parameter list is available by running: Rscript cfDNAanalyzer/Epic-seq/code/runEPIC.R -h.[optional]

----- Options for TSS Coverage (TSSC) -----
  -u                    INT   Number of base pairs upstream of TSS used for calculating TSSC. Default: [1000]
  -d                    INT   Number of base pairs downstream of TSS used for calculating TSSC. Default: [1000]
  -S                    FILE  A BED6 file specifying the coordinates of TSSs used for calculating TSSC. 
                              If not provided, all TSSs of refSeq genes will be used.
  -n                    STR   Methods to normalize the number of reads per bin for the "bamCoverage" command of deeptools (RPKM,CPM,BPM,RPGC) . Default: [RPKM]
  --bamCoverage         STR   Additional parameter seting for the "bamCoverage" command of deeptools.
                              The full parameter list is available by running: bamCoverage -h. [optional] 
  --multiBigwigSummary  STR   Additional parameter seting for the "multiBigwigSummary" command of deeptools.
                              The full parameter list is available by running: multiBigwigSummary -h.[optional].
                          

-------------------- Options for downstream analysis --------------------
  --noDA                LOGIC  Skip all the downstream analysis.
  --noML                LOGIC  Skip machine learning model building, only feature processing and selection will be conducted.
                               Methods in feature selection will leverage the label information of all samples.
  --labelFile           FILE   Label information file for all the samples, seperated by comma.
                               This file must have two columns. The sample column contains the basename of all the samples. The label column contains all the samples' label information (like 0,1 for two class and 0,1,2 for three class).

----- Options for feature processing and selection -----
  --standardM                    STR    The standardization method (Zscore, MMS, RS, or QT). Default:   [Zscore]
                                        The detailed description of each method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --filterMethod                 STR    Filter methods employed for feature selection (IG CHI FS FCBF PI CC LVF MAD DR MI RLF SURF MSURF TRF).
                                        Methods should be set as a string separated by space in quotes, e.g., 'MR IG CHI'. 
                                        Default: [IG].
                                        The detailed description of each filter method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --filterNum                   FLOAT   Number of features to retain when employing the filter method. 
                                        If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default: [0.2]
  --wrapperMethod                STR    Wrapper methods employed for feature selection (FFS BFS EFS RFS BOR).
                                        Methods should be set as a string separated by space in quotes, e.g., 'BOR RFS'. 
                                        Default: [BOR].
                                        The detailed description of each wrapper method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --wrapperNum                  FLOAT   Number of features to retain when employing the wrapper method.
                                        If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default: [0.2]
  --embeddedMethod               STR    Embedded methods employed for feature selection (LASSO RIDGE ELASTICNET RF).
                                        Methods should be set as a string separated by space in quotes, e.g., 'LASSO RIDGE'. 
                                        Default: [LASSO]
                                        The detailed description of each embedded method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --embeddedNum                 FLOAT   Number of features to retain when employing the embedded method. 
                                        If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default: [0.2]
  --hybridType                   STR    Two methods used for hybrid method (FW/FE). Default: [FE]
  --hybridMethod1,hybridMethod2  STR    Subtype methods designated for method 1 and method 2 in the "--hybridType".
                                        Methods should be set as a string separated by space in quotes, e.g., 'BOR RFS'.
                                        Default: [IG,BOR] for FW; [IG,LASSO] for FE.
  --hybridNum1,hybridNum2       FLOAT   Number of features to retain for method 1 and method 2 in the "--hybridType". 
                                        If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default: [0.2,0.2]
  
----- Options for machine learning model building -----
  --classNum          INT   Number of classification categories (2 or more). Default: [2]
  --cvSingle          STR   Cross-validation method applied in single modality, options include leave-one-out (LOO) or K-fold (KFold). Default: [LOO]
  --nsplitSingle      INT   Number of folds designated for K-fold cross-validation in single modality. Default: [5]
  --classifierSingle  STR   Classifiers employed in single modality (KNN SVM RandomForest GaussianNB LogisticRegression XGB).
                            Classifiers should be set as a string separated by space in quotes, e.g., 'KNN XGB'.
                            Default: All available classifiers will be applied.
  --cvMulti           STR   Cross-validation method applied in multiple modalities, options include leave-one-out (LOO) or K-fold (KFold). Default: [LOO]
  --nsplitMulti       INT   Number of folds designated for K-fold cross-validation in multiple modalities. Default: [5]
  --classifierMulti   STR   Classifiers employed in multiple modalities (KNN SVM RandomForest GaussianNB LogisticRegression XGB).
                            Classifiers should be set as a string separated by space in quotes, e.g., 'KNN XGB'.
                            Default: All available classifiers will be applied.
  --modelMethod       STR   Model-based integration methods employed in multiple modality. (average weighted majority stack) 
                            Methods should be set as a string separated by space in quotes, e.g., 'average weighted'. 
                            Default: [average]
                            The detailed description of each method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --transMethod       STR   Transformation-based integration methods employed in multiple modality. (pca linear polynomial rbf sigmoid snf) 
                            Methods should be set as a string separated by space in quotes, e.g., 'pca linear'. 
                            Default: [pca]
                            The detailed description of each method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.     
EOF
    exit 0
}

if [ $# -eq 0 ] || ([ $# -eq 1 ] && { [ "$1" == "-h" ] || [ "$1" == "--help" ]; }); then
    usage
fi

script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"

# feature extraction
inputBamList=""  # -I
outputDir=$(pwd)   # -o
Feature=CNA,EM,FP,NOF,NP,WPS,OCF,EMR,FPR,PFE,TSSC  # -F
genomeVer=hg38  # -g
fasta=""   # -f
sequencing=pair  # -s
bedfiles=""  # -b
threads=1  # -t
binsize=1000  # -B
upstream=1000  # -u
downstream=1000  # -d
TSSfile=""  # -S
Minlength_long=120  # -x
Maxlength_long=180  # -X
window_long=120  # -w
Minlength_short=35  # -m
Maxlength_short=80  # -M
window_short=16  # -W
sitespath=""  # -l
tssinfo=""  # -T
normalization="RPKM"  # -n 
addCNA=""  # --CNA
addNOF=""  # --NOF
addbamCoverage=""  # --bamCoverage
addmultiBigwigSummary=""  # --multiBigwigSummary
addEpicSeq=""   # --PFE

# downstream analysis
noDA=0
noML=0
labelFile=""
standardM=Zscore
filterMethod="IG"
filterNum=0.2
wrapperMethod="BOR"
wrapperNum=0.2
embeddedMethod="LASSO"
embeddedNum=0.2
hybridType="FE"
hybridMethod1=""
hybridMethod2=""
hybridNum1=0.2
hybridNum2=0.2
classNum=2
cvSingle=LOO
nsplitSingle=5
classifierSingle="KNN SVM RandomForest GaussianNB LogisticRegression XGB"
cvMulti=LOO
nsplitMulti=5
classifierMulti="KNN SVM RandomForest GaussianNB LogisticRegression XGB"
modelMethod=average
transMethod=pca

# Use getopt to handle long options
TEMP=$(getopt -o hI:o:F:g:f:s:b:t:B:u:d:S:x:X:w:m:M:W:l:T:n: \
               --long help,CNA:,NOF:,bamCoverage:,multiBigwigSummary:,PFE:,noDA,noML,labelFile:,standardM:,filterMethod:,filterNum:,wrapperMethod:,wrapperNum:,embeddedMethod:,embeddedNum:,hybridType:,hybridMethod1:,hybridMethod2:,hybridNum1:,hybridNum2:,classNum:,cvSingle:,nsplitSingle:,classifierSingle:,cvMulti:,nsplitMulti:,classifierMulti:,modelMethod:,transMethod: \
               -- "$@")

# Check whether getopt succeeded
if [ $? != 0 ]; then
    echo "Error in command line parsing." >&2
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
        -t)                       threads=$2; shift 2;;        
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
        -n)                       normalization=$2; shift 2;;        
        --CNA)                    addCNA=$2; shift 2;;
        --NOF)                    addNOF=$2; shift 2;;
        --bamCoverage)            addbamCoverage=$2; shift 2;;
        --multiBigwigSummary)     addmultiBigwigSummary=$2; shift 2;;
        --PFE)                    addEpicSeq=$2; shift 2;;
        --noDA)                   noDA=1; shift;;
        --noML)                   noML=1; shift;;
        --labelFile)              labelFile=$2; shift 2;;
        --standardM)              standardM=$2; shift 2;;
        --filterMethod)           filterMethod=$2; shift 2;;
        --filterNum)              filterNum=$2; shift 2;;
        --wrapperMethod)          wrapperMethod=$2; shift 2;;
        --wrapperNum)             wrapperNum=$2; shift 2;;
        --embeddedMethod)         embeddedMethod=$2; shift 2;;
        --embeddedNum)            embeddedNum=$2; shift 2;;
        --hybridType)             hybridType=$2; shift 2;;
        --hybridMethod1)          hybridMethod1=$2; shift 2;;
        --hybridMethod2)          hybridMethod2=$2; shift 2;;
        --hybridNum1)             hybridNum1=$2; shift 2;;
        --hybridNum2)             hybridNum2=$2; shift 2;;
        --classNum)               classNum=$2; shift 2;;
        --cvSingle)               cvSingle=$2; shift 2;;
        --nsplitSingle)           nsplitSingle=$2; shift 2;;
        --classifierSingle)       classifierSingle=$2; shift 2;;
        --cvMulti)                cvMulti=$2; shift 2;;
        --nsplitMulti)            nsplitMulti=$2; shift 2;;
        --classifierMulti)        classifierMulti=$2; shift 2;;
        --modelMethod)            modelMethod=$2; shift 2;;
        --transMethod)            transMethod=$2; shift 2;;
        --)                       shift; break;;
        *)                        echo "Unexpected option: $1"; exit 1;;
    esac
done

# Check if jq is installed
if ! command -v jq &> /dev/null
then
    echo "Error: jq is required but not installed. Please install jq and rerun the script."
    exit 1
fi
# You can download ./jq with `sudo apt-get install jq`
# Reference: https://jqlang.github.io/jq/download/ 

# Main function
run_analysis() {
  # Ensure input file list is provided
  if [ -z "$inputBamList" ]; then
    echo "Error: Input bam file list is not provided."
    exit 1
  fi

  # Create output directory if it doesn't exist
  mkdir -p "$outputDir"

  # Validate reference type
  if [[ "$genomeVer" != "hg19" ]] && [[ "$genomeVer" != "hg38" ]]; then
    echo "Error: Incorrect type of reference fasta file."
    exit 1
  fi
  
  if [[ "$Feature" == *"NOF"* ]] || [[ "$Feature" == *"NP"* ]] || [[ "$Feature" == *"WPS"* ]] || [[ "$Feature" == *"OCF"* ]] || [[ "$Feature" == *"EMR"* ]] || [[ "$Feature" == *"FPR"* ]] && [ -z "$bedfiles" ]; then
    echo "Error: A BED3 file specifying the regions (-b) is not provided."
    exit 1
  fi
  
  # Convert relative paths to absolute paths
  inputBamList=$(realpath "$inputBamList")
  outputDir=$(realpath "$outputDir")
  if [ -n "$fasta" ]; then
    fasta=$(realpath "$fasta")
  fi
  
  if [ -n "$bedfiles" ]; then
    bedfiles=$(realpath "$bedfiles")
  fi  
    
  # Process each bam file
  for inputBam in $(cat "$inputBamList"); do
      if [ ! -f "$inputBam" ]; then
        echo "Error: Input bam file was not found: $inputBam"
        exit 1
      fi
      
      start_time=$(date +%s)
      
      samtools index -@ $threads "$inputBam"

      local filename=$(basename "$inputBam" .bam)
      local bam_output_dir="$outputDir/Samples/$filename"
      mkdir -p "$bam_output_dir"
      
      # Features for whole genome
      if [[ ",$Feature," == *",CNA,"* ]]; then calculate_cna "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",EM,"* ]]; then calculate_em "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",FP,"* ]]; then calculate_fp "$inputBam" "$bam_output_dir"; fi

      # Features for specific regions
      if [[ ",$Feature," == *",NOF,"* ]]; then calculate_nof "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",NP,"* ]]; then calculate_np "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",WPS,"* ]]; then calculate_wps "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",OCF,"* ]]; then calculate_ocf "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",EMR,"* ]]; then calculate_emr "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",FPR,"* ]]; then calculate_fpr "$inputBam" "$bam_output_dir"; fi

      # Features for transcription start sites
      if [[ ",$Feature," == *",PFE,"* ]]; then calculate_pfe "$inputBam" "$bam_output_dir"; fi
      if [[ ",$Feature," == *",TSSC,"* ]]; then calculate_tssc "$inputBam" "$bam_output_dir"; fi
      echo "
      Analysis for file: $inputBam completed."
      
      end_time=$(date +%s)
      elapsed_time=$(( end_time - start_time ))
      echo "This tool took $elapsed_time seconds for $inputBam"
  done

  # Generate config file for selected features
  config_all_file="$script_dir/Feature_Processing/config_all_feature.json"
  output_config_file="$outputDir/Samples/config_feature.json"
  generate_config "$Feature" "$config_all_file" "$output_config_file"

  wait
}


# Function to generate config file for selected features
generate_config() {
  local features_list=$1  # comma-separated list of features
  local config_all_file=$2  # path to config_all_feature.json
  local output_config_file=$3  # path to output config file

  # Convert features_list to an array
  IFS=',' read -r -a features <<< "$features_list"

  # Declare associative array for feature to keys mapping
  declare -A feature_to_keys
  feature_to_keys=( 
    [CNA]="CNA"
    [EM]="EM_motifs_frequency EM_motifs_mds"
    [EMR]="EMR_aggregated_mds EMR_aggregated_motif_frequency EMR_region_mds EMR_region_motif_frequency"
    [FP]="FP_fragmentation_profile"
    [FPR]="FPR_fragmentation_profile_regions"
    [NOF]="NOF_meanfuziness NOF_occupancy"
    [NP]="NP_mean_coverage NP_central_coverage NP_amplitude"
    [OCF]="OCF"
    [PFE]="PFE"
    [TSSC]="TSSC_average_coverage"
    [WPS]="WPS_long WPS_short"
  )

  # Build a list of keys to extract from the config_all_feature.json
  keys_to_extract=()
  for feature in "${features[@]}"; do
    keys="${feature_to_keys[$feature]}"
    if [ -z "$keys" ]; then
      echo "Warning: No config keys found for feature '$feature'"
      continue
    fi
    for key in $keys; do
      keys_to_extract+=("$key")
    done
  done

  if [ ${#keys_to_extract[@]} -eq 0 ]; then
    echo "Error: No config keys to extract"
    exit 1
  fi

  # Build the jq filter to extract the desired keys
  jq_filter="{"
  for key in "${keys_to_extract[@]}"; do
    jq_filter+="\"$key\": .\"$key\", "
  done
  jq_filter="${jq_filter%, }"  # Remove trailing comma and space
  jq_filter+="}"

  # Use jq to extract the desired keys
  if ! [ -f "$config_all_file" ]; then
    echo "Error: config_all_feature.json file not found at $config_all_file"
    exit 1
  fi

  jq "$jq_filter" "$config_all_file" > "$output_config_file"

  echo "Generated config file with selected features: $output_config_file"
}


# Calculate Copy Number Alterations
calculate_cna() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating CNA for file: '$inputBam'*****'
  if [[ "$binsize" != "10" ]] && [[ "$binsize" != "50" ]] && [[ "$binsize" != "500" ]] && [[ "$binsize" != "1000" ]]; then
    echo "Error: incorrect bin size."
    exit 1
  fi
  mkdir -p $bam_output_dir/CNA
  $script_dir/ichorCNA/hmmcopy_utils/bin/readCounter --window $((binsize * 1000)) --quality 20 --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $inputBam > $bam_output_dir/CNA/CNA.wig
  if [[ "$genomeVer" == "hg19" ]];then
    Rscript $script_dir/ichorCNA/ichorCNA/scripts/runIchorCNA.R --genomeBuild hg19 --WIG $bam_output_dir/CNA/CNA.wig --gcWig $script_dir/ichorCNA/ichorCNA/inst/extdata/gc_hg19_${binsize}kb.wig --mapWig $script_dir/ichorCNA/ichorCNA/inst/extdata/map_hg19_${binsize}kb.wig --centromere $script_dir/ichorCNA/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --outDir $bam_output_dir/CNA  $addCNA 
  elif [[ "$genomeVer" == "hg38" ]];then
    Rscript $script_dir/ichorCNA/ichorCNA/scripts/runIchorCNA.R --genomeBuild hg38 --WIG $bam_output_dir/CNA/CNA.wig --gcWig $script_dir/ichorCNA/ichorCNA/inst/extdata/gc_hg38_${binsize}kb.wig --mapWig $script_dir/ichorCNA/ichorCNA/inst/extdata/map_hg38_${binsize}kb.wig --centromere $script_dir/ichorCNA/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt --outDir $bam_output_dir/CNA $addCNA
  fi
  cut -f1,2,3,4 $bam_output_dir/CNA/test.cna.seg > $bam_output_dir/CNA/CNA.txt
  rm -r $bam_output_dir/CNA/test
  rm $bam_output_dir/CNA/test.cna.seg
  rm $bam_output_dir/CNA/test.correctedDepth.txt
  rm $bam_output_dir/CNA/test.params.txt
  rm $bam_output_dir/CNA/test.RData
  rm $bam_output_dir/CNA/test.seg
  rm $bam_output_dir/CNA/test.seg.txt
  rm $bam_output_dir/CNA/CNA.wig
}

# Calculate End motif frequency and diversity
calculate_em() {
  if [[ "$sequencing" == "single" ]] ; then
    echo "Error: Feature EM is specifically designed for paired-end sequencing data."
    exit 1
  fi
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating EM for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/EM
  samtools sort -@ $threads -n -o $bam_output_dir/EM/sorted_by_name.bam $inputBam
  bedtools bamtobed -bedpe -i $bam_output_dir/EM/sorted_by_name.bam > $bam_output_dir/EM/${filename}.bedpe
  awk ' BEGIN {OFS="\t"}
  {
    min_2_5 = $2
    if ($5 < min_2_5) min_2_5 = $5

    max_3_6 = $3
    if ($6 > max_3_6) max_3_6 = $6

    print $1, min_2_5, max_3_6
  }
  ' $bam_output_dir/EM/${filename}.bedpe > $bam_output_dir/EM/${filename}.bed
  python $script_dir/End_motif_frequency/endMotifFreq.py -b $bam_output_dir/EM/${filename}.bed -f $fasta -o1 $bam_output_dir/EM/motifs_frequency.txt -o2 $bam_output_dir/EM/motifs_mds.txt
  rm $bam_output_dir/EM/sorted_by_name.bam
  rm $bam_output_dir/EM/${filename}.bedpe
  rm $bam_output_dir/EM/${filename}.bed
}

# Calculate Fragmentation Profile
calculate_fp() {
  if [[ "$sequencing" == "single" ]] ; then
    echo "Error: Feature FP is specifically designed for paired-end sequencing data"
    exit 1
  fi
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating FP for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/FP
  if [[ "$genomeVer" == "hg19" ]];then
    Rscript $script_dir/Fragmentation_profile/fragmentProfile_delfi.R -i=$inputBam -o=$bam_output_dir/FP/FragmentationProfile.txt $script_dir/Fragmentation_profile/hg19_100kb_WG.txt $genomeVer $script_dir
  elif [[ "$genomeVer" == "hg38" ]];then
    Rscript $script_dir/Fragmentation_profile/fragmentProfile_delfi.R -i=$inputBam -o=$bam_output_dir/FP/FragmentationProfile.txt $script_dir/Fragmentation_profile/hg38_100kb_WG.txt $genomeVer $script_dir
  fi
  cut -f1,2,3,6,7,8 $bam_output_dir/FP/FragmentationProfile.txt > $bam_output_dir/FP/Fragmentation_Profile.txt
  rm $bam_output_dir/FP/FragmentationProfile.txt
}


# Calculate Nucleosome Occupancy and Fuzziness
calculate_nof() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating NOF for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/NOF
  python $script_dir/DANPOS3/danpos.py dpos $inputBam -o $bam_output_dir/NOF $addNOF
  if [[ "$genomeVer" == "hg19" ]];then
    $script_dir/DANPOS3/wigToBigWig -clip $bam_output_dir/NOF/pooled/*.wig $script_dir/DANPOS3/hg19.chrom.sizes $bam_output_dir/NOF/NOF.bw
  elif [[ "$genomeVer" == "hg38" ]];then
    $script_dir/DANPOS3/wigToBigWig -clip $bam_output_dir/NOF/pooled/*.wig $script_dir/DANPOS3/hg38.chrom.sizes $bam_output_dir/NOF/NOF.bw
  fi
  awk 'BEGIN {OFS="\t"} {print $0, "region_" NR}' $bedfiles > $bam_output_dir/NOF/regions_with_ids.bed
  $script_dir/DANPOS3/bigWigAverageOverBed $bam_output_dir/NOF/NOF.bw  $bam_output_dir/NOF/regions_with_ids.bed $bam_output_dir/NOF/occupancy.tsv -bedOut=$bam_output_dir/NOF/occupancy.bed
  sed -i '1i\chr\tstart\tend\tID\toccupancy' $bam_output_dir/NOF/occupancy.bed
  awk -v OFS=$'\t' 'NR>1{print $1,$4-1,$4,$5,$6}' $bam_output_dir/NOF/pooled/*.xls |grep -vE "Lambda_NEB|pUC19"| bedtools intersect -a $bedfiles -b stdin -wo > $bam_output_dir/NOF/intersect.tsv 
  awk -v OFS=$'\t' '{sum[$1$2$3]+=$8; count[$1$2$3]++; total[$1$2$3]=$1" "$2" "$3} END {for (key in sum) print total[key], sum[key]/count[key]}' $bam_output_dir/NOF/intersect.tsv | sort -k1,1V -k2,2n -k3,3n | awk 'BEGIN {print "chr\tstart\tend\tmeanfuziness"} {print}' > $bam_output_dir/NOF/meanfuziness.tsv
  rm $bam_output_dir/NOF/intersect.tsv
  rm $bam_output_dir/NOF/NOF.bw
  rm $bam_output_dir/NOF/regions_with_ids.bed
  rm $bam_output_dir/NOF/occupancy.tsv
}

# Calculate Nucleosome Profile
calculate_np() {
  if [[ "$sequencing" == "single" ]] ; then
    echo "Error: Feature NP is specifically designed for paired-end sequencing data."
    exit 1
  fi
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating NP for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/NP
  # Defining path variables
  griffin_path="$script_dir/Griffin"
  if [[ "$genomeVer" == "hg19" ]];then
    Ref_dir="$griffin_path/Ref_hg19"
    if [[ -z "$sitespath" ]]; then
      sitespath=$script_dir/Griffin/Ref_hg19/sites
    fi
  elif [[ "$genomeVer" == "hg38" ]];then
    Ref_dir="$griffin_path/Ref_hg38"
    if [[ -z "$sitespath" ]]; then
      sitespath=$script_dir/Griffin/Ref_hg38/sites
    fi
  fi

  result_dir="$griffin_path/snakemakes"
  gc_correction_dir="$result_dir/griffin_GC_and_mappability_correction"
  nucleosome_profiling_dir="$result_dir/griffin_nucleosome_profiling"
        
  # Preparing for GC and mappability corrections
  sed -e "s!Griffin_path!$griffin_path!g" -e "s!reference_path!$fasta!g" -e "s!Ref_dir!$Ref_dir!g" \
      -e "s!result_name!$gc_correction_dir/result_${filename}!g" "$gc_correction_dir/config/config.yaml" \
      > "$gc_correction_dir/config/config_${filename}.yaml"
        
  sed "s!input!$inputBam!g" "$gc_correction_dir/config/samples.yaml" \
      > "$gc_correction_dir/config/samples_${filename}.yaml"
        
  sed -e "s!config_samples!$gc_correction_dir/config/samples_${filename}.yaml!g" \
      -e "s!config_config!$gc_correction_dir/config/config_${filename}.yaml!g" \
      -e "s!config_cluster!$gc_correction_dir/config/cluster_slurm.yaml!g" \
      "$gc_correction_dir/griffin_GC_and_mappability_correction.snakefile" \
      > "$gc_correction_dir/griffin_GC_and_mappability_correction_${filename}.snakefile"
        
  snakemake -s "$gc_correction_dir/griffin_GC_and_mappability_correction_${filename}.snakefile" --cores $threads
        
  rm "$gc_correction_dir/config/config_${filename}.yaml" "$gc_correction_dir/config/samples_${filename}.yaml"
  rm "$gc_correction_dir/griffin_GC_and_mappability_correction_${filename}.snakefile"
  mv "$gc_correction_dir/result_${filename}/samples.GC.yaml" "$nucleosome_profiling_dir/config/samples.GC_${filename}.yaml"
        
  # Nucleosome profiling analysis
  sed -e "s!Griffin_path!$griffin_path!g" -e "s!reference_path!$fasta!g" -e "s!Ref_dir!$Ref_dir!g" \
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
        
  snakemake -s "$nucleosome_profiling_dir/griffin_nucleosome_profiling_${filename}.snakefile" --cores $threads --unlock
  snakemake -s "$nucleosome_profiling_dir/griffin_nucleosome_profiling_${filename}.snakefile" --cores $threads
        
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
  cp -f "$nucleosome_profiling_dir/result_${filename}/NucleosomeProfile.txt" "$bam_output_dir/NP/NucleosomeProfile.txt"
  cp -f -r "$nucleosome_profiling_dir/result_${filename}/plots" "$bam_output_dir/NP/plots"
  cp -f -r "$nucleosome_profiling_dir/result_${filename}/sample_name_1" "$bam_output_dir/NP/sample_name_1"
  rm -r "$nucleosome_profiling_dir/result_${filename}"
}

# Calculate Windowed Protection Score
calculate_wps() {
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating WPS for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/WPS
  echo -e "chr\tstart\tend\tlong_WPS\tshort_WPS" >$bam_output_dir/WPS/WPS.txt
  # Extract the first three columns of each line of the file (chr,start,end)
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
  rm $bam_output_dir/WPS/*.wig
}

# Calculate Orientation-aware CfDNA Fragmentation
calculate_ocf() {
  if [[ "$sequencing" == "single" ]] ; then
    echo "Error: Feature OCF is specifically designed for paired-end sequencing data"
    exit 1
  fi
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating OCF for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/OCF
  sort -k1,1 -k2,2n $bedfiles | awk 'BEGIN {OFS="\t"} {print $0, "region_" NR}' > $bam_output_dir/OCF/regions_with_ids.bed
  samtools sort -@ $threads -n -o $bam_output_dir/OCF/sorted_by_name.bam $inputBam
  bedtools bamtobed -bedpe -i $bam_output_dir/OCF/sorted_by_name.bam > $bam_output_dir/OCF/${filename}.bedpe
  awk ' BEGIN {OFS="\t"}
  {
    min_2_5 = $2
    if ($5 < min_2_5) min_2_5 = $5
  
    max_3_6 = $3
    if ($6 > max_3_6) max_3_6 = $6
  
    print $1, min_2_5, max_3_6
  }
  ' $bam_output_dir/OCF/${filename}.bedpe | sort -k1,1 -k2,2n > $bam_output_dir/OCF/sorted_${filename}.bed
  $script_dir/OCF/OCF.sh $bam_output_dir/OCF/regions_with_ids.bed $bam_output_dir/OCF/sorted_${filename}.bed sample $bam_output_dir/OCF
  awk -v OFS="\t" '{n=split($1,a,".");print a[n-2],$2}' $bam_output_dir/OCF/sample.sorted_${filename}.bed.OCF | \
  awk '{print substr($1, 8), $0}' | sort -n | cut -d' ' -f2- | \
  awk 'NR==FNR{a[$1]=$0; next} $4 in a {print a[$4], $0}' - $bam_output_dir/OCF/regions_with_ids.bed | \
  awk '{print $3 "\t" $4 "\t" $5 "\t" $2}' | sed '1i\chr\tstart\tend\tOCF' > $bam_output_dir/OCF/OCF.txt
  rm $bam_output_dir/OCF/regions_with_ids.bed
  rm $bam_output_dir/OCF/sorted_by_name.bam
  rm $bam_output_dir/OCF/${filename}.bedpe
  rm $bam_output_dir/OCF/sorted_${filename}.bed
  rm $bam_output_dir/OCF/sample.sorted_${filename}.bed.OCF
  rm $bam_output_dir/OCF/sample.ol.sorted_${filename}.bed
}

# Calculate End Motif frequency and diversity for Regions
calculate_emr() {
  if [[ "$sequencing" == "single" ]] ; then
    echo "Error: Feature EMR is specifically designed for paired-end sequencing data"
    exit 1
  fi
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating EMR for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/EMR
  
  # aggregated
  bedtools intersect -abam $inputBam -b $bedfiles > $bam_output_dir/EMR/intersect.bam
  samtools sort -@ $threads -n -o $bam_output_dir/EMR/sorted_intersect.bam $bam_output_dir/EMR/intersect.bam 
  bedtools bamtobed -bedpe -i $bam_output_dir/EMR/sorted_intersect.bam > $bam_output_dir/EMR/${filename}.bedpe
  awk ' BEGIN {OFS="\t"}
  {
    min_2_5 = $2
    if ($5 < min_2_5) min_2_5 = $5
  
    max_3_6 = $3
    if ($6 > max_3_6) max_3_6 = $6
  
    print $1, min_2_5, max_3_6
  }
  ' $bam_output_dir/EMR/${filename}.bedpe > $bam_output_dir/EMR/${filename}.bed
  python $script_dir/End_motif_frequency/endMotifFreq.py -b  $bam_output_dir/EMR/${filename}.bed -f $fasta -o1 $bam_output_dir/EMR/aggregated_motif_frequency.txt -o2 $bam_output_dir/EMR/aggregated_mds_noheader.txt
  sed -i '1iMotif\tFrequency' $bam_output_dir/EMR/aggregated_motif_frequency.txt
  sed -n '2p' $bam_output_dir/EMR/aggregated_mds_noheader.txt | awk 'NR==1{print "MDS"} {print}' > $bam_output_dir/EMR/aggregated_mds.txt
  
  
  # for all the regions
  python $script_dir/End_motif_frequency/endMotifFreq_all_regions.py -b1 $bedfiles -b2 $bam_output_dir/EMR/${filename}.bed -f $fasta -o1  $bam_output_dir/EMR/region_mds.txt -o2 $bam_output_dir/EMR/region_motif_frequency.txt
  sed -i '1ichr\tstart\tend\tMDS' $bam_output_dir/EMR/region_mds.txt
  sed -i '1ichr\tstart\tend\tMotif\tFrequency'  $bam_output_dir/EMR/region_motif_frequency.txt
  rm $bam_output_dir/EMR/aggregated_mds_noheader.txt
  rm $bam_output_dir/EMR/intersect.bam
  rm $bam_output_dir/EMR/sorted_intersect.bam
  rm $bam_output_dir/EMR/${filename}.bedpe
  rm $bam_output_dir/EMR/${filename}.bed
  
}

# Calculate Fragmentation Profile for Regions
calculate_fpr() {
  if [[ "$sequencing" == "single" ]] ; then
    echo "Error: Feature FPR pecifically designed for paired-end sequencing data"
    exit 1
  fi
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating FPR for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/FPR
  awk 'BEGIN {print "chr\tstart\tend"} {print}' $bedfiles > $bam_output_dir/FPR/regions_with_header.bed
  Rscript $script_dir/Fragmentation_profile/fragmentProfile_delfi.R -i=$inputBam -o=$bam_output_dir/FPR/FragmentationProfile.txt $bam_output_dir/FPR/regions_with_header.bed $genomeVer $script_dir
  cut -f1,2,3,6,7,8 $bam_output_dir/FPR/FragmentationProfile.txt > $bam_output_dir/FPR/Fragmentation_Profile_regions.txt
  rm $bam_output_dir/FPR/FragmentationProfile.txt
  rm $bam_output_dir/FPR/regions_with_header.bed
}

# Calculate Promoter Fragmentation Entropy
calculate_pfe() {
  if [[ "$sequencing" == "single" ]] ; then
    echo "Error: Feature PFE is specifically designed for paired-end sequencing data"
    exit 1
  fi
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating PFE for file: '$inputBam'*****'
  mkdir -p $bam_output_dir/middle
  
  if [[ "$genomeVer" == "hg19" ]]; then
      control_file=$script_dir/Epic-seq/code/priordata/control_hg19.txt
      if [[ -z "$tssinfo" ]]; then
        tssinfo=$script_dir/Epic-seq/code/priordata/sample_hg19.txt
      fi
  elif [[ "$genomeVer" == "hg38" ]]; then
      control_file=$script_dir/Epic-seq/code/priordata/control_hg38.txt
      if [[ -z "$tssinfo" ]]; then
        tssinfo=$script_dir/Epic-seq/code/priordata/sample_hg38.txt
      fi
  fi 
  
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "sample", $4, $5}' $tssinfo  > $bam_output_dir/middle/tssinfo_5col.txt
  cat "$bam_output_dir/middle/tssinfo_5col.txt" "$control_file" | sort -k1,1V -k2,2n | sed '1i\CHR\tTSS\tGene-Symbol\tCategory\tStrand\tTS_ID' > $bam_output_dir/middle/tssinfo.txt
  awk 'NR>1 {print $1 "\t" $2-1000 "\t" $2+1000}' $bam_output_dir/middle/tssinfo.txt > $bam_output_dir/middle/uncorrect_1_panelbed.txt
  if [[ "$genomeVer" == "hg19" ]];then
    awk 'BEGIN {OFS="\t"} NR==FNR {sizes[$1]=$2; next} {if ($3 > sizes[$1]) $3=sizes[$1]; print}' $script_dir/DANPOS3/hg19.chrom.sizes $bam_output_dir/middle/uncorrect_1_panelbed.txt > $bam_output_dir/middle/uncorrect_2_panelbed.txt
  elif [[ "$genomeVer" == "hg38" ]];then
    awk 'BEGIN {OFS="\t"} NR==FNR {sizes[$1]=$2; next} {if ($3 > sizes[$1]) $3=sizes[$1]; print}' $script_dir/DANPOS3/hg38.chrom.sizes $bam_output_dir/middle/uncorrect_1_panelbed.txt > $bam_output_dir/middle/uncorrect_2_panelbed.txt
  fi
  awk 'BEGIN {OFS="\t"} {if ($2 < 0) $2 = 0; print}' $bam_output_dir/middle/uncorrect_2_panelbed.txt > $bam_output_dir/middle/panelbed.txt

  bedtools intersect -a $inputBam -b $bam_output_dir/middle/panelbed.txt -wa > $bam_output_dir/middle/intersected.bam
  samtools index $bam_output_dir/middle/intersected.bam
  Rscript $script_dir/Epic-seq/code/runEPIC.R --epic_dir $script_dir/Epic-seq/code --bamdir $bam_output_dir/middle --tssinfo $bam_output_dir/middle/tssinfo.txt --panelbed $bam_output_dir/middle/panelbed.txt --outdir $bam_output_dir/middle $addEpicSeq
  mkdir -p $bam_output_dir/PFE
  cut -f2,5 $bam_output_dir/middle/pfe.matrix.merged.by.0.txt > $bam_output_dir/PFE/PFE.txt
  rm -r $bam_output_dir/middle
}

# Calculate TSS coverage 
calculate_tssc() { 
  local inputBam=$1
  local bam_output_dir=$2
  echo '
  *****Calculating TSSC for file: '$inputBam'*****'
  
  if [[ "$normalization" != "RPKM" ]] && [[ "$normalization" != "CPM" ]] && [[ "$normalization" != "BPM" ]] && [[ "$normalization" != "RPGC" ]]; then
    echo "Error: Incorrect method to normalize the number of reads per bin for the bamCoverage command of deeptools."
    exit 1
  fi
  mkdir -p $bam_output_dir/TSSC
  
  output_file=$bam_output_dir/TSSC/TSS.bed
  

  if [[ "$genomeVer" == "hg19" ]]; then
    chrom_size="$script_dir/DANPOS3/hg19.chrom.sizes"
    if [[ -z "$TSSfile" ]]; then
      TSSfile=$script_dir/TSScoverage/TSS_hg19_uniq.bed
    fi
  elif [[ "$genomeVer" == "hg38" ]]; then
    chrom_size="$script_dir/DANPOS3/hg38.chrom.sizes"
    if [[ -z "$TSSfile" ]]; then
      TSSfile=$script_dir/TSScoverage/TSS_hg38_uniq.bed
    fi
  fi


  # Create or clear the output file
  echo -n > "$output_file"

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
  # correct the output bed file
  awk 'BEGIN {OFS="\t"} NR==FNR {sizes[$1]=$2; next} {if ($3 > sizes[$1]) $3=sizes[$1]; print}' $chrom_size $bam_output_dir/TSSC/TSS.bed | \
  awk 'NF >= 6' | awk -v OFS='\t' '{$1=$1; print}' > $bam_output_dir/TSSC/correct_6col_TSS.bed

  if [[ "$sequencing" == "pair" ]];then
    bamCoverage -b $inputBam -o $bam_output_dir/TSSC/coverage.bw --extendReads $addbamCoverage --numberOfProcessors $threads --normalizeUsing $normalization --minFragmentLength 180 --maxFragmentLength 200
  elif [[ "$sequencing" == "single" ]];then
    bamCoverage -b $inputBam -o $bam_output_dir/TSSC/coverage.bw --extendReads 167 $addbamCoverage --numberOfProcessors $threads --normalizeUsing $normalization --minFragmentLength 180 --maxFragmentLength 200
  fi
  multiBigwigSummary BED-file -b $bam_output_dir/TSSC/coverage.bw -o $bam_output_dir/TSSC/coverage.npz --BED $bam_output_dir/TSSC/correct_6col_TSS.bed --outRawCounts $bam_output_dir/TSSC/coverage.txt --numberOfProcessors $threads $addmultiBigwigSummary
  sed '1d' $bam_output_dir/TSSC/coverage.txt | sort -k1,1V -k2,2n -k3,3n | awk 'BEGIN {print "chr\tstart\tend\tcoverage"} {print}' > $bam_output_dir/TSSC/average_coverage.txt
  rm $bam_output_dir/TSSC/coverage.txt
  rm $bam_output_dir/TSSC/coverage.bw
  rm $bam_output_dir/TSSC/coverage.npz
  rm $bam_output_dir/TSSC/TSS.bed
  rm $bam_output_dir/TSSC/correct_6col_TSS.bed
}

# Run the analysis
run_analysis

if [[ -z "$labelFile" ]];then
    echo "Error: Please provide the label information file, seperated by comma.."
    exit 1
else
    labelFile=$(realpath "$labelFile")
fi
cp $labelFile $outputDir/Samples/label.txt

mkdir -p $outputDir/Features
python $script_dir/Data_transformation.py --config $outputDir/Samples/config_feature.json --root $outputDir/Samples --output $outputDir/Features

if [[ ",$Feature," == *",PFE,"* ]]; then
   Rscript $script_dir/Epic-seq/code/standard_PFE.R $outputDir/Features/PFE.csv $outputDir/Features/PFE_standard.csv
   rm $outputDir/Features/PFE.csv
   mv $outputDir/Features/PFE_standard.csv $outputDir/Features/PFE.csv
fi

if [[ ",$Feature," == *",NP,"* ]]; then
   mkdir -p $outputDir/Features/NP_site_list
   Rscript $script_dir/Griffin/merge_NP_sites.R $outputDir/Samples $outputDir/Features/NP_site_list  
fi

rm -r $outputDir/Samples

# downstream analysis
if [[ "$noDA" == 0 ]];then
  echo '
  *****Feature Processing is in execution*****'
  mkdir -p $outputDir/Feature_Processing_and_Selection

  # Feature_Processing
  mkdir -p $outputDir/Feature_Processing_and_Selection/Feature_Processing

  # Missing value filter
  mkdir -p $outputDir/Feature_Processing_and_Selection/Feature_Processing/filter
  python $script_dir/Feature_Processing/Missing_data_removal.py \
  --input $outputDir/Features \
  --output $outputDir/Feature_Processing_and_Selection/Feature_Processing/filter

  # standardization
  mkdir -p $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization
  python $script_dir/Feature_Processing/Standardization.py \
  --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing/filter \
  --output_dir  $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization \
  --method $standardM
  
  Rscript $script_dir/Feature_Processing/merge_features.R $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization $outputDir/Feature_Processing_and_Selection/Feature_Processing $Feature
  if [[ ",$Feature," == *",CNA,"* ]]; then 
    cp  $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization/CNA.csv  $outputDir/Feature_Processing_and_Selection/Feature_Processing/CNA.csv
  fi

  if [[ ",$Feature," == *",FP,"* ]]; then 
    cp  $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization/FP_fragmentation_profile.csv  $outputDir/Feature_Processing_and_Selection/Feature_Processing/FP.csv 
  fi
  
  if [[ ",$Feature," == *",OCF,"* ]]; then 
    cp  $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization/OCF.csv  $outputDir/Feature_Processing_and_Selection/Feature_Processing/OCF.csv 
  fi
  
  if [[ ",$Feature," == *",FPR,"* ]]; then 
    cp  $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization/FPR_fragmentation_profile_regions.csv  $outputDir/Feature_Processing_and_Selection/Feature_Processing/FPR.csv 
  fi
  
  if [[ ",$Feature," == *",PFE,"* ]]; then 
    cp  $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization/PFE.csv  $outputDir/Feature_Processing_and_Selection/Feature_Processing/PFE.csv 
  fi
  
  if [[ ",$Feature," == *",TSSC,"* ]]; then 
    cp  $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization/TSSC_average_coverage.csv  $outputDir/Feature_Processing_and_Selection/Feature_Processing/TSSC.csv 
  fi
      
    rm -r $outputDir/Feature_Processing_and_Selection/Feature_Processing/filter
    rm -r $outputDir/Feature_Processing_and_Selection/Feature_Processing/standardization

  if [[ "$noML" == 1 ]];then
    # feature_selection
    echo '
    *****Feature Selection is in execution*****'
    mkdir -p $outputDir/Feature_Processing_and_Selection/Feature_Selection/

    # filter method
    python $script_dir/Feature_Selection/filter_methods.py \
    --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
    --output_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection/ \
    --methods $filterMethod --percentages $filterNum

    # wrapper method
    python $script_dir/Feature_Selection/wrapper_methods.py \
    --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
    --output_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection/ \
    --methods $wrapperMethod --percentage $wrapperNum

    # embeded method
    python $script_dir/Feature_Selection/embedded_methods.py \
    --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
    --output_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection/ \
    --methods $embeddedMethod --percentages $embeddedNum

    # hybrid method
    if [[ "$hybridType" != "FW" ]] && [[ "$hybridType" != "FE" ]];then
      echo "Error: Wrong hybrid method."
      exit 1
    fi

    if [[ -z "$hybridMethod1" ]] ;then
      hybridMethod1="IG"
    fi

    mkdir -p $outputDir/Feature_Processing_and_Selection/Feature_Selection/hybrid_filter_select
    python $script_dir/Feature_Selection/filter_methods.py \
    --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
    --output_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection/hybrid_filter_select \
    --methods $hybridMethod1 --percentages $hybridNum1
    
    if [[ "$hybridType" = "FE" ]] ;then
      if [[ -z "$hybridMethod2" ]] ;then
        hybridMethod2="LASSO"
      fi
      python $script_dir/Feature_Selection/embedded_methods.py \
      --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection/hybrid_filter_select \
      --output_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection \
      --methods LASSO --percentages $hybridNum2      
    elif [[ "$hybridType" = "FW" ]]; then
      if [[ -z "$hybridMethod2" ]] ;then
        hybridMethod2="BOR"
      fi
      python $script_dir/Feature_Selection/wrapper_methods.py \
      --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection/hybrid_filter_select \
      --output_dir $outputDir/Feature_Processing_and_Selection/Feature_Selection/ \
      --methods $hybridMethod2 --percentage $hybridNum2
    fi
    
    rm -r $outputDir/Feature_Processing_and_Selection/Feature_Selection/hybrid_filter_select

  elif [[ "$noML" == 0  ]]; then
    # Machine Learning Model Building
    echo '
    *****Machine Learning Model Building is in execution*****'
    
      # hybrid method
    if [[ "$hybridType" != "FW" ]] && [[ "$hybridType" != "FE" ]];then
      echo "Error: Wrong hybrid method."
      exit 1
    fi

    if [[ "$hybridType" = "FE" ]] ;then
      if [[ -z "$hybridMethod1" ]] ;then
        hybridMethod1="IG"
      fi
      if [[ -z "$hybridMethod2" ]] ;then
        hybridMethod2="LASSO"
      fi
    elif [[ "$hybridType" = "FW" ]]; then
      if [[ -z "$hybridMethod1" ]] ;then
        hybridMethod1="IG"
      fi
      if [[ -z "$hybridMethod2" ]] ;then
        hybridMethod2="BOR"
      fi
    fi


    mkdir -p $outputDir/Machine_Learning
    if [[ "$classNum" -lt 2 ]];then
      echo "Error: Wrong classification category."
      exit 1
    fi

    # Single modality
    mkdir -p $outputDir/Machine_Learning/single_modality

    if [[ "$cvSingle" == "LOO" ]];then
      python $script_dir/Machine_learning/run_single_modality.py \
      --modality single \
      --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
      --DA_output_dir $outputDir/Machine_Learning/single_modality \
      --classNum $classNum \
      --cvSingle LOO \
      --classifierSingle $classifierSingle \
      --filterMethod $filterMethod \
      --filterFrac $filterNum \
      --wrapperMethod $wrapperMethod \
      --wrapperFrac $wrapperNum \
      --embeddedMethod $embeddedMethod \
      --embeddedFrac $embeddedNum \
      --hybridType $hybridType \
      --hybridMethod1 $hybridMethod1 \
      --hybridFrac1 $hybridNum1 \
      --hybridMethod2 $hybridMethod2 \
      --hybridFrac2 $hybridNum2
    elif [[ "$cvSingle" == "KFold" ]];then
      python $script_dir/Machine_learning/run_single_modality.py \
      --modality single \
      --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
      --DA_output_dir $outputDir/Machine_Learning/single_modality \
      --classNum $classNum \
      --cvSingle KFold \
      --nsplitSingle $nsplitSingle \
      --classifierSingle $classifierSingle \
      --filterMethod $filterMethod \
      --filterFrac $filterNum \
      --wrapperMethod $wrapperMethod \
      --wrapperFrac $wrapperNum \
      --embeddedMethod $embeddedMethod \
      --embeddedFrac $embeddedNum \
      --hybridType $hybridType \
      --hybridMethod1 $hybridMethod1 \
      --hybridFrac1 $hybridNum1 \
      --hybridMethod2 $hybridMethod2 \
      --hybridFrac2 $hybridNum2
    fi

    # Multiple modalities
    mkdir -p $outputDir/Machine_Learning/multiple_modality/Concatenation_based
    mkdir -p $outputDir/Machine_Learning/multiple_modality/Model_based
    mkdir -p $outputDir/Machine_Learning/multiple_modality/Transformation_based

    if [[ "$cvMulti" == "LOO" ]];then
      python $script_dir/Machine_learning/run_multi_modality.py \
      --modality multi \
      --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
      --DA_output_dir $outputDir/Machine_Learning/multiple_modality \
      --fusion_type concat model trans \
      --model_method $modelMethod \
      --trans_method $transMethod \
      --classNum $classNum \
      --cvMulti LOO \
      --classifierMulti $classifierMulti \
      --filterMethod $filterMethod \
      --filterFrac $filterNum \
      --wrapperMethod $wrapperMethod \
      --wrapperFrac $wrapperNum \
      --embeddedMethod $embeddedMethod \
      --embeddedFrac $embeddedNum \
      --hybridType $hybridType \
      --hybridMethod1 $hybridMethod1 \
      --hybridFrac1 $hybridNum1 \
      --hybridMethod2 $hybridMethod2 \
      --hybridFrac2 $hybridNum2
    elif [[ "$cvMulti" == "KFold" ]];then
      python $script_dir/Machine_learning/run_multi_modality.py \
      --modality multi \
      --input_dir $outputDir/Feature_Processing_and_Selection/Feature_Processing \
      --DA_output_dir $outputDir/Machine_Learning/multiple_modality \
      --fusion_type concat model trans \
      --model_method $modelMethod \
      --trans_method $transMethod \
      --classNum $classNum \
      --cvMulti KFold \
      --nsplitMulti $nsplitMulti \
      --classifierMulti $classifierMulti \
      --filterMethod $filterMethod \
      --filterFrac $filterNum \
      --wrapperMethod $wrapperMethod \
      --wrapperFrac $wrapperNum \
      --embeddedMethod $embeddedMethod \
      --embeddedFrac $embeddedNum \
      --hybridType $hybridType \
      --hybridMethod1 $hybridMethod1 \
      --hybridFrac1 $hybridNum1 \
      --hybridMethod2 $hybridMethod2 \
      --hybridFrac2 $hybridNum2
    fi

    first_file=true
    find "$outputDir/Machine_Learning/multiple_modality/concat" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
        input_file_1="$dir/all_methods_all_classifiers_metrics.csv"
        input_file_2="$dir/all_methods_all_classifiers_predictions.csv"
        if [[ -f "$input_file_1" ]] && [[ -f "$input_file_2" ]]; then
            if [[ "$first_file" == true ]]; then
                head -n 1 "$input_file_1" > "$outputDir/Machine_Learning/multiple_modality/Concatenation_based/multiple_modality_metrics.csv"
                head -n 1 "$input_file_2" > "$outputDir/Machine_Learning/multiple_modality/Concatenation_based/multiple_modality_probabilities.csv"
                first_file=false
            fi
            tail -n +2 "$input_file_1" >> "$outputDir/Machine_Learning/multiple_modality/Concatenation_based/multiple_modality_metrics.csv"
            tail -n +2 "$input_file_2" >> "$outputDir/Machine_Learning/multiple_modality/Concatenation_based/multiple_modality_probabilities.csv"
        fi
    done

    first_file=true
    find "$outputDir/Machine_Learning/multiple_modality/model" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
        input_file_1="$dir/all_methods_all_classifiers_metrics.csv"
        input_file_2="$dir/all_methods_all_classifiers_predictions.csv"
        if [[ -f "$input_file_1" ]] && [[ -f "$input_file_2" ]]; then
            if [[ "$first_file" == true ]]; then
                head -n 1 "$input_file_1" > "$outputDir/Machine_Learning/multiple_modality/Model_based/multiple_modality_metrics.csv"
                head -n 1 "$input_file_2" > "$outputDir/Machine_Learning/multiple_modality/Model_based/multiple_modality_probabilities.csv"
                first_file=false
            fi
            tail -n +2 "$input_file_1" >> "$outputDir/Machine_Learning/multiple_modality/Model_based/multiple_modality_metrics.csv"
            tail -n +2 "$input_file_2" >> "$outputDir/Machine_Learning/multiple_modality/Model_based/multiple_modality_probabilities.csv"
        fi
    done

    first_file=true
    find "$outputDir/Machine_Learning/multiple_modality/trans" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
        input_file_1="$dir/all_methods_all_classifiers_metrics.csv"
        input_file_2="$dir/all_methods_all_classifiers_predictions.csv"
        if [[ -f "$input_file_1" ]] && [[ -f "$input_file_2" ]]; then
            if [[ "$first_file" == true ]]; then
                head -n 1 "$input_file_1" > "$outputDir/Machine_Learning/multiple_modality/Transformation_based/multiple_modality_metrics.csv"
                head -n 1 "$input_file_2" > "$outputDir/Machine_Learning/multiple_modality/Transformation_based/multiple_modality_probabilities.csv"
                first_file=false
            fi
            tail -n +2 "$input_file_1" >> "$outputDir/Machine_Learning/multiple_modality/Transformation_based/multiple_modality_metrics.csv"
            tail -n +2 "$input_file_2" >> "$outputDir/Machine_Learning/multiple_modality/Transformation_based/multiple_modality_probabilities.csv"
        fi
    done

    rm -r "$outputDir/Machine_Learning/multiple_modality/concat"
    rm -r "$outputDir/Machine_Learning/multiple_modality/model"
    rm -r "$outputDir/Machine_Learning/multiple_modality/trans"

  fi
fi