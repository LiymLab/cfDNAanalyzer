#!/bin/bash

usage(){
  cat << EOF
This is a script for processing raw sequencing data (FASTQ files) into compatible formats (BAM files).
Please ensure that the following software is installed: Trim-galore, Bowtie2, BWA, SAMtools.

Usage:
   bash Raw_data_processing.sh -i <InputFile> -o <OutputDirectory> -t <threads> [Options]
 
-----General options -----  
  -I  FILE      A text file containing all input FASTQ files with one FASTQ file per line.
                FASTQ files must end with "fastq.gz".
                IF input is pair-end FASTQ files, files should have the same prefix like 'sample_1.fastq.gz, sample_2.fastq.gz'
  -o  DIR       Output directory for all the results. Default: [./]
  -s  STR       Sequencing method of input BAM files (single/pair). Default: [pair]
  -t  INT       Number of threads to use. Default: [1]
  
----- Options specific Trim-galore -----
  -q  INT       Trim low-quality ends from reads in addition to adapter removal. Default Phred score: [20].
  -l  INT       Discard reads that became shorter than length INT (bp) because of either quality or adapter trimming. A value of '0' effectively disables this behaviour. Default: [20].

----- Options specific for aligning -----
  -a  STR       Aligning software to use, Bowtie2 or BWA. Default: [BWA].
  -r  FILE      Reference genome in FASTA format to bulid index for Bowtie2 or BWA.
  
----- Options specific for SAMtools -----
  -F  INT       only include reads with none of the FLAGS in INT present. Default: [1796].
  -f  INT       only include reads with all of the FLAGs in INT present. Default: [3].
  -Q  INT       only include reads with mapping quality >= INT. Default: [20].

EOF
    exit 0
}

if [ $# -eq 0 ] || ([ $# -eq 1 ] && { [ "$1" == "-h" ] || [ "$1" == "--help" ]; }); then
    usage
fi

# Default
inputList=""
output_dir=$(pwd)
threads=1
sequencing=pair
quality=20
fasta=""
length=20
map_soft=BWA
Flag=1796
flag=3
Quality=20

TEMP=$(getopt -o hI:o:t:s:q:r:l:a:F:f:Q: \
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
        -I)                       inputList=$2; shift 2;;
        -o)                       output_dir=$2; shift 2;;
        -t)                       threads=$2; shift 2;;
        -s)                       sequencing=$2; shift 2;;
        -q)                       quality=$2; shift 2;;
        -r)                       fasta=$2; shift 2;;
        -l)                       length=$2; shift 2;;
        -a)                       map_soft=$2; shift 2;;
        -F)                       Flag=$2; shift 2;;
        -f)                       flag=$2; shift 2;;
        -Q)                       Quality=$2; shift 2;;
        --)                       shift; break;;
        *)                        echo "Unexpected option: $1"; exit 1;;
    esac
done

if [ -z "$inputList" ]; then
  echo "Error: Input file list is not provided."
  exit 1
fi

if [ -z "$fasta" ]; then
  echo "Error: Input reference fasta file list is not provided."
  exit 1
fi

inputList=$(realpath "$inputList")
fasta=$(realpath "$fasta")
output_dir=$(realpath "$output_dir")

mkdir -p $output_dir
mkdir $output_dir/index
if [[ "$map_soft" == "BWA" ]]; then
    bwa index $fasta -p $output_dir/index/index
elif [[ "$map_soft" == "Bowtie2" ]]; then
    bowtie2-build $fasta $output_dir/index/index
fi


if [[ "$sequencing" == "single" ]]; then
    for input in $(cat "$inputList"); do
      if [ ! -f "$input" ]; then
        echo "Error: Input file was not found: $input"
        exit 1
      fi
      input=$(realpath "$input")
      filename=$(basename "$input" .fastq.gz)
      
      mkdir $output_dir/trim_galore
      trim_galore $input \
                --quality $quality \
                --length $length \
                -j $threads \
                -o $output_dir/trim_galore

      if [[ "$map_soft" == "BWA" ]]; then
        bwa mem -t $threads \
              $output_dir/index/index \
              $output_dir/trim_galore/${filename}_trimmed.fq.gz | samtools view -S -bF 4 - > $output_dir/${filename}.bam
      elif [[ "$map_soft" == "Bowtie2" ]]; then
        bowtie2 -p $threads \
              -x $output_dir/index/index \
              -U $output_dir/trim_galore/${filename}_trimmed.fq.gz | samtools view -S -bF 4 - > $output_dir/${filename}.bam
      fi
      
      samtools view $output_dir/${filename}.bam \
                  -h -b \
                  -@ $threads \
                  -q $Quality \
                  -F $Flag \
                  -f $flag  > $output_dir/filtered_${filename}.bam
    done
    
    rm -r $output_dir/trim_galore
    rm $output_dir/${filename}.bam
    
elif [[ "$sequencing" == "pair" ]]; then
    grep "_1\.fastq\.gz$" $inputList | while read -r filename_1; do
      if [ ! -f "$filename_1" ]; then
        echo "Error: Input file was not found: $filename_1"
        exit 1
      fi
      filename_2=$(echo $filename_1 | sed 's/1.fastq.gz/2.fastq.gz/')
      if [ ! -f "$filename_2" ]; then
        echo "Error: Could not find mate file for: $filename_1"
        exit 1
      fi
      
      filename_1=$(realpath "$filename_1")
      filename_2=$(realpath "$filename_2")
      prefix=$(basename "$filename_1" _1.fastq.gz)
      
      mkdir $output_dir/trim_galore
      trim_galore $filename_1 $filename_2 \
                  --paired \
                  --quality $quality \
                  --length $length \
                  -j $threads \
                  -o $output_dir/trim_galore
                  
      if [[ "$map_soft" == "BWA" ]]; then
        bwa mem -t $threads \
                $output_dir/index/index \
                $output_dir/trim_galore/${prefix}_1_val_1.fq.gz $output_dir/trim_galore/${prefix}_2_val_2.fq.gz \
                 | samtools view -S -bF 4 - > $output_dir/${prefix}.bam
      elif [[ "$map_soft" == "Bowtie2" ]]; then
        bowtie2 -1 $output_dir/trim_galore/${prefix}_1_val_1.fq.gz -2 $output_dir/trim_galore/${prefix}_2_val_2.fq.gz \
                -x $output_dir/index/index \
                -p $threads | samtools view -S -bF 4 - > $output_dir/${prefix}.bam
      fi
      
      samtools view $output_dir/${prefix}.bam \
                  -h -b \
                  -@ $threads \
                  -q $Quality \
                  -F $Flag \
                  -f $flag  > $output_dir/filtered_${prefix}.bam
    done

    rm -r $output_dir/trim_galore
    rm $output_dir/${prefix}.bam
fi

rm -r $output_dir/index

