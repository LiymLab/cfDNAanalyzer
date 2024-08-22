#!/bin/bash

# To use: ./calculate_wps.sh <BAMFILE> <OUTPUT> <REGION> <EXTREGION>

BAMFILE=$1
OUTPUT=$2
REGION=$3
EXTREGION=$4
Minlength_long=$5 
Maxlength_long=$6
window_long=$7 
Minlength_short=$8
Maxlength_short=$9
window_short=${10}
script_dir=${11}

# if [ "$#" -ne 11 ]; then
#     echo "Usage: $0 <BAMFILE> <OUTPUT> <REGION> <EXTREGION> <Minlength_long> <Maxlength_long> <window_long> <Minlength_short> <Maxlength_short> <window_short> "
#     exit 1
# fi

# Running the extraction of coverage, read starts and WPS for the short fraction
# echo "Extracting WPS for long fragments..."
$script_dir/WPS/samtools view -u -m $Minlength_long -M $Maxlength_long $BAMFILE $EXTREGION | $script_dir/WPS/FilterUniqueBAM.py -p | $script_dir/WPS/extractReadStartsFromBAM2Wig.py -p -r $REGION -w $window_long -c $OUTPUT/COVERAGE_long.wig.gz -n $OUTPUT/WPS_long.wig.gz -s $OUTPUT/STARTS_long.wig.gz 

# Running the extraction of coverage, read starts and WPS for the long fraction
# echo "Extracting WPS for short fragments..."
$script_dir/WPS/samtools view -u -m $Minlength_short -M $Maxlength_short $BAMFILE $EXTREGION | $script_dir/WPS/FilterUniqueBAM.py -p | $script_dir/WPS/extractReadStartsFromBAM2Wig.py -p -r $REGION -w $window_short -c $OUTPUT/COVERAGE_short.wig.gz -n $OUTPUT/WPS_short.wig.gz -s $OUTPUT/STARTS_short.wig.gz

# echo "WPS extraction completed."
