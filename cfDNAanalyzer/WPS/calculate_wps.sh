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

$script_dir/WPS/samtools view -u -m $Minlength_long -M $Maxlength_long $BAMFILE $EXTREGION | $script_dir/WPS/FilterUniqueBAM.py -p| $script_dir/WPS/extractReadStartsFromBAM2Wig.py -p -r $REGION -w $window_long -n $OUTPUT/WPS_long.wig.gz -c $OUTPUT/coverage_long.wig.gz -s $OUTPUT/starts_long.wig.gz

$script_dir/WPS/samtools view -u -m $Minlength_short -M $Maxlength_short $BAMFILE $EXTREGION | $script_dir/WPS/FilterUniqueBAM.py -p | $script_dir/WPS/extractReadStartsFromBAM2Wig.py -p -r $REGION -w $window_short -n $OUTPUT/WPS_short.wig.gz -c $OUTPUT/coverage_short.wig.gz -s $OUTPUT/starts_short.wig.gz

