#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <BigWig file> <BED file> <Output file name> <Region>"
    exit 1
fi

BIGWIG_FILE=$1
BED_FILE=$2
OUTPUT_BASE=$3
REGION=$4

OUTPUT_FILE="${OUTPUT_BASE}_${REGION}.txt"

BIGWIGAVERAGEOVERBED="./TSScoverage/bigWigAverageOverBed"

$BIGWIGAVERAGEOVERBED $BIGWIG_FILE $BED_FILE $OUTPUT_FILE

echo "$REGION coverage calculation completed: $OUTPUT_FILE"
