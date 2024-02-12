#!/bin/bash

if [ "$#" -ne 3 ]
then
    echo "This script subsets BAM to CBs specified in a file"
    echo "Run as: $0 input-bam barcode-file output-file"
    exit 1
fi

echo "Software"
echo "----------"
samtools --version
echo "----------"

tmpdir="$3___$(date +"%Y-%m-%d_%T")"
mkdir -p $tmpdir

bamf=$1
bcfile=$2
outbam=$3

# Save the header lines
echo "Saving SAM header .. $subbam .."
samtools view -H $bamf > $tmpdir/SAM_header
echo "Extracting reads .."
samtools view $bamf | LC_ALL=C grep -F -f $bcfile > $tmpdir/filtered_SAM_body
echo "Reconstructing BAM file .."
cat $tmpdir/SAM_header $tmpdir/filtered_SAM_body > $tmpdir/filtered.sam
samtools view -b $tmpdir/filtered.sam > $outbam
samtools index $outbam

# remove not-needed files
#/bin/rm SAM_header filtered.sam filtered_SAM_body

