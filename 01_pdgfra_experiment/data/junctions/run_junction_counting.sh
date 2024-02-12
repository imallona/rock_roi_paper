#!/bin/bash

if [ "$#" -ne 4 ]
then
    echo "This script subsets BAM to specific region, splits by CB, runs featureCounts for junctions on each split"
    echo "'tag' will be prepended to the dir name of featureCounts outputs (can be used to label sample)"
    echo "Run as: $0 input-bam region-bed input-gtf tag"
    exit 1
fi

#echo "Software"
#echo "----------"
#samtools --version
#bamtools --version
#featureCounts -v
#echo "----------"

tmpdir="$4___$(date +"%Y-%m-%d_%T")"
mkdir -p $tmpdir

bamf=$1
regionbed=$2
gtff=$3

subbam="$tmpdir/subset.bam"

# subset to region of interest
echo "Subsetting $bamf according to $regionbed -> $subbam .."
samtools view -b -L $regionbed -h $bamf > $subbam
echo "done."

# split subset-BAM by CB
echo "Splitting $subbam by CB .."
bamtools split -in $subbam -stub "$tmpdir/split_" -tag CB
echo "done. (Got $(ls $tmpdir/split_* | wc -w) files.)"

## run featureCounts on each
echo "Running featureCounts"
for F in $(ls $tmpdir/split_*)
do
    echo "- working on $F"
    featureCounts -a $gtff -o $F.featurecounts $F -F GTF -t exon -g transcript_id -f -O -M -T 2 --fraction
done
echo "done."

