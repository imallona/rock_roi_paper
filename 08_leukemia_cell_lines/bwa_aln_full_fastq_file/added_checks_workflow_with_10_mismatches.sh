#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##


## running bwa aln with standard settings on full fusion file

mkdir -p $WD/test/genome
cd $WD/test/genome

# transcript reference generation for bwa aln. Added -q to wget to suppress the progress bar.

wget -q http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$TRANSCRIPTOME".gz

pigz --decompress *gz

# ENSG00000097007: ABL1
# ENSG00000186716: BCR

# extracting all transcripts for ABL1 and BCR from the .fa transcriptome

grep "ENSG00000097007" "$TRANSCRIPTOME"| sed 's/>//g' > ABL1_human.txt
grep "ENSG00000186716" "$TRANSCRIPTOME"| sed 's/>//g' >  BCR_human.txt

for i in $(cat ABL1_human.txt); do
    samtools faidx $TRANSCRIPTOME "$i" >> ABL1_human.fa
done

for i in $(cat BCR_human.txt) ; do
    samtools faidx $TRANSCRIPTOME "$i"  >> BCR_human.fa
done

cat "$CUSTOM_FA" ABL1_human.fa BCR_human.fa > combined.fa

for file in "$CUSTOM_FA" ABL1_human.fa BCR_human.fa; do
  cat "$file"
  echo
done > combined.fa

# index reference

# bwa version: bwa 0.7.18

bwa index -p indexed combined.fa

# running bwa aln

mkdir -p $WD/test/output
cd $WD/test/output

# -k: Maximum edit distance in the seed [2]
# -n: Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]

bwa aln $WD/test/genome/indexed -t $NTHREADS $cdna_cell_line -n 10 -k 4 > bwa_aln_alignments.sai

bwa samse -f bwa_aln_alignments.sam $WD/test/genome/indexed bwa_aln_alignments.sai $cdna_cell_line

samtools view -o out.bam bwa_aln_alignments.sam

samtools view -H out.bam > header.txt

samtools sort -@ $NTHREADS out.bam -o output.sorted.bam

samtools index -@ $NTHREADS output.sorted.bam

rm bwa_aln_alignments.sam bwa_aln_alignments.sai

# subsetting .bam file for unique alignments and only mapped alignments

samtools view -@ $NTHREADS -b -F 4 output.sorted.bam > mapped.bam

samtools view -@ $NTHREADS -H mapped.bam > header.txt

samtools view -@ $NTHREADS mapped.bam | grep -v "XA:" | cat header.txt - | samtools view -Sbh > no_xa_bwa_aln_sorted.bam

# checking length of the .bam file to know number of entries

echo 'Number of aligments without XA tag'

samtools view no_xa_bwa_aln_sorted.bam| wc -l

rm header.txt

# extracting read ids from bwa aln .bam file

samtools view -@ $NTHREADS no_xa_bwa_aln_sorted.bam | cut -f1 > read_ids.txt

# generating subsetted starsolo .bam file with reads found from bwa aln

samtools view -H -@ $NTHREADS $starsolo_bam > starsolo_header.txt

samtools view -@ $NTHREADS $starsolo_bam | grep -w -f read_ids.txt | cat starsolo_header.txt - | samtools view -Sbh > subsetted_starsolo.bam

# checking number of entries in subsetted_starsolo.bam

echo 'Number of entries extracted by STARsolo .bam file'

samtools view subsetted_starsolo.bam | wc -l

# more entries --> multialigner reads. They are removed in a later step

rm starsolo_header.txt

# extracting read id, CB and UMI from the subsetted_starsolo.bam. Needs to be done separately for mapped and unmapped and first need to separate between mapped and unmapped otherwise for the unmapped getting additional fields at positions 22 and 23 for the mapped

echo 'Mapped'

samtools view -@ $NTHREADS -F 4 subsetted_starsolo.bam | cut -f1,27,28 >> starsolo_reads_with_cb_ub.txt

echo 'Unmapped'

samtools view -@ $NTHREADS -f 4 subsetted_starsolo.bam | cut -f1,22,23 >> starsolo_reads_with_cb_ub.txt

# sort the read ids and unique them 

sort -k 1,1 starsolo_reads_with_cb_ub.txt -u > sorted_starsolo_reads_with_cb_ub.txt

# checking that now they have the same read length as the bwa aln ones

echo 'Number of unique entries extracted by STARsolo .bam file'

wc -l sorted_starsolo_reads_with_cb_ub.txt

# extracting read id, detected fusion/non-fused and read position

samtools view -@ $NTHREADS no_xa_bwa_aln.sorted.bam | cut -f1,3,4 > bwa_aln_reads.txt

sort -k 1,1 bwa_aln_reads.txt > sorted_bwa_aln_reads.txt 

# generate file with read id, detected fusion/non-fused, CB and UMI (to be imported into R)

join -1 1 -2 1 -t $'\t' sorted_bwa_aln_reads.txt sorted_starsolo_reads_with_cb_ub.txt > bwa_aln_annotated_fusion.txt

# checking length of file

echo "Length of joint file"

wc -l bwa_aln_annotated_fusion.txt

rm sorted_starsolo_reads_with_cb_ub.txt sorted_bwa_aln_reads.txt starsolo_reads_with_cb_ub.txt bwa_aln_reads.txt read_ids.txt

# to make file smaller only keep alignments with CB and UB

grep -v "CB:Z:-" bwa_aln_annotated_fusion.txt | grep -v "UB:Z:-" > bwa_aln_with_cb_ub_annotated_fusion.txt

rm bwa_aln_annotated_fusion.txt

# spot checks on barcodes: extracting read id and then going back to the starsolo .bam file to check if they are the same. Picking line 100

echo 'Spot check: line 100 of newly generated file'

awk 'NR==100' bwa_aln_with_cb_ub_annotated_fusion.txt
awk 'NR==100' bwa_aln_with_cb_ub_annotated_fusion.txt | cut -f1 > read_id_100.txt

echo 'Spot check: line from .bam file with same read id as line 100 of newly generated file'

samtools view -@ $NTHREADS $starsolo_bam | grep -w -f read_id_100.txt
