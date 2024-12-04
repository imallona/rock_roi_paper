#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##


## running bwa aln with standard settings on full fusion file

mkdir -p $WD/bwa_aln/genome
cd $WD/bwa_aln/genome

# transcript reference generation for bwa aln

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$TRANSCRIPTOME".gz

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

# index reference

# bwa version: bwa 0.7.18

bwa index -p indexed combined.fa

mkdir -p $WD/bwa_aln/output
cd $WD/bwa_aln/output

# -d: Disallow a long deletion within INT bp towards the 3’-end
# -i: Disallow an indel within INT bp towards the ends
# -k: Maximum edit distance in the seed
# -O: gap open penalty

bwa aln $WD/bwa_aln/genome/indexed -t $NTHREADS $cdna_cell_line > bwa_aln_alignments.sai

bwa samse -f bwa_aln_alignments.sam $WD/bwa_aln/genome/indexed bwa_aln_alignments.sai $cdna_cell_line

samtools view -o out.bam bwa_aln_alignments.sam

samtools view -H out.bam > header.txt

samtools sort out.bam -o output.sorted.bam

rm bwa_aln_alignments.sam bwa_aln_alignments.sai
