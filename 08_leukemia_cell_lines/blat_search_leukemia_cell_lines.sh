#!/bin/bash
##
## Giulia Moro
##
## BLAT search for fusion cell lines
##

## Working directory: /home/gmoro/blat_leukemia_cell_lines in sherborne

## Installation of BLAT

#conda install bioconda::blat

## Generating .bam file based on starsolo .bam file with alignments which are not detected by bwa aln but only by grep

samtools view -H -@ 30 ~/cell_lines_leukemia/rock_roi_repo_for_cell_lines/output/align_tso/leukemia_cell_lines/Aligned.sortedByCoord.out.bam > header.txt

samtools view -@ 30 ~/cell_lines_leukemia/rock_roi_repo_for_cell_lines/output/align_tso/leukemia_cell_lines/Aligned.sortedByCoord.out.bam | grep -f n_10_loose_mapping_bwa_aln_only | cat header.txt - | samtools view -@ 30 -Sbh > n_10_grep_only.bam

rm header.txt

## Generating .fa file from .bam file

samtools fasta -@ 30 n_10_grep_only.bam > n_10_grep_only.fa

## Running BLAT, based on same genome we used for STARsolo mapping 

human_fa=GRCh38.p13.genome.fa
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_fa".gz
pigz --decompress GRCh38.p13.genome.fa.gz

blat GRCh38.p13.genome.fa n_10_grep_only.fa leukemia_cell_lines_grep_only.psl
