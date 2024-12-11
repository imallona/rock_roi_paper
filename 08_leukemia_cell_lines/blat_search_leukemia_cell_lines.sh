#!/bin/bash
##
## Giulia Moro
##
## BLAT search for fusion cell lines 
##

## Installation of BLAT

#conda install bioconda::blat

## Generating .bam file with alignments which are not detected by bwa aln but only by grep

samtools view -@ 30 ~/cell_lines_leukemia/rock_roi_repo_for_cell_lines/output/align_tso/leukemia_cell_lines/Aligned.sortedByCoord.out.bam | grep -f n_10_loose_mapping_bwa_aln_only | cat header.txt - | samtools view -@ 30 -Sbh > n_10_grep_only.bam

## Converting .bam file to .fastq 

