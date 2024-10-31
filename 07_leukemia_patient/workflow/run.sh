#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH #path to STAR for STARsolo and STAR fusion
export WD=~/test_leukemia_simulated_reads #working directory
export COMBINED_INDEXED_GENOME=~/mapping_leukemia/data/index #STAR indexed genome 
export CUSTOM_FA=~/leukemia_bwamem2/genome/BCR_ABL.fa #file containing fusion cDNAs
export NTHREADS=30 #number of threads
export TRANSCRIPTOME=gencode.v38.pc_transcripts.fa #transcriptome to download
export STARSOLO_BAM=$WD/starsolo/Aligned.sortedByCoord.out.bam #path to .bam file, generated during script

./workflow.sh
