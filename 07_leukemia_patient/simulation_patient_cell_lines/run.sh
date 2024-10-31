#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH #path to STAR for STARsolo and STAR fusion
export PATH=/home/gmoro/bwa/:$PATH #path to bwa
export WD=~/test_leukemia_simulated_reads #working directory
export COMBINED_INDEXED_GENOME=~/mapping_leukemia/data/index #STAR indexed genome 
export COMBINED_GTF_GENOME=~/mapping_leukemia_data/genome/combined.gtf #genome .gtf
export r1=/home/gmoro/simulated_leukemia_data/combined_r1.fastq.gz #simulated reads
export r2=/home/gmoro/simulated_leukemia_data/combined_r2.fastq.gz #simulated reads
export CL1=/home/gmoro/whitelist_96x3/BD_CLS1.txt #first whitelist for cell label 1
export CL2=/home/gmoro/whitelist_96x3/BD_CLS2.txt #second whitelist for cell label 2
export CL3=/home/gmoro/whitelist_96x3/BD_CLS3.txt #second whitelist for cell label 3
export CUSTOM_FA=~/leukemia_bwamem2/genome/BCR_ABL.fa #file containing fusion cDNAs
export NTHREADS=5 #number of threads
export TRANSCRIPTOME=gencode.v38.pc_transcripts.fa #transcriptome to download
export STARSOLO_BAM=$WD/starsolo/Aligned.sortedByCoord.out.bam #path to .bam file, generated during script

./workflow.sh
