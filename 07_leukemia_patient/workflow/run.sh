#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH #path to STAR for STAR fusion
export PATH=/home/gmoro/bwa/:$PATH #path to bwa
export WD=~/data_mapping_leukemia_patient/workflow_aln #working directory
export COMBINED_INDEXED_GENOME=~/mapping_leukemia/data/index # STAR indexed genome for STAR fusion
export STARSOLO_BAM=/home/gmoro/kiel_leukemia_data/mapping_patient_data/align_tso/leukemia_patient/Aligned.sortedByCoord.out.bam # TSO .bam file generated with rock_roi>
export CUSTOM_FA=~/leukemia_bwamem2/genome/BCR_ABL.fa #file containing fusion cDNAs
export NTHREADS=30 #number of threads
export TRANSCRIPTOME=gencode.v38.pc_transcripts.fa #transcriptome to download
export REGIONS=~/rock_roi_paper/07_leukemia_patient/simulation_patient_cell_lines/regions.bed

./workflow.sh
