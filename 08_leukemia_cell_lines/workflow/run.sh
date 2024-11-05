#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH #path to STAR for STAR fusion
export PATH=/home/gmoro/bwa/:$PATH #path to bwa
export WD=~/cell_lines_leukemia/rock_roi_repo_for_cell_lines/workflow_aln #working directory
export COMBINED_INDEXED_GENOME=~/mapping_leukemia/data/index # STAR indexed genome for STAR fusion
export STARSOLO_BAM=/home/gmoro/cell_lines_leukemia/rock_roi_repo_for_cell_lines/output/align_tso/leukemia_cell_lines/Aligned.sortedByCoord.out.bam # TSO .bam file generated with rock_roi_method 
export CUSTOM_FA=~/leukemia_bwamem2/genome/BCR_ABL.fa #file containing fusion cDNAs
export NTHREADS=30 #number of threads
export TRANSCRIPTOME=gencode.v38.pc_transcripts.fa #transcriptome to download
export REGIONS=~/rock_roi_paper/07_leukemia_patient/simulation_patient_cell_lines/regions.bed

# installing conda and packages 

#mkdir -p ~/miniconda3
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
#bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

source ~/miniconda3/bin/activate

conda create --name fusion_detection
conda activate fusion_detection

conda config --add channels bioconda
conda config --add channels conda-forge

conda install -c conda-forge -c bioconda star=2.7.10b
conda install -c conda-forge -c bioconda samtools=1.21 # v1.21
conda install -c conda-forge -c bioconda bwa=0.7.18 #0.7.18

./workflow.sh
