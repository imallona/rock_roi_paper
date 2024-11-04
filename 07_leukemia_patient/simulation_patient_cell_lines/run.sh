#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##

## The code needs to be run from the cloned repository in the 07_leukemia_patient_workflow/simulation_patient_cell_lines

export WD=.

export COMBINED_GTF_GENOME=./genome/gencode.v38.basic.annotation.gtf #genome .gtf
export r1=./combined_r1.fastq #simulated reads
export r2=./combined_r2.fastq #simulated reads
export CL1=./BD_CLS1.txt #first whitelist for cell label 1
export CL2=./BD_CLS2.txt #second whitelist for cell label 2
export CL3=./BD_CLS3.txt #second whitelist for cell label 3
export CUSTOM_FA=./BCR_ABL.fa #file containing fusion cDNAs

export NTHREADS=5 #number of threads. If they are higher than 5 STARSolo does not run
export NTHREADS_genome=30 #number of threads for genome generation

export STARSOLO_BAM=./starsolo/Aligned.sortedByCoord.out.bam #path to .bam file, generated during script

export TRANSCRIPTOME=gencode.v38.pc_transcripts.fa #transcriptome to download
export human_fa=GRCh38.p13.genome.fa
export human_gtf=gencode.v38.basic.annotation.gtf

# software installation

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

conda deactivate fusion_detection
