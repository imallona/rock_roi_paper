#!/bin/bash
##
## Generates a genome for our leukemia experiment, including sampletags
##
## 21st Feb 2024
## Giulia Moro

WD=~/kiel_leukemia_data/mapping_patient_data
mkdir -p $WD/genomes
cd $WD/genomes

tab="$(printf '\t')"

human_fa=GRCh38.p13.genome.fa
human_gtf=gencode.v38.basic.annotation.gtf

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_fa".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_gtf".gz

pigz --decompress *gz

rm $human_fa $human_gtf
