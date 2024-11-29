#!/bin/bash
##
## Generating .bam file using STARsolo
##
## Started 29th of November 2024
##
## Giulia Moro

#conda install -c conda-forge -c bioconda star=2.7.10b
#conda install -c conda-forge -c bioconda samtools=1.21 # v1.21

gzip ~/rock_roi_paper/09_fusion_simulations/data/fusion_simulations_cbumi.fq
gzip ~/rock_roi_paper/09_fusion_simulations/data/fusion_simulations_cdna.fq

STAR --runThreadN 5 \
     --genomeDir ~/rock_roi_paper/07_leukemia_patient/simulation_patient_cell_lines/genome \
     --readFilesCommand zcat \
     --outFileNamePrefix ./starsolo/ \
     --readFilesIn ~/rock_roi_paper/09_fusion_simulations/data/fusion_simulations_cdna.fq.gz ~/rock_roi_paper/09_fusion_simulations/data/fusion_simulations_cbumi.fq.gz  \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence AATGNNNNNNNNNCCAC \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloUMIlen 8 \
     --soloCellReadStats Standard \
     --soloCBwhitelist ~/new_whitelist/line_extended_cell_label_1.txt ~/new_whitelist/line_extended_cell_label_2.txt ~/new_whitelist/line_extended_cell_label_3.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter None \
     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY\
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile ~/rock_roi_paper/07_leukemia_patient/simulation_patient_cell_lines/genome/gencode.v38.basic.annotation.gtf \
     --sjdbOverhang 179 \
     --limitBAMsortRAM 20000 * 1024 * 1024 \
     --outSAMunmapped Within

samtools index ./starsolo/Aligned.sortedByCoord.out.bam

gunzip ~/rock_roi_paper/09_fusion_simulations/data/fusion_simulations_cbumi.fq.gz
gunzip ~/rock_roi_paper/09_fusion_simulations/data/fusion_simulations_cdna.fq.gz
