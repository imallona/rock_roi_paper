#!/bin/bash
##

singularity exec -e -B `pwd` -B ~/STAR-Fusion/ctat-genome-lib-builder \
       ~/star-fusion.v1.13.0.simg \
       STAR-Fusion \
       --left_fq ~/simulated_leukemia_data/combined_r2.fastq.gz \
       --genome_lib_dir ~/test_starfusion/genome/ctat_genome_lib_build_dir \
       --output_dir ~/test_starfusion/simulated_data > log.txt
