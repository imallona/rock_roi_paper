#!/bin/bash
##
## bwamem2 mapping

~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem ~/leukemia_bwamem2/indexed_genome/indexed ~/simulated_leukemia_data/combined_r2.fastq.gz -t 20 -k 80 -w 12 | samtools view -bS | samtools sort -o output.sorted.bam

samtools view -b -f 4 output.sorted.bam > unmapped.bam
samtools view -b -F 4 output.sorted.bam > mapped.bam
