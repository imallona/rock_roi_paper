#!/bin/bash
##

WD=/home/gmoro/simulated_kiel_mapping/
GENOME=/home/gmoro/indices
EXPERIMENT=bulk_small_genome_standard

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

NTHREADS=10
nice -n 19 STAR --runThreadN $NTHREADS \
     --genomeDir $GENOME/bcr_abl_downsampled \
     --readFilesCommand zcat \
     --outFileNamePrefix $EXPERIMENT \
     --readFilesIn $WD/pos_neg_r2.fastq.gz \
     --outSAMattributes All \
     --sjdbOverhang 61 \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbGTFfile $GENOME/bcr_abl_downsampled.gtf \
     --outFilterMismatchNmax 999 \
     --seedSearchStartLmax 50

# to get the information on the CIGAR

samtools view $EXPERIMENT"Aligned.sortedByCoord.out.bam" | cut -f 1,6 > $EXPERIMENT"_CIGAR.txt"
