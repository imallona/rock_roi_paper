#!/bin/bash
##

WD=/home/gmoro/simulated_kiel_mapping
GENOME=/home/gmoro/indices

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

NTHREADS=10

for value in 0 1 1000000
do
	mkdir $WD/test_alignMatesGapMax/$EXPERIMENT
	cd $WD/test_alignMatesGapMax/$EXPERIMENT

	EXPERIMENT=outFilterMismatchNmax_"$value"_bulk_small_genome_standard
	echo $EXPERIMENT

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
     	--seedSearchStartLmax 50 \
        --alignMatesGapMax $value

	samtools view $EXPERIMENT"Aligned.sortedByCoord.out.bam" | cut -f 1,3,6 > $EXPERIMENT"_CIGAR.txt"
done
