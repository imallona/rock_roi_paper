#!/bin/bash
##

WD=/home/gmoro/simulated_kiel_mapping/
GENOME=/home/gmoro/indices
EXPERIMENT=bulk_small_genome_standard
FEATURECOUNTS=/home/imallona/soft/subread-2.0.4-source/bin/featureCounts

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
     --sjdbGTFfile ~/kiel_leukemia_data/mapping_patient_data/captured_bcr_abl.gtf \
     --outFilterMismatchNmax 10 \
     --seedSearchStartLmax 50 \
     --outFilterMatchNmin 0

# to get the information on the CIGAR

samtools view $EXPERIMENT"Aligned.sortedByCoord.out.bam" | cut -f 1,3,6 > $EXPERIMENT"_CIGAR.txt"

# counting with custom gtf

$FEATURECOUNTS -a ~/kiel_leukemia_data/mapping_patient_data/captured_bcr_abl.gtf \
                -o "$item" \
                $EXPERIMENT"Aligned.sortedByCoord.out.bam" \
                -F GTF \
                -t exon \
                -g gene_id \
                -f \
                -O \
                -M  \
                -T "$NTHREADS" \
                --byReadGroup
