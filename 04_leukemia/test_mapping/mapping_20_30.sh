#!/bin/bash
##

WD=/home/gmoro/kiel_leukemia_data/
GENOME=/home/gmoro/indices
CB_dir=/home/gmoro/whitelist_96x3

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

NTHREADS=10
nice -n 19 STAR --runThreadN $NTHREADS \
     --genomeDir $GENOME/GRCh38_gencode_38 \
     --readFilesCommand zcat \
     --outFileNamePrefix test_mapping \
     --readFilesIn $WD/jf_bastian_rock_ro915_hjyltdrx3/simulated_r2.fastq.gz $WD/jf_bastian_rock_ro915_hjyltdrx3/simulated_r1.fastq.gz \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence AATGNNNNNNNNNCCAC \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist $CB_dir/BD_CLS1.txt $CB_dir/BD_CLS2.txt $CB_dir/BD_CLS3.txt  \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter None \
     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY \
     --soloCellReadStats Standard \
     --sjdbOverhang 61 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode - \
     --soloUMIlen 8 \
     --sjdbGTFfile $GENOME/gencode.v38.basic.annotation.gtf \
     --outFilterMismatchNmax 20 \
     --seedSearchStartLmax 30
