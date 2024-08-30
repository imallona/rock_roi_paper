#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##

## this script is to be run after the rock_roi_method workflow: https://github.com/imallona/rock_roi_method/tree/main
# only on TSO data, since the fusions are expected to be there 

# for STARsolo

COMBINED_FA_GENOME=~/mapping_leukemia_data/genome/combined.fa
COMBINED_GTF_GENOME=~/mapping_leukemia_data/genome/combined.gtf
WD=~/test_leukemia_simulated_reads
COMBINED_INDEXED_GENOME=~/mapping_leukemia/data/index
STARSOLO_BAM=$WD/starsolo/Aligned.sortedByCoord.out.bam
r1=~/simulated_leukemia_data/combined_r1.fastq.gz
r2=~/simulated_leukemia_data/combined_r2.fastq.gz
r2_no_gz=~/simulated_leukemia_data/combined_r2.fastq

# for bwa mem2

TRANSCRIPTOME=gencode.v38.pc_transcripts.fa

CUSTOM_FA=~/leukemia_bwamem2/genome/BCR_ABL.fa
CUSTOM_GTF=~/leukemia_bwamem2/genome/BCR_ABL.gtf

# for STAR fusion

STAR_FUSION_GENOME=~/test_starfusion/genome/ctat_genome_lib_build_dir
PATH_STAR_FUSION=~/STAR-Fusion/ctat-genome-lib-builder
PATH_SIMG_FILE=~/star-fusion.v1.13.0.simg

# threads

NTHREADS=10

#mkdir -p $WD/starsolo
#cd $WD/starsolo

#STAR --runThreadN $NTHREADS \
#     --genomeDir $COMBINED_INDEXED_GENOME \
#     --readFilesCommand zcat \
#     --outFileNamePrefix $WD/ \
#     --readFilesIn  $r2 $r1  \
#     --soloType CB_UMI_Complex \
#     --soloAdapterSequence AATGNNNNNNNNNCCAC \
#     --soloCBwhitelist ~/test_leukemia_simulated_reads/rock_roi_method/data/whitelist_96x3/BD_CLS1.txt ~/test_leukemia_simulated_reads/rock_roi_method/data/whitelist_96x3/BD_CLS2.txt ~/test_leukemia_simulated_reads/rock_roi_method/data/whitelist_96x3/BD_CLS3.txt \
#     --soloCBmatchWLtype 1MM \
#     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
#     --soloUMIposition 3_10_3_17 \
#     --soloUMIlen 8 \
#     --soloCellReadStats Standard \
#     --soloCellFilter None \
#     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY\
#     --outSAMtype BAM SortedByCoordinate \
#     --quantMode GeneCounts \
#     --sjdbGTFfile $COMBINED_GTF_GENOME \
#     --sjdbOverhang 179 \
#     --outSAMunmapped Within

#samtools index -@ $NTHREADS Aligned.sortedByCoord.out.bam

# get unmapped reads and reads that mapped to BCR and to ABL and append their name to the readname

#cd $WD

# for mapped reads: they will have the gx tag, unmapped won't so will have different column number

# getting reads in BCR and ABL region

#samtools view $STARSOLO_BAM "chr22:23179704-23318037" | grep -v "CB:Z:-" | grep -v "UB:Z:-" | grep -v "NH:i:0" | awk 'BEGIN{OFS="\t"}{print $1";"$27";"$28}' >> R2_names.txt
#samtools view $STARSOLO_BAM "chr9:130713016-130887675" | grep -v "CB:Z:-" | grep -v "UB:Z:-" | grep -v "NH:i:0" | awk 'BEGIN{OFS="\t"}{print $1";"$27";"$28}' >> R2_names.txt

# for unmapped reads: they don't have the gx tag, so the CB and UB will be different

#samtools view $STARSOLO_BAM | grep "NH:i:0" | grep -v "CB:Z:-" | grep -v "UB:Z:-" | awk 'BEGIN{OFS="\t"}{print $1";"$22";"$23}' >> R2_names.txt

# get the start of the readname and export the fastq record into new file

#gunzip $r2

#while read -r line; do
#    x=$(echo "$line" | awk '{split($0,a,/[;,]/); print a[1]}')
#    echo "@$line" >> r2.fastq
#    grep -A 3 "$x" "$r2_no_gz" | tail -n +2 >> r2.fastq
#done < R2_names.txt

#gzip $r2_no_gz

#rm R2_names.txt

# deduplicate based on read name

#cat r2.fastq | paste - - - - | sort -k1,1 -t " " | uniq | tr "\t" "\n" > sorted_duplicate_r2.fastq

#gzip sorted_duplicate_r2.fastq
#rm r2.fastq

#echo "R2 file generated"

# run bwa_mem2

#mkdir -p $WD/bwa_mem2/genome
#cd $WD/bwa_mem2/genome

# genome generation

#wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$TRANSCRIPTOME".gz

#pigz --decompress *gz

# problem -> featureCounts does not accept the long identifier, need to slice it 

#grep "ABL1" "$TRANSCRIPTOME"| sed 's/>//g' > ABL1_human.gtf
#grep "BCR" "$TRANSCRIPTOME"| sed 's/>//g' >  BCR_human.gtf

#touch ABL1_human.fa
#touch BCR_human.fa

#for i in $(cat ABL1_human.gtf); do
#    samtools faidx $TRANSCRIPTOME "$i" >> ABL1_human.fa
#done

#for i in $(cat BCR_human.gtf) ; do
#    samtools faidx $TRANSCRIPTOME "$i"  >> BCR_human.fa
#done

#cat "$CUSTOM_FA" ABL1_human.fa BCR_human.fa > combined.fa
#cat "$CUSTOM_GTF" ABL1_human.gtf BCR_human.gtf > combined.gtf

# index

#~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p indexed combined.fa

#echo "BWAMEM2 finished"

# run

#mkdir -p $WD/bwa_mem2/output
#cd $WD/bwa_mem2/output

#~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem $WD/bwa_mem2/genome/indexed $WD/sorted_duplicate_r2.fastq.gz -t 20 -k 80 -w 12 | samtools view -bS | samtools sort -o bwa_mem2.sorted.bam

# run STAR fusion --> the genome takes too long to generate (3 days), use premade genome

#mkdir -p $WD/star_fusion/output
#cd $WD/star_fusion/output

#singularity exec -e -B `pwd` -B $PATH_STAR_FUSION \
#       $PATH_SIMG_FILE \
#       STAR-Fusion \
#       --left_fq $WD/sorted_duplicate_r2.fastq.gz \
#       --genome_lib_dir $STAR_FUSION_GENOME \
#       --output_dir $WD/star_fusion/output > log.txt

#echo "STAR_fusion finished"

# adding the CB and UB tag to the bwa mem2 bam file and generating count table 

# only want the mapped reads --> they will have non zero XS field 

cd $WD/bwa_mem2/output

samtools view -H bwa_mem2.sorted.bam > header.txt
samtools view bwa_mem2.sorted.bam | grep -v "XS:i:0" | awk 'BEGIN {OFS="\t"} {split($1, a, /[;,]/); print a[1], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, a[2], a[3]}' | cat header.txt - | samtools view -Sbh > annotated_bwa_mem2.sorted.bam

rm header.txt

# count table for bwa_mem2 --> based on counting the occurences for each gene, cell and count. The .fastqs are already deduplicated so don't need to repeat. 
# column 1: gene id, column 2: barcode id, column 3: counts
# only count the primary alignment, even if multimapped --> one alignment for mapped read

mkdir -p $WD/count_tables/bwamem2
cd $WD/count_tables/bwamem2

echo -e "Counts\tGene_id\tBarcode" > counts_bwa_mem2.txt

samtools view $WD/bwa_mem2/output/annotated_bwa_mem2.sorted.bam | cut -f3,16 | sort | uniq -c >> counts_bwa_mem2.txt

# do the same for STAR fusion based on the Chimeric.out.junction
# only count the ones on the + strand
# here we are reporting multiple alignments for the same read --> depending if for the same read multiple possible fusions were detected
# we are also reporting things other than BCR:ABL if we don't filter for the chromosomes --> only take chr22 and chr9

mkdir -p $WD/count_tables/star_fusion
cd $WD/count_tables/star_fusion

# just want the ones from chr9 and chr22 

bcr_start=23179704
bcr_end=23318037
abl_start=130713016
abl_end=130887675

awk '$2 >23179704  || $2 <23318037 || $4 > 130713016 || $4 < 130887675' $WD/star_fusion/output/Chimeric.out.junction > $WD/star_fusion/output/sub_Chimeric.out.junction

echo -e "Counts\tfusion_position\tBarcode" > counts_star_fusion.txt
less $WD/star_fusion/output/sub_Chimeric.out.junction | grep "+" | grep "chr9" | grep "chr22" | awk 'BEGIN {OFS="\t"} {split($10, a, /[;,]/); print $1""$4"_"$1"_"$2"__"$4"_"$5, a[2]}' | sort | uniq -c >> p_counts_star_fusion.txt

grep -v "chr9chr22" p_counts_star_fusion.txt > counts_star_fusion.txt

rm p_counts_star_fusion.txt

# problem --> STARSolo will most likely have counted the unique reads in the .bam file as WT BCR or ABL, so they are being counted twice. How to deal with this?
# not an issue --> for bwa mem2 we are counting the wt separately, so remove the wt BCR and wt ABL counts from the matrix and use the ones detected from bwa mem2
# not an issue --> for STARfusion they won't be detected 
