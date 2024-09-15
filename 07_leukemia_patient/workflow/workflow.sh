#!/bin/bash
##
## Giulia Moro
##
## Fusion analysis for the leukemia dataset
##

## this script is to be run after the rock_roi_method workflow: https://github.com/imallona/rock_roi_method/tree/main
# only on TSO data, since the fusions are expected to be there

# STAR path 

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

# for STARsolo

COMBINED_FA_GENOME=~/mapping_leukemia_data/genome/combined.fa
COMBINED_GTF_GENOME=~/mapping_leukemia_data/genome/combined.gtf
WD=~/test_leukemia_simulated_reads
COMBINED_INDEXED_GENOME=~/mapping_leukemia/data/index
STARSOLO_BAM=$WD/starsolo/Aligned.sortedByCoord.out.bam
r1=/home/gmoro/simulated_leukemia_data/combined_r1.fastq.gz
r2=/home/gmoro/simulated_leukemia_data/combined_r2.fastq.gz
r2_no_gz=/home/gmoro/simulated_leukemia_data/combined_r1.fastq

# for bwa mem2

TRANSCRIPTOME=gencode.v38.pc_transcripts.fa

CUSTOM_FA=~/leukemia_bwamem2/genome/BCR_ABL.fa
CUSTOM_GTF=~/leukemia_bwamem2/genome/BCR_ABL.gtf

# for STAR fusion

STAR_FUSION_GENOME=~/test_starfusion/genome/ctat_genome_lib_build_dir
PATH_STAR_FUSION=~/STAR-Fusion/ctat-genome-lib-builder
PATH_SIMG_FILE=~/star-fusion.v1.13.0.simg

# threads

NTHREADS=20

# running starsolo

#STAR --runThreadN $NTHREADS \
#     --genomeDir ~/mapping_leukemia/data/index \
#     --readFilesCommand zcat \
#     --outFileNamePrefix $WD/starsolo/ \
#     --readFilesIn $r2 $r1  \
#     --soloType CB_UMI_Complex \
#     --soloAdapterSequence AATGNNNNNNNNNCCAC \
#     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
#     --soloUMIposition 3_10_3_17 \
#     --soloUMIlen 8 \
#     --soloCellReadStats Standard \
#     --soloCBwhitelist ~/whitelist_96x3/BD_CLS1.txt ~/whitelist_96x3/BD_CLS2.txt ~/whitelist_96x3/BD_CLS3.txt \
#     --soloCBmatchWLtype 1MM \
#     --soloCellFilter None \
#     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY\
#     --outSAMtype BAM SortedByCoordinate \
#     --quantMode GeneCounts \
#     --sjdbGTFfile $COMBINED_GTF_GENOME \
#     --sjdbOverhang 179 \
#     --limitBAMsortRAM 20000 * 1024 * 1024 \
#     --outSAMunmapped Within

samtools index $STARSOLO_BAM

# get unmapped reads and reads that mapped to BCR and to ABL and append their name to the readname

## Step 1: get unmapped reads, reads in the BCR ABL slice, reads with no gx tag from .bam file and output in new .fastq file. Append the CB and UB tag to their read id. Deduplicate CB and UB.  

cd $WD

# getting reads in BCR and ABL region and adding 1M bp to each side

bcr_start=23179704
bcr_end=23318037
abl_start=130713016
abl_end=130887675

echo "getting read names and deduplicating"

# deduplication is done based on first sorting by UB (28), CB (27) and 1 (read id). If the combination of UB and CB is unique as well as the read name, then the read is appended. 
# if first if loop check if read is already present in .fastq file or not, otherwise not added. 

samtools view -H $STARSOLO_BAM -@ $NTHREADS > header.txt

samtools view "$STARSOLO_BAM" "chr22:22179704-24318037" -@ "$NTHREADS" \
| cat header.txt - \
| samtools view -@ "$NTHREADS" \
| grep -v "CB:Z:-" \
| grep -v "UB:Z:-" \
| sort -k28,28 -k27,27 -k1,1 \
| parallel --pipe --block 10M -j "$NTHREADS" \
'awk -v seen="" -v current="" -v current_ub_cb="" -v seen_ub_cb="" -v current_ub="" -v seen_ub="" \
  "BEGIN { OFS=\"\t\" } \
  { \
    current = \$1; \
    current_ub_cb = \$28 \";\" \$27; \
    current_ub = \$28; \
    cmd = \"grep -q \" current \" r2.fastq\"; \
    exit_code = system(cmd); \
    if (exit_code == 0) { next; } \
    if (current_ub != seen_ub && current_ub_cb != seen_ub_cb && current != seen) { \
      print \"@\" current \";\" \$27 \";\" \$28 \"\\n\" \$10 \"\\n+\\n\" \$11; \
    } \
    seen = current; \
    seen_ub_cb = current_ub_cb; \
    seen_ub = current_ub; \
  }"' >> r2.fastq

samtools view "$STARSOLO_BAM" "chr9:129713016-131887675" -@ "$NTHREADS" \
| cat header.txt - \
| samtools view -@ "$NTHREADS" \
| grep -v "CB:Z:-" \
| grep -v "UB:Z:-" \
| sort -k28,28 -k27,27 -k1,1 \
| parallel --pipe --block 10M -j "$NTHREADS" \
'awk -v seen="" -v current="" -v current_ub_cb="" -v seen_ub_cb="" -v current_ub="" -v seen_ub="" \
  "BEGIN { OFS=\"\t\" } \
  { \
    current = \$1; \
    current_ub_cb = \$28 \";\" \$27; \
    current_ub = \$28; \
    cmd = \"grep -q \" current \" r2.fastq\"; \
    exit_code = system(cmd); \
    if (exit_code == 0) { next; } \
    if (current_ub != seen_ub && current_ub_cb != seen_ub_cb && current != seen) { \
      print \"@\" current \";\" \$27 \";\" \$28 \"\\n\" \$10 \"\\n+\\n\" \$11; \
    } \
    seen = current; \
    seen_ub_cb = current_ub_cb; \
    seen_ub = current_ub; \
  }"' >> r2.fastq

# for unmapped reads: they don't have the gx tag, so the CB and UB will be at a different position (UB (23), CB (22) and 1 (read id)

samtools view -b -f 4 "$STARSOLO_BAM" -@ "$NTHREADS" \
| samtools view -@ "$NTHREADS" \
| grep -v "CB:Z:-" \
| grep -v "UB:Z:-" \
| sort -k23,23 -k22,22 -k1,1 \
| parallel --pipe --block 10M -j "$NTHREADS" \
'awk -v seen="" -v current="" -v current_ub_cb="" -v seen_ub_cb="" -v current_ub="" -v seen_ub="" \
  "BEGIN { OFS=\"\t\" } \
  { \
   current = \$1; \
   current_ub_cb = \$23 \";\" \$22; \
   current_ub = \$23; \
   cmd = \"grep -q \" current \" r2.fastq\"; \
   exit_code = system(cmd); \
   if (exit_code == 0) { next; } \
   if (current_ub != seen_ub && current_ub_cb != seen_ub_cb && current != seen) { \
     print \"@\" current \";\" \$22 \";\" \$23 \"\\n\" \$10 \"\\n+\\n\" \$11; \
   } \
   seen = current; \
   seen_ub_cb = current_ub_cb; \
   seen_ub = current_ub; \
  }"' >> r2.fastq

# also need to get the reads which don't have a gx tag as the minor fusions seem to map to a genomic region where no gene is present. Only do this based on the mapped reads

samtools view -b -F 4 "$STARSOLO_BAM" -@ "$NTHREADS" \
| samtools view -@ $NTHREADS \
| grep "gx:Z:-" \
| grep -v "CB:Z:-" \
| grep -v "UB:Z:-" \
| sort -k28,28 -k27,27 -k1,1 \
| parallel --pipe --block 10M -j "$NTHREADS" \
'awk -v seen="" -v current="" -v current_ub_cb="" -v seen_ub_cb="" -v current_ub="" -v seen_ub="" \
  "BEGIN { OFS=\"\t\" } \
  { \
    current = \$1; \
    current_ub_cb = \$28 \";\" \$27; \
    current_ub = \$28; \
    cmd = \"grep -q \" current \" r2.fastq\"; \
    exit_code = system(cmd); \
    if (exit_code == 0) { next; } \
    if (current_ub != seen_ub && current_ub_cb != seen_ub_cb && current != seen) { \
      print \"@\" current \";\" \$27 \";\" \$28 \"\\n\" \$10 \"\\n+\\n\" \$11; \
    } \
    seen = current; \
    seen_ub_cb = current_ub_cb; \
    seen_ub = current_ub; \
  }"' >> r2.fastq

echo "finished getting read names and deduplicating"

rm header.txt

## Step 2: additonal CB / UB deduplication

# Additional CB and UB deduplication if reads appended in different steps have duplication

grep "@" r2.fastq | awk '{split ($0, a, /[;]/); print a[2]";"a[3]}' | sort | uniq -d > duplicated.txt

echo 'remaining reads to deduplicate'

wc -l duplicated.txt

# obtain the entry in the .fastq file for the duplicated CB and UB. Only taking the second line as the unmapped will always be at that position

grep -Ff duplicated.txt r2.fastq | awk 'NR % 2 == 1' > unmapped_reads.txt

# remove the lines from that entry and generate new file --> only works if the wc -l of the unmapped_reads.txt is larger than 0, so need to add that condition

awk 'NF {exit 1}' unmapped_reads.txt && mv r2.fastq deduplicated_r2.fastq || mv r2.fastq deduplicated_r2.fastq; while read -r id; do sed -i "/$id/,+3d" deduplicated_r2.fastq; done < unmapped_reads.txt

rm duplicated.txt unmapped_reads.txt

mv deduplicated_r2.fastq sorted_duplicate_r2.fastq

pigz sorted_duplicate_r2.fastq

echo "R2 file generated"

## Step 3: running bwa_mem2

mkdir -p $WD/bwa_mem2/genome
cd $WD/bwa_mem2/genome

# transcript reference generation for bwa_mem2

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$TRANSCRIPTOME".gz

pigz --decompress *gz

# ENSG00000097007: ABL1
# ENSG00000186716: BCR

grep "ENSG00000097007" "$TRANSCRIPTOME"| sed 's/>//g' > ABL1_human.gtf
grep "ENSG00000186716" "$TRANSCRIPTOME"| sed 's/>//g' >  BCR_human.gtf

for i in $(cat ABL1_human.gtf); do
    samtools faidx $TRANSCRIPTOME "$i" >> ABL1_human.fa
done

 for i in $(cat BCR_human.gtf) ; do
    samtools faidx $TRANSCRIPTOME "$i"  >> BCR_human.fa
done

cat "$CUSTOM_FA" ABL1_human.fa BCR_human.fa > combined.fa
cat "$CUSTOM_GTF" ABL1_human.gtf BCR_human.gtf > combined.gtf

# index reference

~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p indexed combined.fa

mkdir -p $WD/bwa_mem2/output
cd $WD/bwa_mem2/output

# -k 80: seed length
# -h 1: if there are <INT hits with score >80% of the max score, output all in XA. Only outputting 1 extra to see if there are multialigners
# -L: penalty for soft clipping
# -A: matching score
# -B: mismatch penalty

~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem $WD/bwa_mem2/genome/indexed $WD/sorted_duplicate_r2.fastq.gz -t 30 -k 80 -w 12 -h 1 -L 1000 -A 100 -B 100 | samtools view -bS | samtools sort -o output.sorted.bam

## Step 4: run STAR fusion

# genome was previously generated based on the .fa and .gtf files also used for STARsolo

mkdir -p $WD/star_fusion/output
cd $WD/star_fusion/output

STAR --genomeDir "$COMBINED_INDEXED_GENOME" \
--outReadsUnmapped None \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 8 \
--chimOutJunctionFormat 1 \
--alignSJDBoverhangMin 10 \
--alignMatesGapMax 100000 \
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--runThreadN 4 \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--alignInsertionFlush Right \
--alignSplicedMateMapLminOverLmate 0 \
--alignSplicedMateMapLmin 30 \
--outSAMtype BAM Unsorted \
--readFilesIn $WD/sorted_duplicate_r2.fastq.gz \
--outSAMattrRGline ID:GRPundef \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimMultimapNmax 20 \
--chimOutType Junctions WithinBAM \
--chimNonchimScoreDropMin 10 \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1 \
--genomeLoad NoSharedMemory \
--twopassMode None \
--readFilesCommand zcat \
--quantMode GeneCounts \
--outFilterMismatchNmax 1

echo "STAR_fusion finished"

## Step 5: count bwa mem2 data

# adding the CB and UB tag to the bwa mem2 bam file and generating count table for bwa mem2 after extracting mapped reads

cd $WD/bwa_mem2/output

samtools view -b -F 4 output.sorted.bam > mapped.bam

samtools view -H mapped.bam > header.txt

# only some reads will have the XA tag for secondary alignments (multialigners), others will be empty so need to split between the two (unique)
# also want to remove the reads that have an SA tag, which means that they are chimeric alignments which are split into two. In the simulated data this was the case for the ABL_BCR. 

# extracting the entries that don't have an XA tag

samtools view mapped.bam | grep -v "SA:" | awk 'BEGIN {FS="\t"} $16=="" {print}'| awk 'BEGIN {OFS="\t"} {split($1, a, /[;,]/); print a[1], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, a[2], a[3]}' | cat header.txt - | samtools view -Sbh > no_xa_annotated_bwa_mem2.sorted.bam

# extracting the entries with an XA tag

samtools view mapped.bam | grep -v "SA:" |  awk 'BEGIN {FS="\t"} $16!="" {print}'| awk 'BEGIN {OFS="\t"} {split($1, a, /[;,]/); print a[1], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, a[2], a[3]}' | cat header.txt - | samtools view -Sbh > xa_annotated_bwa_mem2.sorted.bam

rm mapped.bam

echo '# bwa mem2 alignments with xa tag (multimapped)'

samtools view xa_annotated_bwa_mem2.sorted.bam | cut -f 1 | wc -l

echo '# bwa mem2 alignments with no xa tag (unique)'

samtools view no_xa_annotated_bwa_mem2.sorted.bam | cut -f 1 | wc -l

# bwa mem2 does not require deduplication as we already deduplicated the .fastq files and bwa stores multialigned reads in the XA tag, which we are handling while counting. 

# count table for bwa_mem2 --> based on counting the occurences for each gene, cell and count. 
# column 1: gene id, column 2: barcode id, column 3: counts
# primary alignments are counted as one 

mkdir -p $WD/count_tables/bwamem2
cd $WD/count_tables/bwamem2

echo -e "Counts\tGene_id\tBarcode" > counts_bwa_mem2.txt

# sort based on k2 which is the CB field and then k1 which is the gene name and then count 

samtools view $WD/bwa_mem2/output/no_xa_annotated_bwa_mem2.sorted.bam | cut -f3,16 | sort -k2,2 -k1,1 | uniq -c >> counts_bwa_mem2.txt

# secondary alignments are handled as fractonary counts

# the ; at the end of the line will also be counted during the split command. Not subtracting 1 as in the XA tag the primary alignment is not counted and want the total number of reads (primary + multimapped), so would have to add 1 later
# uniq -f -1 -c: don't take into account the number of reported alignments for uniq counting
# sort based on CB field which here is k3 and then k2 for the gene id

# IMPORTANT: the wt BCR will always multimap. There are two annotated transcripts and they have the same exact sequence at the fusion junction. We are summing the multimappers and recovering both, so the count will anyway sum up to 1. 

samtools view $WD/bwa_mem2/output/xa_annotated_bwa_mem2.sorted.bam | awk 'BEGIN {OFS="\t"} {n=split($16, a, /[;]/); print n,$3,$17}' | sort -k3,3 -k2,2 | uniq -f 1 -c | awk '{print $1/$2,$3,$4}' >> counts_bwa_mem2.txt

rm $WD/bwa_mem2/output/header.txt

# need to sum the multimappers

less counts_bwa_mem2.txt | awk '{print $1,$2";"$3}' | sort -k2,2 | awk -v seen_gene_cb="sgc" -v current_gene_cb="cgc" -v counts_current="cc" -v counts_seen="cs" 'BEGIN{OFS="\t"} {
 counts_current = $1;
 current_gene_cb = $2";"$3;
  if (current_gene_cb==seen_gene_cb) {
  {print counts_current+counts_seen"\t"$2"\t"$3};
}
 else {
 {print counts_current "\t"$2"\t"$3};
 } seen_gene_cb=current_gene_cb; counts_seen=counts_current;
}' > combined_counts_bwa_mem2.txt

rm counts_bwa_mem2.txt

## Step 6: count STAR fusion data

# do the same for STAR fusion based on the Chimeric.out.junction
# only count the ones on the + strand
# here we are reporting multiple alignments for the same read --> depending if for the same read multiple possible fusions were detected
# we are also reporting things other than BCR:ABL if we don't filter for the chromosomes --> only take chr22 and chr9 in the correct order

cd $WD/star_fusion/output

bcr_start=23179704
bcr_end=23318037
abl_start=130713016
abl_end=130887675

awk '$2 >23179704  || $2 <23318037 || $4 > 130713016 || $4 < 130887675' Chimeric.out.junction > sub_Chimeric.out.junction

# $1: first chr, $4: second chr, $2: position in first chromosome, $5: position in second chromosome

grep -v "-" sub_Chimeric.out.junction | grep "chr9" | grep "chr22" | awk 'BEGIN {OFS="\t"} {split($10, a, /[;,]/); print a[1],$1"_"$4"_"$2"__"$5,a[2],a[3]}' > annotated_sub_Chimeric.out.junction

echo "# chimeric alignments starfusion"
wc -l annotated_sub_Chimeric.out.junction

less annotated_sub_Chimeric.out.junction | cut -f2,3 | sort | uniq -c >> p_counts_star_fusion.txt

rm sub_Chimeric.out.junction

# counts

mkdir -p $WD/count_tables/star_fusion
cd $WD/count_tables/star_fusion

echo -e "Counts\tfusion_position\tBarcode" > counts_star_fusion_complete.txt

# remove chr9 and chr22 as this would be the inversion of the fusion

grep -v "chr9chr22" $WD/star_fusion/output/p_counts_star_fusion.txt > counts_star_fusion.txt

rm $WD/star_fusion/output/p_counts_star_fusion.txt counts_star_fusion_complete.txt
