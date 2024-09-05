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
WD=~/test_leukemia_downsampled_cell_line_experiment
COMBINED_INDEXED_GENOME=~/mapping_leukemia/data/index
STARSOLO_BAM=$WD/align_tso/downsampled_cell_line/Aligned.sortedByCoord.out.bam
r1=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/downsampled_R1_5M.fastq.gz
r2=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/downsampled_R2_5M.fastq.gz
r2_no_gz=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/downsampled_R2_5M.fastq

# for bwa mem2

TRANSCRIPTOME=gencode.v38.pc_transcripts.fa

CUSTOM_FA=~/leukemia_bwamem2/genome/BCR_ABL.fa
CUSTOM_GTF=~/leukemia_bwamem2/genome/BCR_ABL.gtf

# for STAR fusion

STAR_FUSION_GENOME=~/test_starfusion/genome/ctat_genome_lib_build_dir
PATH_STAR_FUSION=~/STAR-Fusion/ctat-genome-lib-builder
PATH_SIMG_FILE=~/star-fusion.v1.13.0.simg

# get unmapped reads and reads that mapped to BCR and to ABL and append their name to the readname

cd $WD

# for mapped reads: they will have the gx tag, unmapped won't so will have different column number

# getting reads in BCR and ABL region and adding 1000 bp on each side of the region

bcr_start=23179704
bcr_end=23318037
abl_start=130713016
abl_end=130887675

echo "getting read names and deduplicating"

# deduplication is done based on first sorting by UB (28), CB (27) and 1 (read id). If the combination of UB and CB is unique as well as the read name, then the read is appended. 

samtools view -H $STARSOLO_BAM > header.txt

samtools view $STARSOLO_BAM "chr22:22179704-24318037" | cat header.txt - | samtools view | grep -v "CB:Z:-" | grep -v "UB:Z:-" | sort -k23,23 -k22,22 -k1,1 | awk -v seen='sid' -v current="cid" -v current_ub_cb="cubcb" -v seen_ub_cb="subcb" -v current_ub="cub" -v seen_ub="sub" 'BEGIN{OFS="\t"} {
 current = $1;
 current_ub_cb = $28";"$27;
 current_ub = $23;
 if (current_ub != seen_ub && current_ub_cb != seen_ub_cb && current != seen) {
  {print "@" current ";"$27";"$28"\n"$10"\n+\n"$11};
 } seen=current; seen_ub_cb=current_ub_cb; seen_ub=current_ub;
}' >> r2.fastq

samtools view $STARSOLO_BAM "chr9:129713016-131887675" | cat header.txt - | samtools view | grep -v "CB:Z:-" | grep -v "UB:Z:-" | sort -k23,23 -k22,22 -k1,1 | awk -v seen='sid' -v current="cid" -v current_ub_cb="cubcb" -v seen_ub_cb="subcb" -v current_ub="cub" -v seen_ub="sub" 'BEGIN{OFS="\t"} {
 current = $1;
 current_ub_cb = $28";"$27;
 current_ub = $23;
 if (current_ub != seen_ub && current_ub_cb != seen_ub_cb && current != seen) {
  {print "@" current ";"$27";"$28"\n"$10"\n+\n"$11};
 } seen=current; seen_ub_cb=current_ub_cb; seen_ub=current_ub;
}' >> r2.fastq

# for unmapped reads: they don't have the gx tag, so the CB and UB will be at a different position (UB (23), CB (22) and 1 (read id)

samtools view -f 4 $STARSOLO_BAM | samtools view | grep -v "CB:Z:-" | grep -v "UB:Z:-" | sort -k23,23 -k22,22 -k1,1 | awk -v seen='sid' -v current="cid" -v current_ub_cb="cubcb" -v seen_ub_cb="subcb" -v current_ub="cub" -v seen_ub="sub" 'BEGIN{OFS="\t"} {
 current = $1; # $1 is the read id
 current_ub_cb = $23";"$22;
 current_ub = $23;
 if (current_ub != seen_ub && current_ub_cb != seen_ub_cb && current != seen) {
  {print "@" current ";"$22";"$23"\n"$10"\n+\n"$11};
 } seen=current; seen_ub_cb=current_ub_cb; seen_ub=current_ub;
}' >> r2.fastq

echo "finished getting read names and deduplicating"

rm header.txt

# problem -> read with different read name but same UB and CB can be unmapped but also in BCR / ABL region (happened with a read) --> need to check and remove the duplicate, unmapped one --> still to solve 

# first obtain the duplicated CB and UB

grep "@" r2.fastq | awk '{split ($0, a, /[;]/); print a[2]";"a[3]}' | sort | uniq -d > duplicated.txt

# obtain the entry in the .fastq file for the duplicated CB and UB. Only taking the second line as the unmapped will always be at that position

grep -Ff duplicated.txt r2.fastq | awk 'NR % 2 == 1' > unmapped_reads.txt

# remove the lines from that entry and generate new file

while read -r id; do sed -e "/$id/,+3d" r2.fastq; done < unmapped_reads.txt > deduplicated_r2.fastq

rm r2.fastq duplicated.txt unmapped_reads.txt

# deduplicateing

mv deduplicated_r2.fastq sorted_duplicate_r2.fastq

pigz sorted_duplicate_r2.fastq

echo "R2 file generated"

# run bwa_mem2

mkdir -p $WD/bwa_mem2/genome
cd $WD/bwa_mem2/genome

# genome generation

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$TRANSCRIPTOME".gz

pigz --decompress *gz

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

# index

~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p indexed combined.fa

echo "BWAMEM2 finished"

# run

mkdir -p $WD/bwa_mem2/output
cd $WD/bwa_mem2/output

~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem $WD/bwa_mem2/genome/indexed $WD/sorted_duplicate_r2.fastq.gz -t 20 -k 80 -w 12 | samtools view -bS | samtools sort -o bwa_mem2.sorted.bam

mkdir -p $WD/star_fusion/output
cd $WD/star_fusion/output

singularity exec -e -B `pwd` -B $PATH_STAR_FUSION \
       $PATH_SIMG_FILE \
       STAR-Fusion \
       --left_fq $WD/sorted_duplicate_r2.fastq.gz \
       --genome_lib_dir $STAR_FUSION_GENOME \
       --output_dir $WD/star_fusion/output &> log.txt

echo "STAR_fusion finished"

# adding the CB and UB tag to the bwa mem2 bam file and generating count table for bwa mem2

cd $WD/bwa_mem2/output

samtools view -b -F 4 bwa_mem2.sorted.bam > mapped.bam

# removed unmapped reads

samtools view -H mapped.bam > header.txt

# only some reads will have the XA tag for secondary alignments, others will be empty so need to split between the two

# extracting the entries that don't have an XA tag

samtools view mapped.bam | awk 'BEGIN {FS="\t"} $16=="" {print}'| awk 'BEGIN {OFS="\t"} {split($1, a, /[;,]/); print a[1], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, a[2], a[3]}' | cat header.txt - | samtools view -Sbh > no_xa_annotated_bwa_mem2.sorted.bam

# extracting the entries with an XA tag

samtools view mapped.bam | awk 'BEGIN {FS="\t"} $16!="" {print}'| awk 'BEGIN {OFS="\t"} {split($1, a, /[;,]/); print a[1], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, a[2], a[3]}' | cat header.txt - | samtools view -Sbh > xa_annotated_bwa_mem2.sorted.bam

rm mapped.bam

# UMI deduplication based on CB and UMI and gene name field

echo '# bwa mem2 alignments with xa tag (multimapped)'

samtools view xa_annotated_bwa_mem2.sorted.bam | cut -f 1 | wc -l

echo '# bwa mem2 alignments with no xa tag (unique)'

samtools view no_xa_annotated_bwa_mem2.sorted.bam | cut -f 1 | wc -l

# bwa mem2 does not require deduplication as we already deduplicated the .fastq files and bwa stores multimapped reads in the XA tag, which we are handling while counting. 

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

samtools view $WD/bwa_mem2/output/xa_annotated_bwa_mem2.sorted.bam | awk 'BEGIN {OFS="\t"} {n=split($16, a, /[;]/); print n,$3,$17}' | sort -k3,3 -k2,2 | uniq -f 1 -c | awk '{print $1/$2,$3,$4}' >> counts_bwa_mem2.txt

rm $WD/bwa_mem2/output/header.txt $WD/bwa_mem2/output/xa_annotated_bwa_mem2.sorted.bam $WD/bwa_mem2/output/no_xa_annotated_bwa_mem2.sorted.bam

# it may happen that the multimapping counts have the same barcode and gene as the unique counts and thus must be summed

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

# do the same for STAR fusion based on the Chimeric.out.junction
# only count the ones on the + strand
# here we are reporting multiple alignments for the same read --> depending if for the same read multiple possible fusions were detected
# we are also reporting things other than BCR:ABL if we don't filter for the chromosomes --> only take chr22 and chr9

cd $WD/star_fusion/output

# don't want the truncated reads

# just want the ones from chr9 and chr22 

bcr_start=23179704
bcr_end=23318037
abl_start=130713016
abl_end=130887675

awk '$2 >23179704  || $2 <23318037 || $4 > 130713016 || $4 < 130887675' Chimeric.out.junction > sub_Chimeric.out.junction

# UMI deduplication based on CB UB and start / end position

grep -v "-" sub_Chimeric.out.junction | grep "chr9" | grep "chr22" | awk 'BEGIN {OFS="\t"} {split($10, a, /[;,]/); print a[1],$1""$4"_"$1"_"$2"__"$4"_"$5,a[2],a[3]}' > annotated_sub_Chimeric.out.junction

echo "# chimeric alignments starfusion"
wc -l annotated_sub_Chimeric.out.junction

# pos1: read identifier, pos2: position of fusion, pos3: CB, pos4: UB

less annotated_sub_Chimeric.out.junction | sort -k3,4 | uniq -f 2 | cut -f1 > dedup_read_names.txt # uniq based on cb and ub --> all are already uniq
less annotated_sub_Chimeric.out.junction | grep -Ff dedup_read_names.txt > deduplicated_annotated_sub_Chimeric.out.junction

less deduplicated_annotated_sub_Chimeric.out.junction | cut -f2,3 | sort | uniq -c >> p_counts_star_fusion.txt

rm dedup_read_names.txt annotated_sub_Chimeric.out.junction sub_Chimeric.out.junction

# counts

mkdir -p $WD/count_tables/star_fusion
cd $WD/count_tables/star_fusion

echo -e "Counts\tfusion_position\tBarcode" > counts_star_fusion.txt
grep -v "chr9chr22" $WD/star_fusion/output/p_counts_star_fusion.txt > counts_star_fusion.txt

rm $WD/star_fusion/output/p_counts_star_fusion.txt
