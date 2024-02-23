#!/bin/bash
##
## Optimize 62 nt-long cDNA alignment to 45 nt-long sampletags
##
## 07th Nov 2023
## Izaskun Mallona



WD=~/ebrunner_spectral/leukemia_sampletags
mkdir -p $WD/genomes
cd $WD/genomes


tab="$(printf '\t')"
prepend=GTTGTCAAGATGCTACCGTTCAGAG

cat <<EOF > sampletags.fa
>human_sampletag_1
GTTGTCAAGATGCTACCGTTCAGAGATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG
>human_sampletag_2
GTTGTCAAGATGCTACCGTTCAGAGTGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG
>human_sampletag_3
GTTGTCAAGATGCTACCGTTCAGAGCGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT
>human_sampletag_4
GTTGTCAAGATGCTACCGTTCAGAGATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT
>human_sampletag_5
GTTGTCAAGATGCTACCGTTCAGAGCTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG
>human_sampletag_6
GTTGTCAAGATGCTACCGTTCAGAGTTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC
>human_sampletag_7
GTTGTCAAGATGCTACCGTTCAGAGTGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT
>human_sampletag_8
GTTGTCAAGATGCTACCGTTCAGAGCCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT
>human_sampletag_9
GTTGTCAAGATGCTACCGTTCAGAGGTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG
>human_sampletag_10
GTTGTCAAGATGCTACCGTTCAGAGGCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC
>human_sampletag_11
GTTGTCAAGATGCTACCGTTCAGAGCGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC
>human_sampletag_12
GTTGTCAAGATGCTACCGTTCAGAGGCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG
EOF

cat <<EOF > sampletags.gtf
human_sampletag_1${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_1"; transcript_id "human_sampletag_1";
human_sampletag_2${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_2"; transcript_id "human_sampletag_2";
human_sampletag_3${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_3"; transcript_id "human_sampletag_3";
human_sampletag_4${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_4"; transcript_id "human_sampletag_4";
human_sampletag_5${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_5"; transcript_id "human_sampletag_5";
human_sampletag_6${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_6"; transcript_id "human_sampletag_6";
human_sampletag_7${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_7"; transcript_id "human_sampletag_7";
human_sampletag_8${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_8"; transcript_id "human_sampletag_8";
human_sampletag_9${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_9"; transcript_id "human_sampletag_9";
human_sampletag_10${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_10"; transcript_id "human_sampletag_10";
human_sampletag_11${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_11"; transcript_id "human_sampletag_11";
human_sampletag_12${tab}sampletag${tab}exon${tab}1${tab}70${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_12"; transcript_id "human_sampletag_12";
EOF

## no GTF, no splicing later on
STAR --runThreadN 5 \
     --runMode genomeGenerate \
     --genomeSAindexNbases 2 \
     --genomeDir sampletags_genome/ \
     --genomeFastaFiles sampletags.fa
     # --sjdbOverhang 179 \
     # --sjdbGTFfile sampletags.gtf \

cd $WD

ln -s ~/src/rock_roi_method/data/whitelist_96x3/BD_CLS1.txt .
ln -s ~/src/rock_roi_method/data/whitelist_96x3/BD_CLS2.txt .
ln -s ~/src/rock_roi_method/data/whitelist_96x3/BD_CLS3.txt .

cdna=/home/imallona/ebrunner_spectral/data/fastqs_short_leukemia/r2_patient_samples.fastq.gz
cbumi=/home/imallona/ebrunner_spectral/data/fastqs_short_leukemia/r1_patient_samples.fastq.gz


## downsampled
# cd ~/ebrunner_spectral/data/fastqs_short_leukemia/
# zcat $cdna | head -500000 | gzip -c > r2_pat_down.fastq.gz
# zcat $cbumi | head -500000 | gzip -c > r1_pat_down.fastq.gz

# cdna=/home/imallona/ebrunner_spectral/data/fastqs_short_leukemia/r2_pat_down.fastq.gz
# cbumi=/home/imallona/ebrunner_spectral/data/fastqs_short_leukemia/r1_pat_down.fastq.gz


## notice the alignIntronMax 1 and the seedSearchStartLmax
##  and the outFilter* (given the read length v reference's)
ulimit -Sn 1000000

nice -n 19 STAR --runThreadN 80 \
     --genomeDir genomes/sampletags_genome/ \
     --readFilesCommand zcat \
     --outFileNamePrefix patients_wta_vs_sampletags/ \
     --readFilesIn $cdna $cbumi  \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GTGANNNNNNNNNGACA \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist BD_CLS1.txt BD_CLS2.txt BD_CLS3.txt  \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter None \
     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY \
     --soloCellReadStats Standard \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --soloUMIlen 8 \
     --sjdbGTFfile genomes/sampletags.gtf \
     --outFilterScoreMinOverLread 0.5 \
     --outFilterMatchNminOverLread 0.5 \
     --outFilterMismatchNmax 25 \
     --seedSearchStartLmax 30 \
     --alignIntronMax 1 \
     --limitBAMsortRAM 9093676412


# less patients_wta_vs_sampletags/Solo.out/Gene/raw/matrix.mtx
