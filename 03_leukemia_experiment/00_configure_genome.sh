#!/bin/bash
##
## Generates a genome for our leukemia experiment, including sampletags
##
## Started 06 Nov 2023
## Izaskun Mallona


WD=~/ebrunner_spectral/leukemia
mkdir -p $WD/genomes
cd $WD/genomes


tab="$(printf '\t')"


cat <<EOF > sampletags.fa
>human_sampletag_1
ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG
>human_sampletag_2
TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG
>human_sampletag_3
CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT
>human_sampletag_4
ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT
>human_sampletag_5
CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG
>human_sampletag_6
TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC
>human_sampletag_7
TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT
>human_sampletag_8
CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT
>human_sampletag_9
GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG
>human_sampletag_10
GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC
>human_sampletag_11
CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC
>human_sampletag_12
GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG
EOF

cat <<EOF > sampletags.gtf
human_sampletag_1${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_1"; transcript_id "human_sampletag_1";
human_sampletag_2${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_2"; transcript_id "human_sampletag_2";
human_sampletag_3${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_3"; transcript_id "human_sampletag_3";
human_sampletag_4${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_4"; transcript_id "human_sampletag_4";
human_sampletag_5${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_5"; transcript_id "human_sampletag_5";
human_sampletag_6${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_6"; transcript_id "human_sampletag_6";
human_sampletag_7${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_7"; transcript_id "human_sampletag_7";
human_sampletag_8${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_8"; transcript_id "human_sampletag_8";
human_sampletag_9${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_9"; transcript_id "human_sampletag_9";
human_sampletag_10${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_10"; transcript_id "human_sampletag_10";
human_sampletag_11${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_11"; transcript_id "human_sampletag_11";
human_sampletag_12${tab}sampletag${tab}exon${tab}1${tab}45${tab}.${tab}+${tab}.${tab}gene_id "human_sampletag_12"; transcript_id "human_sampletag_12";
EOF


cat <<EOF > captured.gtf
bcr_abl1_e1_a2${tab}captured${tab}exon${tab}1${tab}333${tab}.${tab}+${tab}.${tab}gene_id "e1a2_before_roi"; transcript_id "e1a2_before_roi";
bcr_abl1_e1_a2${tab}captured${tab}exon${tab}334${tab}395${tab}.${tab}+${tab}.${tab}gene_id "e1a2_roi"; transcript_id "e1a2_roi";
bcr_abl1_e1_a2${tab}captured${tab}exon${tab}396${tab}509${tab}.${tab}+${tab}.${tab}gene_id "e1a2_between_roi_rock"; transcript_id "e1a2_between_roi_rock";
bcr_abl1_e1_a2${tab}captured${tab}exon${tab}510${tab}534${tab}.${tab}+${tab}.${tab}gene_id "e1a2_rock"; transcript_id "e1a2_rock";
bcr_abl1_e1_a2${tab}captured${tab}exon${tab}535${tab}799${tab}.${tab}+${tab}.${tab}gene_id "e1a2_after_rock"; transcript_id "e1a2_after_rock";
bcr_abl1_e13_a2${tab}captured${tab}exon${tab}1${tab}383${tab}.${tab}+${tab}.${tab}gene_id "e13a2_before_roi"; transcript_id "e13a2_before_roi";
bcr_abl1_e13_a2${tab}captured${tab}exon${tab}384${tab}445${tab}.${tab}+${tab}.${tab}gene_id "e13a2_roi"; transcript_id "e13a2_roi";
bcr_abl1_e13_a2${tab}captured${tab}exon${tab}446${tab}557${tab}.${tab}+${tab}.${tab}gene_id "e13a2_between_roi_rock"; transcript_id "e13a2_between_roi_rock";
bcr_abl1_e13_a2${tab}captured${tab}exon${tab}558${tab}582${tab}.${tab}+${tab}.${tab}gene_id "e13a2_rock"; transcript_id "e13a2_rock";
bcr_abl1_e13_a2${tab}captured${tab}exon${tab}583${tab}847${tab}.${tab}+${tab}.${tab}gene_id "e13a2_after_rock"; transcript_id "e13a2_after_rock";
bcr_abl1_e14_a2${tab}captured${tab}exon${tab}1${tab}383${tab}.${tab}+${tab}.${tab}gene_id "e14a2_before_roi_e13"; transcript_id "e14a2_before_roi_e13";
bcr_abl1_e14_a2${tab}captured${tab}exon${tab}384${tab}445${tab}.${tab}+${tab}.${tab}gene_id "e14a2_roi_e13"; transcript_id "e14a2_roi_e13";
bcr_abl1_e14_a2${tab}captured${tab}exon${tab}446${tab}457${tab}.${tab}+${tab}.${tab}gene_id "e14a2_between_roi_e13_e14"; transcript_id "e14a2_between_roi_e13_e14";
bcr_abl1_e14_a2${tab}captured${tab}exon${tab}458${tab}519${tab}.${tab}+${tab}.${tab}gene_id "e14a2_roi_e14"; transcript_id "e14a2_roi_e14";
bcr_abl1_e14_a2${tab}captured${tab}exon${tab}520${tab}632${tab}.${tab}+${tab}.${tab}gene_id "e14a2_between_roi_e14_rock"; transcript_id "e14a2_between_roi_e14_rock";
bcr_abl1_e14_a2${tab}captured${tab}exon${tab}633${tab}657${tab}.${tab}+${tab}.${tab}gene_id "e14a2_rock"; transcript_id "e14a2_rock";
bcr_abl1_e14_a2${tab}captured${tab}exon${tab}658${tab}922${tab}.${tab}+${tab}.${tab}gene_id "e14a2_after_rock"; transcript_id "e14a2_after_rock";
EOF


cat <<EOF > captured.fa
>bcr_abl1_e1_a2
CCTGGCCCCGCAGGTCCTACTCCCCCCGGAGTTTTGAGGATTGCGGAGGCGGCTATACCCCGGACTGCAGCTCCAATGAGAACCTCACCTCCAGCGAGGAGGACTTCTCCTCTGGCCAGTCCAGCCGCGTGTCCCCAAGCCCCACCACCTACCGCATGTTCCGGGACAAAAGCCGCTCTCCCTCGCAGAACTCGCAACAGTCCTTCGACAGCAGCAGTCCCCCCACGCCGCAGTGCCATAAGCGGCACCGGCACTGCCCGGTTGTCGTGTCCGAGGCCACCATCGTGGGCGTCCGCAAGACCGGGCAGATCTGGCCCAACGATGGCGAGGGCGCCTTCCATGGAGACGCAGAAGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGCACTGTATGATTTTGTGGCCAGTGGAGATAACACTCTAAGCATAACTAAAGGTGAAAAGCTCCGGGTCTTAGGCTATAATCACAATGGGGAATGGTGTGAAGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTACATCACGCCAGTCAACAGTCTGGAGAAACACTCCTGGTACCATGGGCCTGTGTCCCGCAATGCCGCTGAGTATCTGCTGAGCAGCGGGATCAATGGCAGCTTCTTGGTGCGTGAGAGTGAGAGCAGTCCTGGCCAGAGGTCCATCTCGCTGAGATACGAAGGGAGGGTGTACCATTACAGGAT
>bcr_abl1_e13_a2
TGGAGGCAGTGCCCAACATCCCCCTGGTGCCCGATGAGGAGCTGGACGCTTTGAAGATCAAGATCTCCCAGATCAAGAATGACATCCAGAGAGAGAAGAGGGCGAACAAGGGCAGCAAGGCTACGGAGAGGCTGAAGAAGAAGCTGTCGGAGCAGGAGTCACTGCTGCTGCTTATGTCTCCCAGCATGGCCTTCAGGGTGCACAGCCGCAACGGCAAGAGTTACACGTTCCTGATCTCCTCTGACTATGAGCGTGCAGAGTGGAGGGAGAACATCCGGGAGCAGCAGAAGAAGTGTTTCAGAAGCTTCTCCCTGACATCCGTGGAGCTGCAGATGCTGACCAACTCGTGTGTGAAACTCCAGACTGTCCACAGCATTCCGCTGACCATCAATAAGGAAGAAGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGCACTGTATGATTTTGTGGCCAGTGGAGATAACACTCTAAGCATAACTAAAGGTGAAAAGCTCCGGGTCTTAGGCTATAATCACAATGGGGAATGGTGTGAAGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTACATCACGCCAGTCAACAGTCTGGAGAAACACTCCTGGTACCATGGGCCTGTGTCCCGCAATGCCGCTGAGTATCTGCTGAGCAGCGGGATCAATGGCAGCTTCTTGGTGCGTGAGAGTGAGAGCAGTCCTGGCCAGAGGTCCATCTCGCTGAGATACGAAGGGAGGGTGTACCATTACAGGAT
>bcr_abl1_e14_a2
TGGAGGCAGTGCCCAACATCCCCCTGGTGCCCGATGAGGAGCTGGACGCTTTGAAGATCAAGATCTCCCAGATCAAGAATGACATCCAGAGAGAGAAGAGGGCGAACAAGGGCAGCAAGGCTACGGAGAGGCTGAAGAAGAAGCTGTCGGAGCAGGAGTCACTGCTGCTGCTTATGTCTCCCAGCATGGCCTTCAGGGTGCACAGCCGCAACGGCAAGAGTTACACGTTCCTGATCTCCTCTGACTATGAGCGTGCAGAGTGGAGGGAGAACATCCGGGAGCAGCAGAAGAAGTGTTTCAGAAGCTTCTCCCTGACATCCGTGGAGCTGCAGATGCTGACCAACTCGTGTGTGAAACTCCAGACTGTCCACAGCATTCCGCTGACCATCAATAAGGAAGATGATGAGTCTCCGGGGCTCTATGGGTTTCTGAATGTCATCGTCCACTCAGCCACTGGATTTAAGCAGAGTTCAAAAGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGCACTGTATGATTTTGTGGCCAGTGGAGATAACACTCTAAGCATAACTAAAGGTGAAAAGCTCCGGGTCTTAGGCTATAATCACAATGGGGAATGGTGTGAAGCCCAAACCAAAAATGGCCAAGGCTGGGTCCCAAGCAACTACATCACGCCAGTCAACAGTCTGGAGAAACACTCCTGGTACCATGGGCCTGTGTCCCGCAATGCCGCTGAGTATCTGCTGAGCAGCGGGATCAATGGCAGCTTCTTGGTGCGTGAGAGTGAGAGCAGTCCTGGCCAGAGGTCCATCTCGCTGAGATACGAAGGGAGGGTGTACCATTACAGGAT
EOF


human_fa=GRCh38.p13.genome.fa
human_gtf=gencode.v38.basic.annotation.gtf


wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_fa".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_gtf".gz

pigz --decompress *gz


cat $human_fa captured.fa sampletags.fa  > leukemia.fa
cat $human_gtf captured.gtf sampletags.gtf > leukemia.gtf

rm $human_fa $human_gtf
