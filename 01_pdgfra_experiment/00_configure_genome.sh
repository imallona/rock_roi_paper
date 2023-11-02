#!/bin/bash
##
## Generates a mouse + human + gfp + mouse genome to analyze the pdfgfra experiment
##
## Started 24th Oct 2023
## Izaskun Mallona


WD=~/ebrunner_spectral/pdgfra
mkdir -p $WD/genomes
cd $WD/genomes


tab="$(printf '\t')"

## the 'exon_id' refers to splicing junctions I think (?)
cat <<EOF > captured.gtf
egfp${tab}captured${tab}exon${tab}1${tab}378${tab}.${tab}+${tab}.${tab}gene_id "h2b"; transcript_id "h2b";
egfp${tab}captured${tab}exon${tab}379${tab}396${tab}.${tab}+${tab}.${tab}gene_id "linker"; transcript_id "linker";
egfp${tab}captured${tab}exon${tab}397${tab}708${tab}.${tab}+${tab}.${tab}gene_id "egfp_before_roi"; transcript_id "egfp_before_roi";
egfp${tab}captured${tab}exon${tab}709${tab}769${tab}.${tab}+${tab}.${tab}gene_id "roi_egfp"; transcript_id "roi_egfp";
egfp${tab}captured${tab}exon${tab}770${tab}1049${tab}.${tab}+${tab}.${tab}gene_id "between_roi"; transcript_id "between_roi";
egfp${tab}captured${tab}exon${tab}1050${tab}1074${tab}.${tab}+${tab}.${tab}gene_id "capture_egfp"; transcript_id "capture_egfp";
egfp${tab}captured${tab}exon${tab}1075${tab}1113${tab}.${tab}+${tab}.${tab}gene_id "3_to_capture"; transcript_id "3_to_capture";
egfp${tab}captured${tab}exon${tab}1114${tab}1425${tab}.${tab}+${tab}.${tab}gene_id "3_utr_egfp"; transcript_id "3_utr_egfp"
chr5${tab}captured${tab}exon${tab}75166756${tab}75166771${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_4"; exon_id "roi_4_4";
chr5${tab}captured${tab}exon${tab}75167837${tab}75167882${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_4"; exon_id "roi_4_5";
chr5${tab}captured${tab}exon${tab}75167954${tab}75167967${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_5"; exon_id "roi_5_5";
chr5${tab}captured${tab}exon${tab}75170495${tab}75170544${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_5"; exon_id "roi_5_6";
chr5${tab}captured${tab}exon${tab}75170652${tab}75170666${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_6"; exon_id "roi_6_6";
chr5${tab}captured${tab}exon${tab}75170762${tab}75170808${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_6"; exon_id "roi_6_7";
chr5${tab}captured${tab}exon${tab}75180257${tab}75180272${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi 13"; exon_id "roi_13_13";
chr5${tab}captured${tab}exon${tab}75180999${tab}75181044${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_13"; exon_id "roi_13_14";
chr5${tab}captured${tab}exon${tab}75181096${tab}75181110${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_14"; exon_id "roi_14_14";
chr5${tab}captured${tab}exon${tab}75181522${tab}75181569${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_14"; exon_id "roi_14_15";
chr5${tab}captured${tab}exon${tab}75181661${tab}75181675${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_15"; exon_id "roi_15_15";
chr5${tab}captured${tab}exon${tab}75182976${tab}75183022${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_15"; exon_id "roi_15_16";
chr5${tab}captured${tab}exon${tab}75183128${tab}75183142${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_16"; exon_id "roi_16_16";
chr5${tab}captured${tab}exon${tab}75185514${tab}75185560${tab}.${tab}+${tab}.${tab}gene_id "ENSMUSG00000029231.15"; transcript_id "roi_16"; exon_id "roi_16_17";
EOF

cat <<EOF > captured.fa
>egfp
ATGCCAGAGCCAGCGAAGTCTGCTCCCGCCCCGAAAAAGGGCTCCAAGAAGGCGGTGACTAAGGCGCAGA
AGAAAGGCGGCAAGAAGCGCAAGCGCAGCCGCAAGGAGAGCTATTCCATCTATGTGTACAAGGTTCTGAA
GCAGGTCCACCCTGACACCGGCATTTCGTCCAAGGCCATGGGCATCATGAATTCGTTTGTGAACGACATT
TTCGAGCGCATCGCAGGTGAGGCTTCCCGCCTGGCGCATTACAACAAGCGCTCGACCATCACCTCCAGGG
AGATCCAGACGGCCGTGCGCCTGCTGCTGCCTGGGGAGTTGGCCAAGCACGCCGTGTCCGAGGGTACTAA
GGCCATCACCAAGTACACCAGCGCTAAGGATCCACCGGTCGCCACCATGGTGAGCAAGGGCGAGGAGCTG
TTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCG
GCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCC
CGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCAC
ATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCA
AGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGA
GCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGC
CACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACA
TCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCT
GCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCAC
ATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCGG
CCAATTCCAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTCTGGTGTCAGGGGATCCCCCCTCGCTGA
TCAGCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCC
TGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAAGAAATTGCATCGCATTGTCTGAGTAGGTG
TCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCAT
GCTGGGGATGCGGTGGGCTCTATGG
EOF


# download a mouse genome and its GTF

mouse_fa=GRCm38.p6.genome.fa
mouse_gtf=gencode.vM25.annotation.gtf

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"$mouse_fa".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"$mouse_gtf".gz

pigz --decompress *gz

cat $mouse_fa captured.fa  > pdgfra.fa
cat $mouse_gtf captured.gtf  > pdgfra.gtf

rm $mouse_fa $mouse_gtf
