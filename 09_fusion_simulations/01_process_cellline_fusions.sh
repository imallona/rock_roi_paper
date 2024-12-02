#!/bin/bash
##
## Evaluates a fusion calling procedure on simulated data with ground truths
##
## Started 11th Nov 2024
##
## Izaskun Mallona

NTHREADS=30
run_id="cell_lines_full"
clines_cdna=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/331131_1-Cell_lines_50_50_S1_R2_001.fastq.gz
star_bam=/home/gmoro/cell_lines_leukemia/rock_roi_repo_for_cell_lines/output/align_tso/leukemia_cell_lines/Aligned.sortedByCoord.out.bam



:<<EOF
Expected reads follow this pattern, being the `/` the splicing site or fusion site
  apparently only 16 nt on each side of the fusion are expected (Giulia dixit, something about
  rois not being efficient)
The actual content is available as ./data/reference_fusions.fa

ENST00000372348.7						GTTTTTGTGGAACATG/AAGCCCTTCAGCGGCC # Human ABL1 across fusion (no fusion)
ENST00000318560.6 						AGCTGTTATCTGGAAG/AAGCCCTTCAGCGGCC # Human ABL1 across fusion (no fusion)
ENST00000359540.7_BCRABLe1a2					TTCCATGGAGACGCAG/ATGGCTCGTTCGGAAC # Human BCR across fusion at position of e1a2 (no fusion)
ENST00000359540.7_BCRABLe13a2 					ACCATCAATAAGGAAG/ATGATGAGTCTCCGGG # Human BCR across fusion at position of e13a2 (no fusion)
ENST00000359540.7_BCRABLe14a2 					TTTAAGCAGAGTTCAA/ATCTGTACTGCACCCT # Human BCR across fusion at position of e14a2 (no fusion)
BCRABLe1a2 							TTCCATGGAGACGCAG/AAGCCCTTCAGCGGCC # BCR-ABL1 fusion e1a2
BCRABL1e13a2 							ACCATCAATAAGGAAG/AAGCCCTTCAGCGGCC # BCR-ABL1 fusion e13a2
BCRABLe14a2							TTTAAGCAGAGTTCAA/AAGCCCTTCAGCGGCC # BCR-ABL1 fusion e14a2
Reciprocal_ENST00000372348.7_ENST00000359540.7_BCRABLe1a2	GTTTTTGTGGAACATG/ATGGCTCGTTCGGAAC # Reciprocal ABL1-BCR e1a2 fusion
Reciprocal_ENST00000318560.6_ENST00000359540.7_BCRABLe1a2	AGCTGTTATCTGGAAG/ATGGCTCGTTCGGAAC # Reciprocal ABL1-BCR e1a2 fusion
Reciprocal_ENST00000372348.7_ENST00000359540.7_BCRABLe13a2	GTTTTTGTGGAACATG/ATGATGAGTCTCCGGG # Reciprocal ABL1-BCR e13a2 fusion
Reciprocal_ENST00000318560.6_ENST00000359540.7_BCRABLe13a2	AGCTGTTATCTGGAAG/ATGATGAGTCTCCGGG # Reciprocal ABL1-BCR e13a2 fusion
Reciprocal_ENST00000372348.7_ENST00000359540.7_BCRABLe14a2 	GTTTTTGTGGAACATG/ATCTGTACTGCACCCT # Reciprocal ABL1-BCR e14a2 fusion
Reciprocal_ENST00000318560.6_ENST00000359540.7_BCRABLe14a2 	AGCTGTTATCTGGAAG/ATCTGTACTGCACCCT # Reciprocal ABL1-BCR e14a2 fusion

To allow for deletions at the fusion site we also search for regexes

TTTAAGCAGAGTTCAA/AAGCCCTTCAGCGGCC with up to six deletions or six insertions symmetrical to the junction
                                  that is, up to 12 "free" nucleotides centering 6:6 from the junction viewpoint
 
EOF

echo 'Uncomment to install using conda/similar the environment env/fusion_conda_env.txt'

# conda install micromamba
# mamba init
# micromamba create -n fusion
# micromamba activate fusion
# micromamba install -f env/fusion_conda_env.yaml

mkdir -p out log

## cannot run the UMI deduplication with UMI tools on the .bam file from STARsolo as we also need the unmapped reads

## run scan on fastq file with cDNA

seqkit locate $clines_cdna \
        --use-regexp \
        --pattern-file ./data/reference_fusions_regex.fa \
        -j $NTHREADS > ./out/"$run_id"_fusion_locate_out_reg.txt

## from .bam file generate .txt file with read id, cb and ub

## need to do it separate for mapped and unmapped as will have different columns

## also need to grep for the read ids which are in the seqkit 

cut -f1 ./out/"$run_id"_fusion_locate_out_reg.txt > read_ids.txt

samtools view -H -@ $NTHREADS $star_bam > header.txt

samtools view -@ $NTHREADS $star_bam | grep -f read_ids.txt | cat header.txt - | samtools view -Sbh > ./out/subsetted.bam

echo 'Mapped'

echo -e 'Read_id\tCB\tUB' > ./out/"$run_id"_reads_with_cb_ub.txt

samtools view -@ $NTHREADS -F 4 ./out/subsetted.bam | cut -f1,27,28 >> ./out/"$run_id"_reads_with_cb_ub.txt

echo 'Unmapped'

samtools view -@ $NTHREADS -f 4 ./out/subsetted.bam | cut -f1,22,23 >> ./out/"$run_id"_reads_with_cb_ub.txt

# need to sort and unique them 

sort -k 1,1 ./out/"$run_id"_reads_with_cb_ub.txt -u > ./out/sorted_"$run_id"_reads_with_cb_ub.txt

rm header.txt read_ids.txt ./out/subsetted.bam

wc -l ./out/sorted_"$run_id"_reads_with_cb_ub.txt
wc -l ./out/"$run_id"_fusion_locate_out_reg.txt

## append the information to the seqtk file

## need to sort the two files and remove the header since it is not part of the information

tail -n+2 ./out/"$run_id"_fusion_locate_out_reg.txt | sort -k 1,1 >> ./out/"$run_id"_sorted_fusion_locate_out_reg.txt
tail -n+2 ./out/sorted_"$run_id"_reads_with_cb_ub.txt | sort -k 1,1 >> ./out/"$run_id"_sorted_reads_with_cb_ub.txt

join -1 1 -2 1 -e '-' ./out/"$run_id"_sorted_fusion_locate_out_reg.txt ./out/"$run_id"_sorted_reads_with_cb_ub.txt > ./out/"$run_id"_annotated_fusion_locate_out_reg.txt

rm ./out/"$run_id"_reads_with_cb_ub.txt ./out/"$run_id"_sorted_reads_with_cb_ub.txt ./out/"$run_id"_sorted_fusion_locate_out_reg.txt ./out/sorted_"$run_id"_reads_with_cb_ub.txt

## removing all alignments with no UB and CB to make the file smaller

less ./out/"$run_id"_annotated_fusion_locate_out_reg.txt | grep -v 'CB:Z:-' | grep -v 'UB:Z:-' > ./out/"$run_id"_cb_ub_annotated_fusion_locate_out_reg.txt

rm ./out/"$run_id"_annotated_fusion_locate_out_reg.txt
