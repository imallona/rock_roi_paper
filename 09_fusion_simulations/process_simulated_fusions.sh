#!/bin/bash
##
## Evaluates a fusion calling procedure on simulated data with ground truths
##
## Started 11th Nov 2024
##
## Izaskun Mallona

NTHREADS=10

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

# micromamba create fusion
# micromamba activate fusion
# micromamba install -f env/fusion_conda_env.txt


echo 'Extract UMIs/CBs'

#regex='^(?P<discard_1>.{0,3})(?P<cell_1>.{9})(?P<discard_2>GTGA)(?P<cell_2>.{9})(?P<discard_3>GACA)(?P<cell_3>.{9})(?P<umi_1>.{8}).*'
regex='^(?P<discard_1>.{0,3})(?P<cell_1>.{9})(?P<cell_2>[GTGA|AATG]{4}){s<=1}(?P<cell_3>.{9})(?<cell_4>[GACA|CCAC]{4}){s<=1}(?P<cell_5>.{9})(?P<umi_1>.{8}).*'


mkdir -p out log

## so we label the cdnas with the CB and UMI
umi_tools extract --extract-method=regex \
          --stdin=./data/fusion_simulations_cbumi.fq \
          --read2-in=./data/fusion_simulations_cdna.fq \
          --read2-out=./out/labelled_umis_cdna.fq.gz \
          --bc-pattern="$regex" --log=log/umitools.log --stdout out/processed.fastq.gz

zcat ./out/labelled_umis_cdna.fq.gz | tail

echo 'Scan for motifs (with mismatches, against a reference)'

## five mm
zcat out/labelled_umis_cdna.fq.gz | \
    seqkit locate  \
           --max-mismatch 5 \
           --pattern-file data/reference_fusions.fa \
           -j $NTHREADS | gzip -c > ./out/fusion_locate_out_mm5.txt.gz

zcat ./out/fusion_locate_out_mm5.txt.gz | head 


# matched
zcat out/fusion_locate_out_mm5.txt.gz | cut -f2 | grep -f - data/fusion_simulations_cdna.fq

# missing
zcat out/fusion_locate_out_mm5.txt.gz | cut -f2 | grep -v -f - data/fusion_simulations_cdna.fq | grep "^@"

## with regular expressions

zcat out/labelled_umis_cdna.fq.gz | \
    seqkit locate  \
           --use-regexp \
           --pattern-file data/reference_fusions_regex.fa \
           -j $NTHREADS | gzip -c > ./out/fusion_locate_out_reg.txt.gz

zcat ./out/fusion_locate_out_reg.txt.gz | head

# matched
zcat out/fusion_locate_out_reg.txt.gz | cut -f2 | grep -f - data/fusion_simulations_cdna.fq

# missing
zcat out/fusion_locate_out_reg.txt.gz | cut -f2 | grep -v -f - data/fusion_simulations_cdna.fq | grep "^@"

# ## deduplicate by CB, UMI, pattern and start - here or in R?


# # quick dirty check cell-lines
# giulia_cdna=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/downsampled_R2_5M.fastq.gz
# giulia_cbumi=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/downsampled_R1_5M.fastq.gz

# umi_tools extract --extract-method=regex \
#           --stdin="$giulia_cbumi" \
#           --read2-in="$giulia_cdna" \
#           --read2-out=./out/cellline_labelled_umis_cdna.fq.gz \
#           --bc-pattern="$regex" --log=log/cellline_umitools.log --stdout out/celllines.fq.gz


# pigz -dc -p 5 ./out/cellline_labelled_umis_cdna.fq.gz | \
#     seqkit locate  \
#            --use-regexp \
#            -G \
#            --pattern-file data/reference_fusions.fa \
#            -j 40 | gzip -c > ./out/cellline_regex_greedy.txt.gzip


