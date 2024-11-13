#!/bin/bash
##
## Evaluates a fusion calling procedure on simulated data with ground truths
##
## Started 11th Nov 2024
##
## Izaskun Mallona

NTHREADS=30

# enhanced beads regex
regex='^(?P<discard_1>.{0,3})(?P<cell_1>.{9})(?P<cell_2>[GTGA|AATG]{4}){s<=1}(?P<cell_3>.{9})(?<cell_4>[GACA|CCAC]{4}){s<=1}(?P<cell_5>.{9})(?P<umi_1>.{8}).*'


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
# micromamba install -f env/fusion_conda_env.yaml

echo 'Extract UMIs/CBs'

mkdir -p out log

## edit here
run_id="cell_lines_full"
clines_cdna=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/331131_1-Cell_lines_50_50_S1_R2_001.fastq.gz
clines_cbumi=/home/gmoro/test_leukemia_downsampled_cell_line_experiment/331131_1-Cell_lines_50_50_S1_R1_001.fastq.gz

umi_tools extract --extract-method=regex \
          --stdin="$clines_cbumi" \
          --read2-in="$lines_cdna" \
          --read2-out=./out/"$run_id"_labelled_umis_cdna.fq.gz \
          --bc-pattern="$regex" --log=log/"$run_id"_umitools.log --stdout out/"$run_id".fq.gz


pigz -dc -p "$NTHREADS" ./out/cellline_labelled_umis_cdna.fq.gz | \
    seqkit locate  \
           --use-regexp \
           --pattern-file data/reference_fusions_regex.fa \
           -j "$NTHREADS" | pigz -c -p "$NTHREADS" > ./out/"$run_id"_regex.txt.gz
