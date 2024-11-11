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
Content matches ./data/reference_fusions.fa

wt_ABL1_1      GTTTTTGTGGAACATG/AAGCCCTTCAGCGGCC
wt_ABL1_2      AGCTGTTATCTGGAAG/AAGCCCTTCAGCGGCC
wt_BCR_e1a2    TTCCATGGAGACGCAG/ATGGCTCGTTCGGAAC
wt_BCR_e13a2   ACCATCAATAAGGAAG/ATGATGAGTCTCCGGG
wt_BCR_e14a2   TTTAAGCAGAGTTCAA/ATCTGTACTGCACCCT
control1       GGGGGGGGGGGGGGGG/CCCCCCCCCCCCCCCC
control2       GGGGGGGGGGGGGGGG/GGGGGGGGGGGGGGGG
BCRABL_e1a2    TTCCATGGAGACGCAG/AAGCCCTTCAGCGGCC
BCRABL_1e13a2  ACCATCAATAAGGAAG/AAGCCCTTCAGCGGCC
BCRABL1_e14a2  TTTAAGCAGAGTTCAA/AAGCCCTTCAGCGGCC 


EOF


echo 'Uncomment to install using conda/similar the environment env/fusion_conda_env.txt'

# micromamba create fusion
# micromamba activate fusion
# micromamba install -f env/fusion_conda_env.txt

echo 'Scan for motifs (with mismatches, against a reference)'

mkdir -p out

cat ./data/fusion_simulations_cdna.fq | \
    seqkit locate  \
           --max-mismatch 2 \
           --pattern-file data/reference_fusions.fa \
           -j $NTHREADS > .out/fusion_simulations_seqkit_locate.out
       
# This cannot be correct, the simulations are not ok because in line 3 we get a fused simulation aligning to
## a nonfused reference?
# seqID	patternName	pattern	strand	start	end	matched
# major_13_fusion__overlaps_fusion=major__bcr=fused__abl=fused__expb=1__exps=1	BCRABL1_e14a2	TTTAAGCAGAGTTCAAAAGCCCTTCAGCGGCC 	+	76	108	TTTAAGCAGAGTTCAAAAGCCCTTCAGCGGCCA
# major_13_fusion__overlaps_fusion=major__bcr=fused__abl=fused__expb=1__exps=1	wt_BCR_e13a2	ACCATCAATAAGGAAGATGATGAGTCTCCGGG	+	1	32	ACCATCAATAAGGAAGATGATGAGTCTCCGGG



# first, grep the reads with the fusion motifs
## https://github.com/fulcrumgenomics/fqgrep

## or perhaps with seqkit locate
## https://bioinf.shenwei.me/seqkit/usage/#grep

# then, collapse by umi seqkit rmdup?

# rather umitools the fastqs, to add the UMI to each read name?
# umi_tools dedup -I example.bam --output-stats=deduplicated -S deduplicated.bam
