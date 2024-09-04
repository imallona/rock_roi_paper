#!/bin/bash
##

human_fa=GRCh38.p13.genome.fa
human_gtf=gencode.v38.basic.annotation.gtf

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_fa".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_gtf".gz
pigz --decompress *gz

singularity exec -e -B `pwd` -B ~/STAR-Fusion/ctat-genome-lib-builder \
	~/star-fusion.v1.13.0.simg \
	~/STAR-Fusion//ctat-genome-lib-builder/prep_genome_lib.pl \
	--genome_fa "$human_fa" \
	--gtf "$human_gtf" \
	--fusion_annot_lib ~/STAR-Fusion-Tutorial/CTAT_HumanFusionLib.mini.dat.gz \
	--dfam_db human \
	--pfam_db current

rm "$human_fa"
rm "$human_gtf"
