#!/bin/bash
##
## Generate combined genome for human and leukemia rois

WD=/home/gmoro/mapping_leukemia_data/genome
cd $WD

human_fa=GRCh38.p13.genome.fa
human_gtf=gencode.v38.basic.annotation.gtf

captured_gtf=roi_bcr_abl.gtf

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_fa".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_gtf".gz

pigz --decompress *gz


cat $human_gtf $captured_gtf > combined.gtf

mv $human_fa combined.fa

rm $human_gtf
