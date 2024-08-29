#!/bin/bash
##
## Generate combined genome for human and leukemia rois

WD=/home/gmoro/leukemia_bwamem2/genome
cd $WD

# just using the transcript sequences as bwamem2 is not splicing aware

human_fa=gencode.v38.pc_transcripts.fa
captured_fa=BCR_ABL.fa

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_fa".gz

pigz --decompress *gz

# getting the ABL1 sequence

grep "ABL1" "$human_fa"| sed 's/>//g' > ABL1_human.txt
grep "BCR" "$human_fa"| sed 's/>//g' >  BCR_human.txt

touch ABL1_human.fa
touch BCR_human.fa

for i in $(cat ABL1_human.txt); do
    samtools faidx $human_fa "$i" >> ABL1_human.fa
done

for i in $(cat BCR_human.txt) ; do
    samtools faidx $human_fa "$i"  >> BCR_human.fa
done

cat "$captured_fa" ABL1_human.fa BCR_human.fa > combined.fa

