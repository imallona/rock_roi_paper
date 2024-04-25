#!/bin/bash
##
##

#cd ~/first_mixing_experiment/align_wta/unmod_n1

echo "mapped reads WTA"
#samtools view Aligned.sortedByCoord.out.bam |cut -f1 |sort |uniq |wc -l 

echo "in cell barcodes WTA"
#samtools view -@ 20 Aligned.sortedByCoord.out.bam | grep -vP "CB:Z:-\t" | cut -f1 | sort | uniq | wc -l

echo "multimappers WTA"
#samtools view -@ 20 Aligned.sortedByCoord.out.bam | grep -vP "CB:Z:-\t" | grep -vP '^@|NH:i:1\b' |  cut -f1 | sort | uniq | wc -l

echo "unique WTA"
#samtools view -@ 20 Aligned.sortedByCoord.out.bam |  grep -vP "CB:Z:-\t" | grep -P '^@|NH:i:1\b'| cut -f1 | sort | uniq | wc -l

echo "unique in genes WTA"
#samtools view -@ 20  Aligned.sortedByCoord.out.bam| grep -vP "CB:Z:-\t" | grep -P '^@|NH:i:1\b' | grep -vP "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 

echo "unique not in genes WTA"
#samtools view -@ 20  Aligned.sortedByCoord.out.bam| grep -vP "CB:Z:-\t" | grep -P '^@|NH:i:1\b' | grep -P "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 


cd ~/first_mixing_experiment/align_tso/unmod_n1
ln -s ~/first_mixing_experiment/align_wta/unmod_n1/Solo.out/Gene/filtered/barcodes.tsv

echo "mapped reads TSO"
#samtools view Aligned.sortedByCoord.out.bam |cut -f1 |sort |uniq |wc -l 

echo "in cell barcodes TSO"
#samtools view -@ 20 Aligned.sortedByCoord.out.bam | grep -vP "CB:Z:-\t" | cut -f1 | sort | uniq | wc -l

echo "TSO reads in WTA filtered barcodes"

samtools view -D CB:barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam | cut -f1 | sort | uniq | wc -l

echo "multimappers TSO"
samtools view -D CB:barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam | grep -vP '^@|NH:i:1\b' |  cut -f1 | sort | uniq | wc -l

echo "unique TSO"
samtools view -D CB:barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam | grep -P '^@|NH:i:1\b'| cut -f1 | sort | uniq | wc -l

echo "unique in genes TSO"
samtools view -D CB:barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam| grep -P '^@|NH:i:1\b' | grep -vP "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 

echo "unique not in genes TSO"
samtools view -D CB:barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam| grep -P '^@|NH:i:1\b' | grep -P "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 
