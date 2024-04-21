#!/bin/bash
##
##

cd /home/imallona/ebrunner_spectral/mixing/align_wta/mixing_rockroi/

echo "mapped reads WTA"
samtools view Aligned.sortedByCoord.out.bam |cut -f1 |sort |uniq |wc -l 

echo "in cell barcodes WTA"
samtools view -@ 20 Aligned.sortedByCoord.out.bam | grep -vP "CB:Z:-\t" | cut -f1 | sort | uniq | wc -l

echo "multimappers WTA"
samtools view -@ 20 Aligned.sortedByCoord.out.bam | grep -vP "CB:Z:-\t" | grep -vP '^@|NH:i:1\b' |  cut -f1 | sort | uniq | wc -l

echo "unique WTA"
samtools view -@ 20 Aligned.sortedByCoord.out.bam |  grep -vP "CB:Z:-\t" | grep -P '^@|NH:i:1\b'| cut -f1 | sort | uniq | wc -l

echo "unique in genes WTA"
samtools view -@ 20  Aligned.sortedByCoord.out.bam| grep -vP "CB:Z:-\t" | grep -P '^@|NH:i:1\b' | grep -vP "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 

echo "unique not in genes WTA"
samtools view -@ 20  Aligned.sortedByCoord.out.bam| grep -vP "CB:Z:-\t" | grep -P '^@|NH:i:1\b' | grep -P "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 




cd /home/imallona/ebrunner_spectral/mixing/align_tso/mixing_rockroi/

echo "mapped reads TSO"
samtools view Aligned.sortedByCoord.out.bam |cut -f1 |sort |uniq |wc -l 

echo "in cell barcodes TSO"
samtools view -@ 20 Aligned.sortedByCoord.out.bam | grep -vP "CB:Z:-\t" | cut -f1 | sort | uniq | wc -l

echo "TSO reads in WTA filtered barcodes"
samtools view -D CB:/home/imallona/ebrunner_spectral/mixing/align_wta/mixing_rockroi/Solo.out/Gene/filtered/barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam | cut -f1 | sort | uniq | wc -l

echo "multimappers TSO"
samtools view -D CB:/home/imallona/ebrunner_spectral/mixing/align_wta/mixing_rockroi/Solo.out/Gene/filtered/barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam | grep -vP '^@|NH:i:1\b' |  cut -f1 | sort | uniq | wc -l

echo "unique TSO"
samtools view -D CB:/home/imallona/ebrunner_spectral/mixing/align_wta/mixing_rockroi/Solo.out/Gene/filtered/barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam | grep -P '^@|NH:i:1\b'| cut -f1 | sort | uniq | wc -l

echo "unique in genes TSO"
samtools view -D CB:/home/imallona/ebrunner_spectral/mixing/align_wta/mixing_rockroi/Solo.out/Gene/filtered/barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam| grep -P '^@|NH:i:1\b' | grep -vP "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 

echo "unique not in genes TSO"
samtools view -D CB:/home/imallona/ebrunner_spectral/mixing/align_wta/mixing_rockroi/Solo.out/Gene/filtered/barcodes.tsv -@ 20 Aligned.sortedByCoord.out.bam| grep -P '^@|NH:i:1\b' | grep -P "gn:Z:-\t" | cut -f1 | sort | uniq | wc -l 

