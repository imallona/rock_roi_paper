#!/bin/bash
##
## Build STAR indices for chr22 chr9 bcr abl
##

cd /home/gmoro/indices

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

GENOME=bcr_abl_downsampled.fa
GTF=bcr_abl_dowsampled.gtf
NTHREADS=10
ID=bcr_abl_downsampled

## installs
cd /home/gmoro/indices

mkdir -p "$ID"

STAR --runMode genomeGenerate \
      --runThreadN "$NTHREADS" \
      --genomeDir "$ID" \
      --genomeFastaFiles "$GENOME" \
      --genomeSAindexNbases 12
