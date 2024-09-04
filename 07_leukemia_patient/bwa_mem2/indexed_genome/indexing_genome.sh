#!/bin/bash
##
## Generate combined genome for human and leukemia rois

ln -s ~/leukemia_bwamem2/genome/combined.fa combined.fa

~/leukemia_bwamem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p indexed combined.fa
