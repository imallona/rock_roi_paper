#!/bin/bash
##
## Analyzes the leukemia exp
##
## Started 6th Nov 2023
## Izaskun Mallona

# mkdir -p ~/src
# cd $_

snakemake --cores 50 \
          -s /home/imallona/rock_roi_method/main/Snakefile \
          --configfile 03_leukemia_simple.yaml -p --use-conda
