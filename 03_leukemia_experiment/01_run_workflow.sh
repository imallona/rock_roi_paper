#!/bin/bash
##
## Analyzes the leukemia exp
##
## Started 6th Nov 2023
## Izaskun Mallona

# mkdir -p ~/src
# cd $_


snakemake --cores 20 \
          -s /home/imallona/src/rock_roi_method/main/Snakefile \
          --configfile ./leukemia.yaml -p



