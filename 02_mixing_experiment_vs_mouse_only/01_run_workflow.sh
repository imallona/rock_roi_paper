#!/bin/bash
##
## Analyzes the mixing exp vs the mouse genome only
##
## Started 27th Oct 2023
## Izaskun Mallona

mkdir -p ~/src
cd $_

# git clone  git@github.com:imallona/rock_roi_method.git

ln -s ~/src/rock_roi_method/main .

## dirty workaround to subset `chr5` too (harcoded within the Snakefile)
snakemake --cores 20 \
          -s /home/imallona/src/rock_roi_method/main/Snakefile \
          --configfile ./config.yaml -p



