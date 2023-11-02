#!/bin/bash
##
## Analyzes the PDGFRA experiment fastqs from scratch
##
## Started 24th Oct 2023
## Izaskun Mallona

# mkdir -p ~/src
# cd $_

# git clone  git@github.com:imallona/rock_roi_method.git


cd /home/imallona/src/rock_roi_paper/01_pdgfra_experiment

## dirty workaround to subset `chr5` too (harcoded within the Snakefile)
snakemake --cores 20 \
          -s /home/imallona/src/rock_roi_method/main/Snakefile \
          --configfile config.yaml -p


# docker build . -t rock

# docker run -it --entrypoint /bin/bash rock


