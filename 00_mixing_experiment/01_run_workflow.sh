#!/bin/bash
##
## Analyzes the mixing experiment fastqs from scratch

mkdir -p ~/src
cd $_

# git clone  git@github.com:imallona/rock_roi_method.git

# cd ~/src/rock_roi_method/main
cd /home/imallona/src/rock_roi_paper/00_mixing_experiment

snakemake --cores 10 --configfile config.yaml \
          -s /home/imallona/src/rock_roi_method/main/Snakefile -p


# docker build . -t rock

# docker run -it --entrypoint /bin/bash rock


