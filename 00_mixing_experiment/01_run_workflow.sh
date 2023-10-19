#!/bin/bash
##
## Analyzes the mixing experiment fastqs from scratch

mkdir -p ~/src
cd $_

git clone  git@github.com:imallona/rock_roi_method.git

cd ~/src/rock_roi_method/main

snakemake --cores 50 --configfile /home/imallona/src/rock_roi_paper/00_mixing_experiment/mixing_conf.yaml


# docker build . -t rock

# docker run -it --entrypoint /bin/bash rock


