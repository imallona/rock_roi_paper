Note: the code in this repository was used for initial analysis but is no longer implemented, please see the folder 09_fusion_simulations

# Experiment

Leukemia patient sample

# Files and folders

Snakefile: used for running rock_roi_method (see https://github.com/imallona/rock_roi_method). The Snakefile was modified to add additional folders in which tmp files are stored and later deleted.<br />
config.yaml: used for running rock_roi_method (see https://github.com/imallona/rock_roi_method) for the patient sample. <br />

qmd: files for analysis of count tables and generation of plots<br />
bwa_aln_full_fastq_file: workflow used for detection of fusions in TSO data<br />

Installs neede to run the workflow:
- STAR: version 2.7.10b was used for the analysis
- bwa: bwa 0.7.18 was used for the analysis
- samtools: samtools 1.21 was used for the analysis <br />
The installs can be done with the run.sh file

Files needed to run the workflow (on real data): 
- custom .fa file: contains the cDNAs covering the fusion region to be detected. The file is found in the simulation_patient_cell_lines folder (BCR_ABL.fa) 
- STAR indexed genome: used for running STAR fusion
- TSO .bam file generated with the rock_roi_method <br />
The script run.sh is used to define the parameters used in the workflow.sh script and to run the script
