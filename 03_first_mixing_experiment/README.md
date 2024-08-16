# Experiment
First  mixing experiment

# Biological sample
eGFP positive mouse cells, tdTomato positive human cells

# Conditions
**unmod:** unmodified beads, standard library generation protocol  
**unmod_n1:** unmodified beads, N1 primer  
**rock:** RoCKseq modified beads, N1 primer  
**rockroi:** RoCKseq modified beads, N1 primer, ROI primers  
 
# .qmd files and corresponding figure panels
**03_first_mixing_experiment_c_figures.qmd:** generation of filtered objects used for plots  
**03_first_mixing_experiment_d-paper_figures.qmd:** code for QC figures / detection of eGFP and tdTomato. Figure 2C-D, Supp Figure 7B-C  
**03_first_mixing_experiment_e-tso_data.qmd:** TSO motif search. Figure 4D
**03_first_mixing_experiment_f-qualimap_coverages.qmd:** coverage across all transcripts derived from rnaqc. Supp Figure 8F  
**03_first_mixing_experiment_g-mapping_stats_tso_data.qmd:** barplots with statistics derived from bamqc and rnaqc  
**03_first_mixing_experiment_h-mito_coverages.qmd:** coverage plots for mitochondrial transcripts with TSO motif. Figure 4E-F  
**03_first_mixing_experiment_i-wta_analysis.qmd:** corrplots with comparison whole transcriptome. Figure 2E-F
**03_first_mixing_experiment_j-downsampling_count_tables.qmd:** generation of downsampled objects
**03_first_mixing_experiment_l-downsampled_figures.qmd:** QC figures using downsampled data. Supp Figure 4D-F
