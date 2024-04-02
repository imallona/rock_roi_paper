# Experiment
Second mixing experiment

# Biological sample
eGFP positive mouse cells, tdTomato positive human cells

# Conditions
unmod: unmodified beads, standard library generation protocol
unmod_roi: unmodified beads, N1 primer, ROI primers
rock: RoCKseq modified beads, N1 primer
rockroi: RoCKseq modified beads, N1 primer, ROI primers
 
# .qmd files and corresponding figure panels

00_mixing_experiment_a-preprocessing-filtering.qmd: preliminary code for filtering and QC plots
00_mixing_experiment_b-dimred-ROIsensitivity.qmd: preliminary code for dimensionality reduction and detection of ROI sequences
00_mixing_experiment_c-figures.qmd: generation of filtered objects used for plots 
00_mixing_experiment_d-coverage_tracks.qmd: preliminary code for coverage plots
00_mixing_experiment_e_ROIs_other_genes.qmd: ROIseq peaks in off-target transcripts. Figure 4H-J
00_mixing_experiment_f-paper_figures.qmd: code for QC figures / detection of eGFP and tdTomato. Supp Figure 6D-F, Figure 3B, Supp Figure 6G, Figure 3E,Figure 4A-B, Supp Figure 7D-E, Figure 3H, Figure 3I, Supp Figure 7C, 
00_mixing_experiment_g-coverage_tracks.qmd: code for coverage tracks for eGFP and tdTomato. Figure 3F-G, Supp Figure 7A-B
00_mixing_h_RNAScope_figures.qmd: barplots indicating number of UMIs and RNAScope spots. Figure 3 J-K
00_mixing_i_wta_analysis.qmd: corrplots with comparison whole transcriptome. Figure 3C-D
00_mixing_j_mapping_stats.qmd: barplots with statistics derived from bamqc and rnaqc
00_mixing_k_qualimap_coverages.qmd: coverage across all transcripts derived from rnaqc. Figure 4C
00_mixing_l_gene_types.qmd: gene types detected in the TSO data and percentage of eGFP and tdTomato positive transcripts. Supp Figure 8A-C

# Figure panels and corresponding .qmd file
Figure 3B: 00_mixing_experiment_f-paper_figures.qmd
Figure 3C-D: 00_mixing_i_wta_analysis.qmd
Figure 3E: 00_mixing_experiment_f-paper_figures.qmd
Figure 3F-G: 00_mixing_experiment_g-coverage_tracks.qmd
Figure 3H: 00_mixing_experiment_f-paper_figures.qmd
Figure 3I: 00_mixing_experiment_f-paper_figures.qmd
Figure 3J-K: 00_mixing_h_RNAScope_figures.qmd
Figure 4A-B: 00_mixing_experiment_f-paper_figures.qmd
Figure 4C: 00_mixing_k_qualimap_coverages.qmd
Figure 4H-J: 00_mixing_experiment_e_ROIs_other_genes.qmd
Supp Figure 6D-F: 00_mixing_experiment_f-paper_figures.qmd
Supp Figure 6G: 00_mixing_experiment_f-paper_figures.qmd
Supp FIgure 7A-B: 00_mixing_experiment_g-coverage_tracks.qmd
Supp Figure 7C: 00_mixing_experiment_f-paper_figures.qmd
Supp Figure 7D-E: 00_mixing_experiment_f-paper_figures.qmd
Supp Figure 8A-B: 00_mixing_l_gene_types.qmd
Supp Figure 8C: 00_mixing_l_gene_types.qmd
