# Experiment
Detection of Pdgfra junctions  

# Biological sample
Epithelial cells and fibroblasts from murine colon  

# Conditions
**unmod:** unmodified beads, standard library generation protocol  
**rockroi:** RoCKseq modified beads, N1 primer, ROI primers  

**unimodal:**  dT data from unmod and rockroi samples  
**multimodal:** dT and TSO data from rockroi sample  
 
# .qmd files and corresponding figure panels
**01_pdgfra_experiment_a-preprocessing-filtering.qmd:** generation of filtered datasets  
**01_pdgfra_experiment_b-dimred-clustering.qmd:** dimensionality reduction (without annotation)  
**01_pdgfra_experiment_c-annotation.qmd:** annotation of clusters  
**01_pdgfra_experiment_d-coverage_tracks.qmd:** preliminary code for Pdgfra coverage tracks  
**01_pdgfra_experiment_d-tuft.qmd:** annotation of Tuft and enteroendocrine clusters  
**01_pdgfra_experiment_e-bring-annotations-together.qm:** bringing all annotations together  
**01_pdgfra_experiment_e-unimodal_vs_multimodal.qmd:** preliminary code for comparison unimodal vs multimodal  
**01_pdgfra_experiment_f-junction-counts.qmd:** counting reads spanning exon junctions
**01_pdgfra_experiment_f-unimodal_multimodal_paper_figures.qmd:** comparison of unimodal vs multimodal data. Supp Figure 10J-K
**01_pdgfra_experiment_g-paper_figures_unmod_multimodal.qmd:** figures for paper. Figure 5B-G, Supp Figure 10D-I, Supp Figure 11A-B
**01_pdgfra_experiment_h-downsampling_files.qmd:** generation of downsampled data objects. 
**01_pdgfra_experiment_r-pdgfra_only_junctions.qmd:** QC figures for downsampled data. Supp Figure 10A-B
