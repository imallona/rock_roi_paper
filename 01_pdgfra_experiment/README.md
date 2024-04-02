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
**01_pdgfra_experiment_a-preprocessing-filtering.qmd:** code from Mark, generation of filtered datasets  
**01_pdgfra_experiment_b-dimred-clustering.qmd:** code from Mark, dimreduction (without annotation)  
**01_pdgfra_experiment_c-annotation.qmd:** code from Mark, annotation of clusters  
**01_pdgfra_experiment_d-coverage_tracks.qmd:** preliminary code for Pdgfra coverage tracks  
**01_pdgfra_experiment_d-tuft.qmd:** code from Mark, annotation of Tuft and enteroendocrine clusters  
**01_pdgfra_experiment_e-bring-annotations-together.qm:** code from Mark, bringing all annotations together  
**01_pdgfra_experiment_e-unimodal_vs_multimodal.qmd:** preliminary code for comparison unimodal vs multimodal  
**01_pdgfra_experiment_f-junction-counts.qmd:** code from Mark, counting over exon junctions  
**01_pdgfra_experiment_f-unimodal_multimodal_paper_figures.qmd:** comparison of unimodal vs multimodal data. Supp Figure 10I-J  
**01_pdgfra_experiment_g-paper_figures_unmod_multimodal.qmd:** figures for paper. All panels Figure 5, Supp Figure 10A-H  

# Figure panels and corresponding .qmd file

**All panels Figure 5, Supp Figure 10A-H:** 01_pdgfra_experiment_g-paper_figures_unmod_multimodal.qmd  
**Supp Figure 10I-J:** 01_pdgfra_experiment_f-unimodal_multimodal_paper_figures.qmd
