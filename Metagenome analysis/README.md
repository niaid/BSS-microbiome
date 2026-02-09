## Metagenome analysis

This directory contains all scripts and code used to analyze the metagenomes and make raw figures for Band et al. (2026). 


### AMR_analysis_final.Rmd

R notebook containing code for data processing and statistical analysis of differential abundance of antibiotic resistance genes for day 0 and other timepoints after dosing. Raw figures related to these anaylses are also generated for manuscript. Figures were finalized (fonts, size, organization of panels) in Adobe Illustrator. 

### AMR_merge_table.py 

This python script works to merge all AMRs from WGSA2 outputs and their abundance (CPM) into one table in tsv format.

To run the script, use the following code if your working directory is the WGSA2 output folder (i.e. PipelineResultXXXXXXXXX).

```
python3 AMR_merge_table.py ./asmbMetaSpades
```

If there are multiple WGSA2 runs, copy all into one directory and include the path that directory.

### BSS_Bacterial-Archaeal_v9.Rmd

R notebook containing code for statistical analysis of alpha diversity, beta diversity, differential abundance, and clustering related to comparisons between day 0 and other timepoints after dosing. Raw figures related to these anaylses are also generated for manuscript. Figures were finalized (fonts, size, organization of panels) in Adobe Illustrator. 

### Eukaryote_analysis_final.Rmd

R notebook containing code for data processing and statistical analysis of differential abundance of eukaryotes for day 0 and other timepoints after dosing. Raw figures related to these anaylses are also generated for manuscript. Figures were finalized (fonts, size, organization of panels) in Adobe Illustrator. 

### Functional_pathways_analysis_final.Rmd

R notebook containing code for data processing and statistical analysis of differential abundance of functional pathways for day 0 and other timepoints after dosing. Raw figures related to these anaylses are also generated for manuscript. Figures were finalized (fonts, size, organization of panels) in Adobe Illustrator. 

### High_and_low_responders_final.Rmd

R notebook containing code for data processing and statistical analysis of differential abundance of genera for day 0 compared to day 2 after dosing in both high and low responders, as well as genera that are differentially abundant at day 8, 14, and 28 in responders. Raw figures related to these anaylses are also generated for manuscript. Figures were finalized (fonts, size, organization of panels) in Adobe Illustrator. 

### Mouse_metagenomes_pathways_batch-corrected.Rmd

R notebook containing code for data processing and statistical analysis of compositional analysis and differential abundance for functional pathways in mouse metagenomes. 

### BSS_mouse_metagenomes_tax_batch-corrected

R notebook containing code for statistical analysis of alpha diversity, beta diversity, and differential abundance related to comparisons between day 0 and other timepoints after dosing.

### Predictive taxa.Rmd

R notebook containing code for data processing and statistical analysis (ROC analysis) for identifying predicting genera and class of response at day 2. 

### Viral_analysis_final.Rmd

R notebook containing code for data processing and statistical analysis of alpha diversity of viral OTUs and differential abundance of viral OTUs classified by hosts between day 0 and other timepoints after dosing. Raw figures related to these anaylses are also generated for manuscript. Figures were finalized (fonts, size, organization of panels) in Adobe Illustrator.

### WGCNA_analysis_final.Rmd

R notebook containing code for WGCNA analysis of genera used with metabolomic network analysis. 

