## Metagenome analysis

This directory contains all scripts and code used to analyze the metagenomes and make raw figures for Band et al. (2026). 

### BSS_Bacterial-Archaeal_v9.Rmd

R notebook containing code for statistical analysis of alpha diversity, beta diversity, differential abundance, and clustering related to comparisons between day 0 and other timepoints after dosing. Raw figures related to these anaylses are also generated for manuscript. Figures were finalized (fonts, size, organization of panels) in Adobe Illustrator. 

### AMR_merge_table.py 

This python script works to merge all AMRs from WGSA2 outputs and their abundance (CPM) into one table in tsv format.

To run the script, use the following code if your working directory is the WGSA2 output folder (i.e. PipelineResultXXXXXXXXX).

```
python3 AMR_merge_table.py ./asmbMetaSpades
```

If there are multiple WGSA2 runs, copy all into one directory and include the path that directory.
