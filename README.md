# Study overview
This folder contains R scripts that were used to analyze single-cell RNA-sequencing data from peripheral blood mononuclear cells (PBMCs) from sepsis patients (N=16) as well as healthy controls (N=6).

# Overview of scripts
'210817_1_Data_processing.R': This script is used to filter, normalize and annotate cells.

'210824_2_Cell_proportions.R': This script calculates and visualizes differences in cell proportions. 

'210824_3_DE_analyses.R': This script calculates differentially expressed genes between conditions per cell type, and visualizes the DEGs in volcano plots.

'210824_4_Pathway_enrichment_analyses.R': This script performs pathway enrichment analyses based on KEGG and GO terms and visualizes the results in dotplots.

'220525_5_Gene_enrichment_analyses.R': This script calculates the enrichment of a set of genes in the various monocyte clusters identified.

'220525_6_Identifying_regulatory_transcription_factors.R': This script identifies the regulatory transcription factors based on a list of DE genes. 
