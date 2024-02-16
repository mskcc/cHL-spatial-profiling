# Spatial profiling of classic Hodgkin lymphoma (cHL)

#### Code supporting manuscript:
#### "Multiplexed spatial profiling of Hodgkin Reed-Sternberg cell neighborhoods in Classic Hodgkin lymphoma unveils distinct immune escape mechanisms"

## Multiplexed IF analysis

- For Halo->RDA used: https://github.com/mskcc/cHL-pre-processing, which is a fork of `HaloX` repo branch `ver/halo_v3.5.1`

*scripts/misc* contains R scripts for cell reassignments, UMAP clustering, HRS aggregate analysis. 

*statistics* contains R scripts for all statistical analyses. 

Refer to *reassignment.md* for Hodgkin cell reassignment hierarchy. 


## Transcriptome

Contains R script for analysis of NanoString data and RNA sequencing data with some plotting for figures. 

Contains *rawdata* files used in R script. 

Install snakemake, install R 4.0.2 and all libraries in source_all.R

To run analysis, run
