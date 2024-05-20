# Spatial profiling of classic Hodgkin lymphoma (cHL)

#### Code supporting manuscript:
#### *Multiplexed spatial profiling of Hodgkin Reed-Sternberg cell neighborhoods in Classic Hodgkin lymphoma*

## Multiplexed IF analysis

- Install `snakemake`, `R 4.0.2`, and all libraries in `statistics/source_all.R`

- For reformatting of HALO output (csv files) to `RDA` files, please reference: https://github.com/mskcc/cHL-pre-processing
- This is a fork of `HaloX` repo branch `ver/halo_v3.5.1`

- `scripts/misc` contains R scripts for cell reassignments, UMAP clustering, and HRS aggregate analysis. `reassignment_rules.md` contains cell reassignment hierarchy. 

- `statistics` contains R scripts for all statistical analyses.

## Bulk transcriptomic analysis

- `NanoString/hodgkin_rna_analysis.Rmd` contains all code for analysis of RNA sequencing data and generation of paper figures. 

## Raw data

- Please reference: https://zenodo.org/records/10659311
- Contains raw multiplexed IF data and files referenced in `NanoString/hodgkin_rna_analysis.Rmd`.
