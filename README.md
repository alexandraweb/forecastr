## Purpose
Estimate STR variant size from SNP haplotype data using STR callset from WGS of a reference population

Reference STR VCFs can be generated by calling STR variants from the WGS datasets using:  
ExpansionHunter - https://github.com/Illumina/ExpansionHunter  
GangSTR - https://github.com/gymreklab/GangSTR   


## Requirements
Snakemake (v7.3.1) - https://snakemake.readthedocs.io/en/stable/getting_started/installation.html  
Singularity (v3.9.8) - https://docs.sylabs.io/guides/3.0/user-guide/installation.html  

## Install
`git clone https://github.com/alexandraweb/forecastr.git`

## Test installation
`snakemake --configfile config/config_test.yml --use-singularity --use-conda --cores all`

## Quick Start
Fill out the configfile at `config.yml` with the approiate file paths and run:
`snakemake --configfile config/config.yml --use-singularity --use-conda --cores all`

## Output
VCF file containing the predicted STR variants  
`{outdir}/forecastr.vcf.gz`


