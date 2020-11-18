# SZS Pipeline

This repository contains code to perform the analyses in the paper ["The SZS is a robust and efficient method to identify regulated splicing events in droplet-based RNA Sequencing" (Olivieri, Dehghannasiri, and Salzman 2020)](https://www.biorxiv.org/content/10.1101/2020.11.10.377572v1.abstract). 

This pipeline takes the output from [SICILIAN](https://github.com/salzmanlab/SICILIAN) and returns the SZS for each gene and cell, as well as analyses of differential alternative splicing.

![Pipeline](pipeline.png)


## Installation

To run the pipeline, first set up the conda environment from the environment.yml file:

`conda env create --name szs_env --file=environments.yml`

and activate it:

` source activate szs_env`

You will need to place the following files in the "data" directory:
* `HLCA4_P2_10x_with_postprocessing_lung.pq`
* `HLCA4_P3_10x_with_postprocessing_lung.pq`

And the following files in the util_files directory:
* `GRCh38_latest_genomic.gtf`
* `ucscGenePfam.txt`

## Running the pipeline

Then run `snakemake -p` in the main directory (I run `snakemake -p --profile slurm` to run on sherlock). You can run `snakemake -np` first to see what jobs will be run.

## References
Olivieri, Dehghannasiri, and Salzman. "The SZS is an efficient statistical method to identify regulated splicing events in droplet-based RNA sequencing." bioRxiv. (2020) [https://www.biorxiv.org/content/10.1101/2020.11.10.377572v1.abstract](https://www.biorxiv.org/content/10.1101/2020.11.10.377572v1.abstract).

Dehghannasiri, Olivieri, and Salzman. "Specific splice junction detection in single cells with SICILIAN," bioRxiv. (2020) [https://www.biorxiv.org/content/10.1101/2020.04.14.041905v1](https://www.biorxiv.org/content/10.1101/2020.04.14.041905v1).
