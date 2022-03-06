# SCLC_project

This repository specifically downloads and proccesses the bioproject [PRJNA608275](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP250470&o=acc_s%3Aa) .\
The different steps involve:
* Download and conversion to fastq files
* Align to Hg38
* Differential splicing analysis using Majiq
* Differential expression analysis using Salmon and limma


### Run as:
```
snakemake -s workflow/snakemake.smk --profile slurm -p --restart-times 3 --use-conda --dryrun
```
