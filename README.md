# Pre-processing workflow for 1000KGP data

This respository is intended to pre-process 1000KGP vcf data to be workable with our workflows

This will:
* Filter SNPs with MAF < 0.001
* Retain only biallelic SNPs
* Perform LD-prunning
* Generate a single VCF file with all chromosomes

## How to run

First install Nextlfor:

```
curl -s "https://get.sdkman.io" | bash
sdk install java 17.0.6-tem
wget -qO- https://get.nextflow.io | bash
```
 
Perform a dry-run to check
`nextflow run main.nf -preview `