# Pre-processing workflow for 1000KGP data

This respository is intended to pre-process 1000KGP vcf data to be workable with our workflows

This will:
* Retain only biallelic SNPs
* Select specific samples
* Filter SNPs with MAF < 0.001
* Perform LD-prunning
* Annotate SNPs with RSIDs from dbSNP or rename the SNP ids with consistent naming
* Update INFO/AC and INFO/AF VCF tags

## How to run

First install Nextlfor:

```
curl -s "https://get.sdkman.io" | bash
sdk install java 17.0.6-tem
wget -qO- https://get.nextflow.io | bash
```
 
Perform a dry-run to check
`nextflow run main.nf -preview `

To run a specific workflow you can do:

```
## Run biallelic and MAF filtering for all 1000G Samples
nextflow run main.nf -entry all \
    --readPaths "/my/datasets/1000G_high_coverage/latest/*.vcf.gz" \
    --keepSamples "." \
    --maf "0.001" \
    --outdir "/my/datasets/1000G_high_coverage/latest/processed"
```
