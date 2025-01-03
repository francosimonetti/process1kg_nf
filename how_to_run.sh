## Run filtering for all 1000G
nextflow run main.nf -entry all \
    --readPaths "/biodata/franco/datasets/1000G_high_coverage/latest/*.vcf.gz" \
    --keepSamples "." \
    --maf "0.001" \
    --outdir "/biodata/franco/datasets/1000G_high_coverage/latest/processed"

## Run annotation only for 1000GKP vcfs
nextflow run main.nf -entry annotate_only \
    --filteredVcfs "/biodata/franco/datasets/1000G_high_coverage/latest/processed/filtered_vcfs/1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.biallelic.MAF_0.001.vcf.gz" \
    --outdir "/biodata/franco/datasets/1000G_high_coverage/latest/processed"

## Run rename_ids for 1000G

nextflow run main.nf -entry rename_snpids_only \
    --filteredVcfs "/biodata/franco/datasets/1000G_high_coverage/latest/processed/1kGP_high_coverage.complete.allMAFs.vcf.gz" \
    --outdir "/biodata/franco/datasets/1000G_high_coverage/latest/processed/"

########

## Run filtering and keep samples for AnswerALS (VERSION V5)
nextflow run main.nf -entry all \
    --readPaths "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/AnswerALS-866-G-v1-release5_joint-vcf-vqsr-annotated.chr*.vcf.gz" \
    --keepSamples "/biodata/franco/datasets/answerALS/metadata_v6/genomic_samples_with_rnaseq_v6.txt" \
    --maf "0.01" \
    --outdir "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed_v6/"

## Run filtering and keep samples for AnswerALS (VERSION V6)
nextflow run main.nf -entry all \
    --readPaths "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/AnswerALS-866-G-v1-release5_joint-vcf-vqsr-annotated.chr*.vcf.gz" \
    --keepSamples "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/samples_with_rnaseq.txt" \
    --maf "0.01" \
    --outdir "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed/"

## Run annotation for AnswerALS
nextflow run main.nf -entry annotate_only \
    --filteredVcfs "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed/filtered_vcfs/AnswerALS-866-G-v1-release5_joint-vcf-vqsr-annotated.chr*.samples.biallelic.MAF_0.01.vcf.gz" \
    --outdir "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed/"

## Run filltags_only for AnswerALS
nextflow run main.nf -entry filltags_only \
    --filteredVcfs "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed/annotated_vcfs/AnswerALS-866-G-v1-release5_joint-vcf-vqsr-annotated.chr*.samples.biallelic.MAF_0.01.dbSNP.vcf.gz" \
    --outdir "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed/"

## Run rename_ids for AnswerALS
nextflow run main.nf -entry rename_snpids_only \
    --filteredVcfs "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed/annotated_vcfs/AnswerALS-866-G-v1-release5_joint-vcf-vqsr-annotated.chr*.samples.biallelic.MAF_0.01.dbSNP.updated.vcf.gz" \
    --outdir "/biodata/franco/datasets/answerALS/genomics/4_JointGenotyping/processed/annotated_vcfs"


## Run all for GTEx 2022
nextflow run main.nf -entry all \
    --readPaths "/biodata/franco/datasets/gtex_v8/genotypes/latest_20122023/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.chr*.vcf.gz" \
    --keepSamples "." \
    --maf "0.01" \
    --outdir "/biodata/franco/datasets/gtex_v8/genotypes/latest_20122023/processed"


## Run all for GTEx 2022 with INDELs!
nextflow run main.nf -entry filter_snps_indels_only \
    --readPaths "/biodata/franco/datasets/gtex_v8/genotypes/latest_20122023/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.chr*.vcf.gz" \
    --keepSamples "." \
    --outdir "/biodata/franco/datasets/gtex_v8/genotypes/latest_20122023/processed"

## Run annotation for GTEx 2022
nextflow run main.nf -entry annotate_only \
    --filteredVcfs "/biodata/franco/datasets/gtex_v8/genotypes/latest_20122023/processed/filtered_vcfs/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.chr*.biallelic.vcf.gz" \
    --outdir "/biodata/franco/datasets/gtex_v8/genotypes/latest_20122023/processed/"


## Run worflow to extract GEUVADIS samples from 1000G dataset
nextflow run main.nf -entry geuvadis \
    --readPaths "/biodata/franco/datasets/1000G_high_coverage/latest/processed/filtered_vcfs/1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.biallelic.vcf.gz" \
    --keepSamples "/biodata/franco/datasets/geuvadis/reprocess/fsimonetti-nfdata-aws/FULL/geuvadis.sample_ids.txt" \
    --maf "0.01" \
    --outdir "/biodata/franco/datasets/geuvadis/genotypes/latest_from1000G_high_coverage/"



## Run all for UCT agusting data, keep snps and indels but only biallelic
nextflow run main.nf -entry filter_snps_indels_only \
    --readPaths "/biodata/franco/agustin/ld_myh7/chr*.UCT.vcf.gz" \
    --keepSamples "." \
    --outdir "/biodata/franco/agustin/ld_myh7"

nextflow run main.nf -entry filter_snps_indels_only \
    --readPaths "/biodata/franco/agustin/genotipado_conjunto_UCT_malbran_hg38_2024_ok.vcf.gz" \
    --keepSamples "." \
    --outdir "/biodata/franco/agustin/"

## Run for MESA dataset
nextflow run main.nf -entry filter_snps_indels_only \
    --readPaths "/biodata/franco/datasets/mesa/downloads/genotype_merge/MESA_phs001416_TOPMed_WGS_freeze.9b.chr*.hg38.merged.vcf.gz" \
    --keepSamples "." \
    --outdir "/biodata/franco/datasets/mesa/downloads/genotype_merge/processed"
