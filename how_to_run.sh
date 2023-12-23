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

## Run filtering and keep samples for AnswerALS
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