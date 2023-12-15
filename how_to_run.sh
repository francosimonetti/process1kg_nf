## Run annotation only for 1000GKP vcfs
nextflow run main.nf -entry annotate_only \
    --filteredVcfs "/biodata/franco/datasets/1000G_high_coverage/latest/process1kg_nf/filtered_vcfs/1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.biallelic.MAF_0.001.vcf.gz" \
    --outdir "/biodata/franco/datasets/1000G_high_coverage/latest/processed"

## Run filtering for a single VCF 1000G
nextflow run main.nf -entry all \
    --readPaths /biodata/franco/datasets/1000G_high_coverage/latest/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

## Run filtering for all 1000G
nextflow run main.nf -entry all \
    --readPaths "/biodata/franco/datasets/1000G_high_coverage/latest/*.vcf.gz" \
    --keepSamples = "" \
    --maf = "0.001" \
    --outdir "/biodata/franco/datasets/1000G_high_coverage/latest/processed"

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