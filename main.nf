#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def helpMessage() {
    log.info """

    Usage:

    Mandatory arguments:
      --readPaths     Tab-seperated file with sample names and path to the fastq files. (Used if --reads not provided.)
      -profile            Configuration profile to use. slurm? / docker / test ???

    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Examples below

include {rename_snp_ids; filltags; filter_vcf_bgz; retain_biallelic; retain_biallelic_snp_indels; ldprune; keep_samples; keep_samples_force; annotate_dbsnp; recompress} from './modules/filters.nf'


workflow filter_snps_indels_only {

    def vcf_files = Channel.fromPath(params.readPaths)
    //check_params()
    dbsnpfile = params.dbsnpfile
    dbsnpfile_index = params.dbsnpfile_index

    if( params.keepSamples == ".") {
        retain_biallelic_snp_indels(vcf_files)    
    } else {
        keep_samples(vcf_files, params.keepSamples)
        retain_biallelic_snp_indels(keep_samples.output)
    }

    rename_snp_ids(retain_biallelic_snp_indels.output[0])
    //annotate_dbsnp(retain_biallelic_snp_indels.output[0], dbsnpfile, dbsnpfile_index)
    //filltags(annotate_dbsnp.output[0])
}

workflow all {

    def vcf_files = Channel.fromPath(params.readPaths)
    //check_params()

    if( params.keepSamples == ".") {
        retain_biallelic(vcf_files)    
    } else {
        keep_samples(vcf_files, params.keepSamples)
        retain_biallelic(keep_samples.output)
    }
    
    filter_vcf_bgz(retain_biallelic.output, params.maf) 
    filltags(filter_vcf_bgz.out.vcf)
    rename_snp_ids(filltags.out.vcf)
    //ldprune(filter_vcf.output)
}

workflow geuvadis {
    def vcf_files = Channel.fromPath(params.readPaths)

    keep_samples_force(vcf_files, params.keepSamples)
    filltags(keep_samples_force.output)
    rename_snp_ids(filltags.out.vcf)
    // Filter MAF at the end, more convenient
    filter_vcf_bgz(rename_snp_ids.out.vcf, params.maf) 
}

workflow annotate_only {

    def vcf_files = Channel.fromPath(params.filteredVcfs)
    dbsnpfile = params.dbsnpfile
    dbsnpfile_index = params.dbsnpfile_index
    
    // recompress(vcf_files)
    // annotate_dbsnp(recompress.output.vcf_bgz, recompress.output.index_vcf, dbsnpfile, dbsnpfile_index)
    annotate_dbsnp(vcf_files, dbsnpfile, dbsnpfile_index)
    filltags(annotate_dbsnp.output)
}

workflow filltags_only {

    def vcf_files = Channel.fromPath(params.filteredVcfs)
    filltags(vcf_files)
}

workflow rename_snpids_only {

    def vcf_file = Channel.fromPath(params.filteredVcfs)
    rename_snp_ids(vcf_file)
}

def check_params() {

    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
    exit 0

    // additional validation here
}