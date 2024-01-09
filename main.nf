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

include {rename_snp_ids; filltags; filter_vcf_bgz; retain_biallelic; ldprune; keep_samples; annotate_dbsnp; recompress} from './modules/filters.nf'


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