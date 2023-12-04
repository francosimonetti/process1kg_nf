#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def helpMessage() {
    log.info """

    Usage:

    Mandatory arguments:
      --readPathsFile     Tab-seperated file with sample names and path to the fastq files. (Used if --reads not provided.)
      -profile            Configuration profile to use. slurm? / docker / test ???

    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Examples below

include {filter_vcf} from './modules/filters.nf'

params.readPaths = "/data/franco/datasets/1000G_high_coverage/20220422_3202_phased_SNV_INDEL_SV/*.vcf.gz"

workflow {

    def vcf_files = Channel.fromPath(params.readPaths)
    filter_vcf(vcf_files)    

    // if (params.run_ge_quant){
    //     count_features(align_reads.out.bam_sorted_by_name)
    // }
}
