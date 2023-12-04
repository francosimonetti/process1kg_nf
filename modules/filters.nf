nextflow.enable.dsl=2

process filter_vcf {

    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 2
    memory "4 GB"

    input:
    path vcf_file

    output:
    path "${vcf_file.getBaseName(2)}.filtered.vcf.gz"


    """
    filter_vcf.py --in $vcf_file --out ${vcf_file.getBaseName(2)}.filtered.vcf.gz --maf 0.001
    """
}

// process retain_biallelic {
//     """
//     bcftools view -m2 -M2 -v snps -Oz -o $output $input
//     """
// }

// process ldprune {
    
//     R=0.5
//     W="200kb"
//     bcftools +prune -m ${R} -w ${W} $input -Oz -o $output.ldpruned_${R}_${W}.vcf.gz

// }

