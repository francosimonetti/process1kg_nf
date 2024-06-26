nextflow.enable.dsl=2

process filter_vcf_bgz {

    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 2
    memory "4 GB"

    input:
    path vcf_file
    val maf

    output:
    path "${vcf_file.getBaseName(2)}.MAF_${maf}.vcf.gz", emit: vcf
    path "${vcf_file.getBaseName(2)}.MAF_${maf}.vcf.gz.tbi", emit: vcf_index

    """
    bcftools view -q ${maf}:minor -Oz -o ${vcf_file.getBaseName(2)}.MAF_${maf}.vcf.gz $vcf_file
    bcftools index -t ${vcf_file.getBaseName(2)}.MAF_${maf}.vcf.gz
    """
}

process keep_samples {
    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 2
    memory "4 GB"

    input:
    path vcf_file
    path sample_file

    output:
    path "${vcf_file.getBaseName(2)}.samples.vcf.gz"

    """
    bcftools view -S ${sample_file} -Oz -o ${vcf_file.getBaseName(2)}.samples.vcf.gz $vcf_file
    bcftools index -t "${vcf_file.getBaseName(2)}.samples.vcf.gz"
    """
}

process keep_samples_force {
    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 2
    memory "4 GB"

    input:
    path vcf_file
    path sample_file

    output:
    path "${vcf_file.getBaseName(2)}.samples.vcf.gz"

    """
    bcftools view -S ${sample_file} --force-samples -Oz -o ${vcf_file.getBaseName(2)}.samples.vcf.gz $vcf_file
    bcftools index -t "${vcf_file.getBaseName(2)}.samples.vcf.gz"
    """
}

process retain_biallelic {

    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 2
    memory "4 GB"

    input:
    path vcf_file

    output:
    path "${vcf_file.getBaseName(2)}.biallelic.vcf.gz"
    path "*.vcf.gz.tbi"

    """
    bcftools view -m2 -M2 -v snps -Oz -o ${vcf_file.getBaseName(2)}.biallelic.vcf.gz $vcf_file
    bcftools index -t "${vcf_file.getBaseName(2)}.biallelic.vcf.gz"
    """
}

process retain_biallelic_snp_indels {

    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 2
    memory "4 GB"

    input:
    path vcf_file

    output:
    path "${vcf_file.getBaseName(2)}.biallelic_snps_indels.vcf.gz"
    path "*.vcf.gz.tbi"

    """
    bcftools view -m2 -M2 -v snps,indels -Oz -o ${vcf_file.getBaseName(2)}.biallelic_snps_indels.vcf.gz $vcf_file
    bcftools index -t "${vcf_file.getBaseName(2)}.biallelic_snps_indels.vcf.gz"
    """
}

process ldprune {

    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 4
    memory "8 GB"

    input:
    path vcf_file    

    output:
    path "${vcf_file.getBaseName(2)}.ldpruned.vcf.gz"

   """
   R=0.5
   W="200kb"
   bcftools +prune -m \$R -w \$W $vcf_file -Oz -o ${vcf_file.getBaseName(2)}.ldpruned.vcf.gz
   """
}

// TODO: doesn't work, unfinished
// process concat_files {
//     publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

//     cpus 4
//     memory "8 GB"

//     input:
//     path vcf_file_list    

//     output:
//     path "${params.study_name}.complete.vcf.gz"

//     """
//     bcftools concat -f $vcf_file_list -Oz -o ${params.study_name}.complete.vcf.gz
//     """
// }

process annotate_dbsnp {
    //publishDir "${params.outdir}/annotated_vcfs/", mode: 'copy'

    cpus 4
    memory "8 GB"

    input:
    path vcf_file
    path dbsnp_file
    path dbsnp_file_index

    output:
    path "${vcf_file.getBaseName(2)}.dbSNP.vcf.gz"
    path "*.vcf.gz.tbi"

    """
    bcftools index -t ${vcf_file}
    bcftools annotate -a ${dbsnp_file} -c ID -Oz -o ${vcf_file.getBaseName(2)}.dbSNP.vcf.gz ${vcf_file}
    """
}

process filltags {

    publishDir "${params.outdir}/annotated_vcfs/", mode: 'copy'

    cpus 4
    memory "8 GB"

    input:
    path vcf_file
    
    output:
    path "${vcf_file.getBaseName(2)}.updated.vcf.gz", emit: vcf
    path "*.vcf.gz.tbi"

    """
    bcftools +fill-tags -Oz -o ${vcf_file.getBaseName(2)}.updated.vcf.gz ${vcf_file} 
    bcftools index -t "${vcf_file.getBaseName(2)}.updated.vcf.gz"
    """
}


process rename_snp_ids {

    publishDir "${params.outdir}/annotated_vcfs", mode: 'copy'

    cpus 4
    memory "8 GB"

    input:
    path vcf_file

    output:
    path "${vcf_file.getBaseName(2)}.renamed_ids.vcf.gz", emit: vcf
    path "*.vcf.gz.tbi"


    """
    bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT\\_b38' -Oz -o ${vcf_file.getBaseName(2)}.renamed_ids.vcf.gz ${vcf_file}
    bcftools index -t "${vcf_file.getBaseName(2)}.renamed_ids.vcf.gz"
    """
}

process recompress {
    publishDir "${params.outdir}/filtered_vcfs/", mode: 'copy'

    cpus 4
    memory "8 GB"

    input:
    path vcf_file

    output:
    path "${vcf_file.getBaseName(2)}.bgz.vcf.gz", emit: vcf_bgz
    path "*.vcf.gz.tbi", emit: index_vcf

    """
    zcat ${vcf_file} | bcftools view -Oz -o "${vcf_file.getBaseName(2)}.bgz.vcf.gz" -
    bcftools index -t "${vcf_file.getBaseName(2)}.bgz.vcf.gz"
    """
}