/* 
Demix and deconvolute each of the mixed samples and 
identify the abundance of COVID variants in each sample
*/ 

process freyja {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/freyja:1.5.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam)
    path(reference)

    output:
    path("variants_files/*.variants.tsv")
    path("depth_files/*.depth")
    path("demix_files/*.output")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    mkdir -p variants_files depth_files demix_files

    freyja variants ${bam} --variants variants_files/${prefix}.variants.tsv --depths depth_files/${prefix}.depth --ref ${reference}
    
    freyja update

    freyja demix variants_files/${prefix}.variants.tsv depth_files/${prefix}.depth --output demix_files/${prefix}.output --confirmedonly --depthcutoff 1

    """
}

process aggregate {
    label 'process_low'
    container 'quay.io/biocontainers/freyja:1.5.1--pyhdfd78af_0'

    input:
    path("demix_files/")

    output: 
    path("covid_variants.tsv"), emit: aggregated_variants

    script:
    """
    freyja aggregate demix_files/ --output covid_variants.tsv
    """
}