/* 
Demix and deconvolute each of the mixed samples and 
identify the abundance of COVID variants in each sample
*/ 

process DEMIX {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/freyja:1.4.8--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam)
    path(covid_ref_ch)

    output:
    path("variants_files/*.variants.tsv")
    path("depth_files/*.depth")
    path("demix_files"), emit: demixed_dir

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    mkdir -p variants_files depth_files demix_files

    freyja variants ${bam} --variants variants_files/${prefix}.variants.tsv --depths depth_files/${prefix}.depth --ref ${covid_ref_ch}
    
    freyja update

    freyja demix variants_files/${prefix}.variants.tsv depth_files/${prefix}.depth --output demix_files/${prefix}.output --confirmedonly --depthcutoff 1

    """
}

process AGGREGATE {
    label 'process_low'
    container 'quay.io/biocontainers/freyja:1.5.1--pyhdfd78af_0'

    input:
    path demix_dirs

    output: 
    path("covid_variants.tsv")

    script:
    """
    mkdir -p merged_demix_files
    for dir in ${demix_dirs}; do
        cp \$dir/* merged_demix_files/
    done
    freyja aggregate merged_demix_files/ --output covid_variants.tsv
    """
}