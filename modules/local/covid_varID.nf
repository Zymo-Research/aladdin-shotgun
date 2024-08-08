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
    path("variants_files/*.tsv")
    path("depth_files/*.depth")
    path("demix_files/*.output"), emit: demix_files

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    mkdir -p variants_files depth_files demix_files/
    freyja update
    freyja variants ${bam} --variants variants_files/${prefix} --depths depth_files/${prefix}.depth --ref ${covid_ref_ch}
    freyja demix variants_files/${prefix}.tsv depth_files/${prefix}.depth --output demix_files/${prefix}.output --confirmedonly --depthcutoff 1

    """
}

process AGGREGATE {
    label 'process_low'
    container 'quay.io/biocontainers/freyja:1.5.1--pyhdfd78af_0'

    input:
    path demix_files

    output: 
    path("*.tsv"), emit: aggregated_tsv
    path("*.json"), emit: freyja_mqc_json

    script:
    """
    mkdir -p demix_dir/
    mv ${demix_files} demix_dir/
    freyja aggregate demix_dir/ --output covid_variants.tsv
    display_freyja_mqc.py covid_variants.tsv
    """
}