/* 
Extract COVID reads from matched samples (and remove any host & baterial reads) 
to reduce computation time during read alignment
*/ 

process COVID_READ_EXTRACTION {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_1'

    input:
    tuple val(meta), path(reads)
    path(kraken_db)

    output:
    tuple val(meta), path("*{1,2}.fq.gz")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    kraken2 --db ${kraken_db} --threads 32 --paired ${reads[0]} ${reads[1]} --classified-out ${prefix}#.fq.gz > /dev/null
    """
}