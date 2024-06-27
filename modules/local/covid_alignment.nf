/* 
Align COVID reads from each matched sample to the ref genome
and generate a sorted and indexed bam file 
*/ 

process BWA {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/bwa:0.7.17--ha92aebf_3'

    input:
    tuple val(meta), path(covid_reads)
    path(covid_ref_ch)

    output:
    tuple val(meta), path("*.sam")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    bwa index ${covid_ref_ch} 
    bwa mem -t ${task.cpus} ${covid_ref_ch} ${covid_reads[0]} ${covid_reads[1]} > ${prefix}.sam
    """
}


process SAMTOOLS {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-rsamtools:1.32.3--r351hfc679d8_0'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.bam")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    samtools view -Sb ${sam} | samtools sort -o "${prefix}.bam" -
    """
