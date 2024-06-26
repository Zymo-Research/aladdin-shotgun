/* 
Align COVID reads from each matched sample to the ref genome
and generate a sorted and indexed bam file 
*/ 

process COVID_ALIGNMENT {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/bwa:0.7.3a--hed695b0_5' 

    input:
    tuple val(sample_id), path(covid_reads)
    path(reference)

    output:
    tuple val(meta), path("*.bam")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    bwa index ${reference} 
    bwa mem -t 32 ${reference} ${covid_reads[0]} ${covid_reads[1]} | samtools view -Sb - | samtools sort -o "${prefix}.bam" -
    """
}
