/* 
Align COVID reads from each matched sample to the ref genome
and generate a sorted and indexed bam file 
*/ 

process COVID_ALIGNMENT {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mulled-v2-48af4db78cb1b6b5f01b94358a9657ea8abe95e4:d92e8cf3f6050b6f7782a5c2e986e220ef5f6bd4'

    input:
    tuple val(meta), path(covid_reads)
    path(covid_ref_ch)

    output:
    tuple val(meta), path("*.bam")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    bwa index ${covid_ref_ch} 
    bwa mem -t 32 ${covid_ref_ch} ${covid_reads[0]} ${covid_reads[1]} | samtools view -Sb - | samtools sort -o "${prefix}.bam" -
    """
}
