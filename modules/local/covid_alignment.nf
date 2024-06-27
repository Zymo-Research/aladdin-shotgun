/* 
Align COVID reads from each matched sample to the ref genome
and generate a sorted and indexed bam file 
*/ 

process COVID_ALIGNMENT {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/bwa:0.7.3a--hed695b0_5' 

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
