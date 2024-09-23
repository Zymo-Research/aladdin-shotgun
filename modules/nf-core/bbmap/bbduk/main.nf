process BBMAP_BBDUK {
    tag "$meta.id"
    label 'process_dynamic_bbduk'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.06--h92535d8_1':
        'quay.io/biocontainers/bbmap:39.06--h92535d8_1' }"

    input:
    tuple val(meta), path(reads)
    path contaminants

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_1.fastq.gz out2=${prefix}_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : ''
    """
    maxmem=\$(echo \"$task.memory\"| cut -d ' ' -f 1 | cut -d '.' -f 1)
    bbduk.sh \\
        -Xmx\${maxmem}g \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        $contaminants_fa \\
        qin=${params.phred_encode} \\
        &> ${prefix}.bbduk.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
