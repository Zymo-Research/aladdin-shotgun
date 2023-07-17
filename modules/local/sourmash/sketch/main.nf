process SOURMASH_SKETCH {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/sourmash:4.8.2--hdfd78af_0'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path('*.sig'), path('*.log'), emit: sketch
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash sketch dna -p k=51,scaled=1000,abund $input -o ${prefix}.sig --name ${prefix} 2> ${prefix}_sourmashsketch.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}
