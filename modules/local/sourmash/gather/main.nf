process SOURMASH_GATHER {
    tag "$meta.id"
    container 'zxl124/sourmash_branchwater:0.9.3'
    if (params.ignore_failed_samples) {
        errorStrategy { task.attempt>1 ? 'ignore' : 'retry' }
    }

    input:
    tuple val(meta), path(sketch), path(sketch_log)
    path "sourmash_database/*"

    output:
    tuple val(meta), path('*with-lineages.csv'), emit: gather
    path "versions.yml", emit: versions 
    path "*.log"

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    DB=`find -L "sourmash_database" -name "*k${params.sourmash_kmersize}.zip"`
    LINEAGE=`find -L "sourmash_database" -name "*.csv"`

    sourmash scripts fastgather \\
        $sketch \\
        \$DB \\
        --ksize ${params.sourmash_kmersize} \\
        --threshold-bp ${params.sourmash_threshold_bp} \\
        -o ${prefix}_fastgather.csv 2> ${prefix}_fastgather.log
    sourmash gather \\
        $sketch \\
        \$DB \\
        --dna \\
        --ksize ${params.sourmash_kmersize} \\
        --threshold-bp ${params.sourmash_threshold_bp} \\
        -o ${prefix}_sourmashgather.csv \\
        --picklist ${prefix}_fastgather.csv:match_name:ident 1> ${prefix}_sourmashgather.log
    sourmash tax annotate -g ${prefix}_sourmashgather.csv -t \$LINEAGE 2> ${prefix}_sourmashannotate.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS

    """
}

