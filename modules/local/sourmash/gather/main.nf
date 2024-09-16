process SOURMASH_GATHER {
    tag "$meta.id"
    container 'quay.io/biocontainers/sourmash:4.8.2--hdfd78af_0'
            
    if (params.sourmash_db_ver == 'reps-rs214') {
        label = 'process_medium'
        maxForks = 20
    } else {
        label = 'process_high_memory'
        maxForks = 10
    }
    
    if (params.ignore_failed_samples) {
        errorStrategy = { task.exitStatus in [1,255] ? 'ignore' : 'retry' }
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
    DB=`find -L "sourmash_database" -name "*${params.sourmash_kmersize}.zip"`
    LINEAGE=`find -L "sourmash_database" -name "*.csv"`

    sourmash gather $sketch \$DB --dna --ksize ${params.sourmash_kmersize} --threshold-bp ${params.sourmash_threshold_bp} -o ${prefix}_sourmashgather.csv 2> ${prefix}_sourmashgather.log
    sourmash tax annotate -g ${prefix}_sourmashgather.csv -t \$LINEAGE 2> ${prefix}_sourmashannotate.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS

    """
}

