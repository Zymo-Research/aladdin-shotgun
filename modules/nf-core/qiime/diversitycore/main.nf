process QIIME_DIVERSITYCORE {
    label 'process_low'
   
    container 'quay.io/qiime2/core:2023.9'
     
    input:
    path(qza)
    val(min_total)
    path(filtered_metadata)

    output:
    path('diversity_core/*_vector.qza')          , emit: vector
    path('diversity_core/*_distance_matrix.qza') , emit: distance
    path('diversity_core/*_pcoa_ordination.txt') , emit: pcoa
    path("diversity_core/*.qzv")                 , emit: qzv
    path "versions.yml"                          , emit: versions

    script:
    """
    qiime diversity core-metrics \
        --i-table $qza \
        --p-sampling-depth ${min_total.toInteger()} \
        --m-metadata-file $filtered_metadata \
        --output-dir diversity_core \
        --p-n-jobs ${task.cpus}
    
    for i in 'bray_curtis' 'jaccard'
    do
        qiime tools export \
            --input-path diversity_core/\${i}_pcoa_results.qza \
            --output-path pcoa_export
        mv pcoa_export/ordination.txt diversity_core/\${i}_pcoa_ordination.txt
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
    END_VERSIONS
    """
}
