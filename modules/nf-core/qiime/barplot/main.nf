process QIIME_BARPLOT {
    label 'process_low'
    
    container 'quay.io/qiime2/core:2023.9'
    
    input:
    path(counts)
    path(taxonomy)
    path(metadata)

    output:
    path('*exported_QIIME_barplot/*')     , emit: barplot_export
    path('*exported_QIIME_barplot/*.csv') , emit: barplot_composition
    path('allsamples_compbarplot.qzv')    , emit: qzv
    path "versions.yml"                   , emit: versions

    script:   
    """
    qiime taxa barplot --i-table $counts --i-taxonomy $taxonomy --m-metadata-file $metadata --o-visualization allsamples_compbarplot.qzv
    qiime tools export --input-path allsamples_compbarplot.qzv --output-path allsamples_exported_QIIME_barplot
    
    array=( \$( seq 1 ${params.taxonomy_collapse_level}) )

    for i in \${array[@]}
    do
        qiimebarplot_tomultiqc.py -d allsamples_exported_QIIME_barplot/level-\$i.csv -l \$i
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
    END_VERSIONS
    """
}
