process QIIME_BARPLOT {
    label 'process_low'
    
    container 'quay.io/qiime2/core:2023.2'
    
    input:
    path(counts)
    path(taxonomy)
    path(metadata)
    val(taxa_max)

    output:
    path('*exported_QIIME_barplot/*')     , emit: barplot_export
    path("level-*.csv")                   , emit: barplot_composition
    path("level-${taxa_max}.csv")             , emit: krona_tsv
    path('allsamples_compbarplot.qzv')    , emit: qzv
    path "versions.yml"                   , emit: versions

    script:   
    """
    mkdir tmpdir
    export TMPDIR=\$PWD/tmpdir
    
    qiime taxa barplot --i-table $counts --i-taxonomy $taxonomy --m-metadata-file $metadata --o-visualization allsamples_compbarplot.qzv
    qiime tools export --input-path allsamples_compbarplot.qzv --output-path allsamples_exported_QIIME_barplot
    
    array=( \$( seq 1 $taxa_max ) )

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
