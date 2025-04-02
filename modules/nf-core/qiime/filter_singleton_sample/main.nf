// Remove samples from Qiime feature table that are not in the metadata file
// Covert to relative freq tsv for heatmap

process QIIME_FILTER_SINGLETON_SAMPLE {

    label 'process_low'

    container 'quay.io/qiime2/core:2023.9'

    input:
    path(qza)
    path(metadata)

    output:
    path('abs_feature_table_for_diversity.qza'), emit: abs_qza
    path('versions.yml')                       , emit: versions

    script:
    """
    qiime feature-table filter-samples \
        --i-table $qza \
        --m-metadata-file $metadata \
        --p-no-exclude-ids \
        --o-filtered-table abs_feature_table_for_diversity.qza
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
        biom: \$(biom --version | sed 's/biom, version //')
    END_VERSIONS 
    """
}