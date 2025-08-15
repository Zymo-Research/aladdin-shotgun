/* 
Merge sample feature table qza -> raw merged qza, not used in other process
Filter sample by no. reads -> filtered merged qza, used in barplot, diversity, ancombc
Export to tsv -> used metadata filtering and delivered to customer
*/ 

process QIIME2_FILTERSAMPLES {
    label 'process_low'

    container 'quay.io/qiime2/core:2023.2'    

    input:
    path(abs_qza)

    output:
    path('merged_filtered_counts.qza')           , emit: filtered_counts_qza
    path('merged_filtered_counts.tsv')           , emit: filtered_counts_tsv
    path('versions.yml')                         , emit: versions

    script:
    """
    mkdir tmpdir
    export TMPDIR=\$PWD/tmpdir

    qiime feature-table merge \
        --i-tables $abs_qza \
        --o-merged-table merged_raw_counts.qza

    qiime feature-table filter-samples \
        --i-table merged_raw_counts.qza \
        --p-min-frequency ${params.lowread_filter} \
        --o-filtered-table merged_counts_filter-samples.qza
    
    qiime feature-table filter-features \
        --i-table merged_counts_filter-samples.qza \
        --p-min-frequency ${params.min_frequency} \
        --p-min-samples ${params.min_samples} \
        --o-filtered-table merged_filtered_counts.qza
    
    qiime tools export \
        --input-path merged_filtered_counts.qza \
        --output-path merged_filtered_counts_export
    biom convert -i merged_filtered_counts_export/feature-table.biom -o merged_filtered_counts.tsv --to-tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
        biom: \$(biom --version | sed 's/biom, version //')
    END_VERSIONS 
    """
}
