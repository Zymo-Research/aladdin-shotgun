/* 
Merge sample feature table qza -> raw merged qza, not used in other process
Filter sample by no. reads -> filtered merged qza, used in barplot
Merge sample taxonomy and import to qiime -> taxonomy qza, used in barplot
Collapse to intended taxonomy level -> filtered collapsed qza, used in diversity analysis
Export to tsv -> used metadata filtering and delivered to customer
*/ 

process QIIME_DATAMERGE {
    label 'process_low'

    container 'quay.io/qiime2/core:2023.2'    

    input:
    path(abs_qza)
    path(taxonomy)

    output:
    path('merged_filtered_counts.qza')           , emit: filtered_counts_qza
    path('merged_taxonomy.qza')                  , optional: true, emit: taxonomy_qza
    path('merged_filtered_counts_collapsed.qza') , optional: true, emit: filtered_counts_collapsed_qza
    path('merged_filtered_counts_collapsed.tsv') , optional: true, emit: filtered_counts_collapsed_tsv
    path('versions.yml')                         , emit: versions
    path('*.qza')
    path('merged_filtered_counts.tsv')           , optional: true

    script:
    """
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
        --output-path merged_filtered_counts_out

    biom summarize-table -i merged_filtered_counts_out/feature-table.biom > biom_table_summary.txt
    SAMPLES=\$(grep 'Num samples' biom_table_summary.txt | sed 's/Num samples: //')
    FEATURES=\$(grep 'Num observations' biom_table_summary.txt | sed 's/Num observations: //')

    if [ "\$SAMPLES" -gt 0 ] && [ "\$FEATURES" -gt 0 ]
    then
        biom convert -i merged_filtered_counts_out/feature-table.biom -o merged_filtered_counts.tsv --to-tsv

        qiime_taxmerge.py $taxonomy
        qiime tools import \
            --input-path merged_taxonomy.tsv \
            --type 'FeatureData[Taxonomy]' \
            --input-format TSVTaxonomyFormat \
            --output-path merged_taxonomy.qza

        qiime taxa collapse \
            --i-table merged_filtered_counts.qza \
            --i-taxonomy merged_taxonomy.qza \
            --p-level ${params.taxonomy_collapse_level} \
            --o-collapsed-table merged_filtered_counts_collapsed.qza

        qiime tools export \
            --input-path merged_filtered_counts_collapsed.qza \
            --output-path merged_filtered_counts_collapsed_out
        biom convert -i merged_filtered_counts_collapsed_out/feature-table.biom -o merged_filtered_counts_collapsed.tsv --to-tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
        biom: \$(biom --version | sed 's/biom, version //')
    END_VERSIONS 
    """
}