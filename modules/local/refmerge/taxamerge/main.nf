process REFMERGE_TAXAMERGE {
    label 'process_low'

    container "quay.io/qiime2/core:2023.9"

    input:
    path(ch_user_table), stageAs: 'user_table.qza'
    path(ch_tax), stageAs: 'user_taxonomy.qza'
    path(ch_ref_table), stageAs: 'ref_table.qza'
    path(ch_ref_tax), stageAs: 'ref_taxonomy.qza'

    output:
    path("refmerged_filtered-table.qza") , emit: merged
    path("refmerge_filtered-table.tsv")  , emit: tsv
    path("min_total.txt")                , emit: min_total

    script:
    """
    # collapse user and ref taxa table to genus level, respectively
    qiime taxa collapse \
        --i-table user_table.qza \
        --i-taxonomy user_taxonomy.qza \
        --p-level ${params.taxonomy_collapse_level} \
        --o-collapsed-table user_table-${params.taxonomy_collapse_level}.qza 
    
    qiime taxa collapse \
        --i-table ref_table.qza \
        --i-taxonomy ref_taxonomy.qza \
        --p-level ${params.taxonomy_collapse_level} \
        --o-collapsed-table ref_table-${params.taxonomy_collapse_level}.qza

    # merge user and ref tables
    qiime feature-table merge \
        --i-tables user_table-${params.taxonomy_collapse_level}.qza ref_table-${params.taxonomy_collapse_level}.qza \
        --o-merged-table refmerged-table.qza

    qiime feature-table filter-features \
        --i-table refmerged-table.qza \
        --p-min-frequency ${params.min_frequency} \
        --p-min-samples ${params.min_samples} \
        --o-filtered-table refmerged_filtered-table.qza
    
    # produce raw count table in biom format "table/feature-table.biom"
    qiime tools export \
        --input-path refmerged_filtered-table.qza  \
        --output-path table
    
    # produce raw count table
    biom convert -i table/feature-table.biom \
        -o table/feature-table.tsv  \
        --to-tsv
    cp table/feature-table.tsv refmerge_filtered-table.tsv

    # get min total counts
    get_min_total_counts.py -c refmerge_filtered-table.tsv -o min_total.txt
    """
}
