process QIIME_METADATAFILTER {
    label 'process_low'
   
    input:
    path( group_metadata ) 
    path( samples_filtered )

    output:
    path('filtered_metadata.tsv')                          , emit: filtered_metadata, optional: true

    script:   
    """
    metadata_pairwise.r $group_metadata $samples_filtered
    """
}