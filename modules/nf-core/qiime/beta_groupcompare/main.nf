process QIIME_BETAGROUPCOMPARE {
    tag "${distance.baseName}"
    label 'process_low'
    container "quay.io/qiime2/core:2023.9"

    when:
    !params.skip_betadiversity

    input:
    path(distance)
    path(metadata)

    output:
    path("beta_diversity/*"), emit: beta
    path("beta_diversity/*.qzv"), emit: qzv
    path("beta_diversity/*.tsv"), optional:true, emit: tsv
    path "versions.yml", emit: versions

    script:
    """ 
    qiime diversity beta-group-significance \
        --i-distance-matrix ${distance} \
        --m-metadata-file ${metadata} \
        --m-metadata-column "group" \
        --o-visualization ${distance.baseName}-group.qzv \
        --p-pairwise
    qiime tools export --input-path ${distance.baseName}-group.qzv \
        --output-path beta_diversity/${distance.baseName}-group
    #rename the output file name
    mv beta_diversity/${distance.baseName}-group/raw_data.tsv beta_diversity/${distance.baseName}-group.tsv
    mv ${distance.baseName}-group.qzv beta_diversity/    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
    END_VERSIONS
    """
}
