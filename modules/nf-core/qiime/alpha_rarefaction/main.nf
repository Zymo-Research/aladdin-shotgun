process QIIME_ALPHARAREFACTION {
    label 'process_low'

    container "quay.io/qiime2/core:2023.9"

    when:
    !params.skip_alpha_rarefaction
    
    input:
    path(metadata)
    path(table)
    val(min_total)

    output:
    path("alpha-rarefaction/*")     , emit: rarefaction
    path("alpha-rarefaction/*.csv") , optional:true, emit: rarefaction_csv
    path("alpha-rarefaction.qzv")   , emit: qzv
    path "versions.yml"             , emit: versions

    script:
    def maxdepth = min_total.toInteger() < 250000 ? min_total.toInteger() : 250000 
    """
    qiime diversity alpha-rarefaction  \
        --i-table $table  \
        --p-max-depth $maxdepth  \
        --m-metadata-file $metadata  \
        --p-steps 50  \
        --p-iterations 5 \
        --o-visualization alpha-rarefaction.qzv
    qiime tools export --input-path alpha-rarefaction.qzv  \
        --output-path alpha-rarefaction

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
    END_VERSIONS
    """
}
