process QIIME2_ANCOMBC_INITIAL {
    tag "${table.baseName} - taxonomic level: ${taxlevel}"
    label 'process_low'
    label 'single_cpu'

    container "quay.io/qiime2/core:2023.2"

    input:
    tuple path(metadata), path(table), path(taxonomy)

    output:
    tuple path("*ancombc.qza"), val(taxlevel)  , emit: ancom, optional: true
    path("*ancombc.qza")
    tuple env(failfilter), val(taxlevel)       , emit: failcheck
    path "versions.yml"                        , emit: versions

    script:
    taxlevel = table.toString() - "table-" - ".qza"
    """
    mkdir tmpdir
    export TMPDIR=\$PWD/tmpdir

    mkdir ancom

    # Extract summarised table and output a file with the number of taxa
    qiime tools export --input-path ${table} --output-path exported/
    biom convert -i exported/feature-table.biom -o ${table.baseName}-level-${taxlevel}.feature-table.tsv --to-tsv

    if [ \$(grep -v '^#' -c ${table.baseName}-level-${taxlevel}.feature-table.tsv) -lt 3 ]; then
        failfilter="true"
    else
        qiime composition ancombc \
            --i-table ${table} \
            --m-metadata-file ${metadata} \
            --p-formula group \
            --o-differentials lvl${taxlevel}-ancombc.qza
        failfilter="false"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
