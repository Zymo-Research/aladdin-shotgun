/* 
Merge and import sample taxonomy to qiime -> taxonomy qza, used in barplot
*/ 

process QIIME2_PREPTAX {
    label 'process_low'

    container 'quay.io/qiime2/core:2023.2'    

    input:
    path(taxonomy)

    output:
    path('merged_taxonomy.qza')                  , emit: taxonomy_qza
    path('merged_taxonomy.tsv')                  , optional: true, emit: taxonomy_tsv
    path('versions.yml')                         , emit: versions

    script:
    """
    mkdir tmpdir
    export TMPDIR=\$PWD/tmpdir

    qiime_taxmerge.py $taxonomy -o "merged_taxonomy.tsv"

    qiime tools import \
        --input-path merged_taxonomy.tsv \
        --type 'FeatureData[Taxonomy]' \
        --input-format TSVTaxonomyFormat \
        --output-path merged_taxonomy.qza


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | sed '2,2d' | sed 's/q2cli version //g')
    END_VERSIONS 
    """
}
