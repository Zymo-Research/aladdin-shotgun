process QIIME2_EXPORT_ABSOLUTE {
    label 'process_low'

    container "quay.io/qiime2/core:2023.2"
    
    input:
    path(table)
    path(taxonomy_qza)
    path(taxonomy_tsv)
    val(tax_agglom_min)
    val(tax_agglom_max)

    output:
    path("feature-table.tsv")        , emit: tsv
    path("feature-table.biom")       , emit: biom
    path("table-[2-7].qza")          , emit: collapse_qza
    path("abs-abund-table-*.tsv")    , emit: abundtable
    path "versions.yml"              , emit: versions

    script:
    """
    mkdir tmpdir
    export TMPDIR=\$PWD/tmpdir

    #produce raw count table in biom format "table/feature-table.biom"
    qiime tools export --input-path ${table}  \
        --output-path table
    cp table/feature-table.biom .

    #produce raw count table "table/feature-table.tsv"
    biom convert -i table/feature-table.biom \
        -o feature-table.tsv  \
        --to-tsv

    biom summarize-table -i feature-table.biom > biom_table_summary.txt
    SAMPLES=\$(grep 'Num samples' biom_table_summary.txt | sed 's/Num samples: //')
    FEATURES=\$(grep 'Num observations' biom_table_summary.txt | sed 's/Num observations: //')

    levels=( "Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species" )

    if [ "\$SAMPLES" -gt 0 ] && [ "\$FEATURES" -gt 0 ]
    then

        ##on several taxa level
        array=(\$(seq ${tax_agglom_min} 1 ${tax_agglom_max}))
        for i in \${array[@]}
        do
            #collapse taxa
            qiime taxa collapse \
                --i-table ${table} \
                --i-taxonomy ${taxonomy_qza} \
                --p-level \$i \
                --o-collapsed-table table-\$i.qza
            #export to biom
            qiime tools export --input-path table-\$i.qza \
                --output-path table-\$i
            #convert to tab separated text file
            biom convert \
                -i table-\$i/feature-table.biom \
                -o abs-abund-table-\${levels[\${i}-1]}.tsv --to-tsv
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
