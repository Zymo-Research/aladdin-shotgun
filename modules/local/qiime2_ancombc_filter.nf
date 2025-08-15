process QIIME2_ANCOMBC_FILTER {
    tag "taxonomic level: ${taxlevel}"
    label 'process_low'
    label 'single_cpu'
    label 'process_long'
    label 'error_ignore'

    container "quay.io/qiime2/core:2023.2"

    input:
    tuple path(metadata), path(qza), val(taxlevel)
    val(taxa_max)
    val(fdr_cutoff)

    output:
    path "to_multiqc/*lvl${taxa_max}_visualization_export/*ancombc-barplot.html", emit: to_mqc, optional: true
    path "refgroup.txt", emit: ref_group
    path "qzv/*visualization.qzv", emit: ancombc_vis
    path "versions.yml" , emit: versions

    script:
    """
    mkdir tmpdir
    export TMPDIR=\$PWD/tmpdir

    qiime tools export --input-path $qza --output-path ancombc_export/

    mkdir to_multiqc/
    mkdir qzv/

    parse_ancom_output.py -l "./ancombc_export/lfc_slice.csv" -s "./ancombc_export/se_slice.csv" -q "./ancombc_export/q_val_slice.csv"  -w "./ancombc_export/w_slice.csv" -p "./ancombc_export/p_val_slice.csv" -f $metadata -c $fdr_cutoff

    REFGROUP=\$( cat refgroup.txt )

    for dir in \$(find ./to_multiqc/ -maxdepth 1 -mindepth 1 -type d) 
    do 
        qiime tools import \
            --input-path \$dir \
            --type 'FeatureData[DifferentialAbundance]' \
            --output-path \${dir}_lvl${taxlevel}_ancombc.qza

        qiime composition da-barplot --i-data \${dir}_lvl${taxlevel}_ancombc.qza --o-visualization \${dir}_lvl${taxlevel}_visualization 

        qiime tools export --input-path \${dir}_lvl${taxlevel}_visualization.qzv --output-path \${dir}_lvl${taxlevel}_visualization_export/

        grouponly=\$( echo \${dir} | sed -e 's/.\\/to_multiqc\\///' )

        [ ! -f \${dir}_lvl${taxlevel}_visualization_export/\${grouponly}-ancombc-barplot.html ] && > \${dir}_lvl${taxlevel}_visualization_export/\${grouponly}-ancombc-barplot.html

        mv \${dir}_lvl${taxlevel}_visualization.qzv ./qzv/\${REFGROUP}_vs_\${grouponly}_lvl${taxlevel}_ancombc_visualization.qzv

    done 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
