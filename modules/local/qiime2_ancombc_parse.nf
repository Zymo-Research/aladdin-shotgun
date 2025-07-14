process QIIME2_ANCOMBC_PARSE {
    label 'process_low'

    input:
    path(html)
    path(reference_group)

    output:
    path"*_mqc.html", emit: ancombc_plot

    script:
    """
    replace_vis.sh "$html"
    """
}
