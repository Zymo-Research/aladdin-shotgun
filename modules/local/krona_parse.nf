process KRONA_PARSE {
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)

    input:
    path(rel_tsv)

    output:
    path "*.txt", emit: krona_input

    script:
    """
    parse_qiime_forkrona.py $rel_tsv    
    """
}
