process KRONA_PARSE {
    label 'process_low'

    input:
    path(rel_tsv)

    output:
    path "*.txt", emit: krona_input

    script:
    """
    parse_qiime_forkrona.py $rel_tsv    
    """
}
