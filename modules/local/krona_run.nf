process KRONA_RUN {
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container 'quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1'

    input:
    path krona_input

    output:
    path "krona_report.html", emit: html
    path "versions.yml", emit: versions


    script:
    """
    ktImportText -o krona_report.html $krona_input   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(ktImportXML 2>&1 | grep -o "KronaTools [0-9]\\.[0-9]\\.[0-9]" | sed -e 's/KronaTools /KronaTools: /' )
    END_VERSIONS

    """
}
