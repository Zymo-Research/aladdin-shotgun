process MULTIQC {
    cache false
    container 'quay.io/biocontainers/multiqc:1.23--pyhdfd78af_0'

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path "multiqc_custom_plugins"
    val warnings

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def rtitle = params.run_name ? "--title \"Shotgun report for ${params.run_name}\"" : ''
    def rfilename = params.run_name ? "--filename " + params.run_name.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    def comment = warnings ? "--comment \"$warnings\"" : '' 
    """
    cd multiqc_custom_plugins
    python setup.py develop
    cd ..

    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        $rtitle $rfilename \\
        $comment \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_data
    touch multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
