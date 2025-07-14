process SOURMASH_QIIMEPREP {
    label 'process_low'
    tag "${meta.id}"

    conda "bioconda::biom-format=2.1.15"
    container 'quay.io/biocontainers/biom-format:2.1.15'

    input:
    tuple val(meta), path(gather), path(sketch_log)
    path host_lineage

    output:
    tuple val(meta), path ("*absabun_profile.biom"), optional: params.ignore_failed_samples, emit: biom
    path "*profile_taxonomy.txt", optional: params.ignore_failed_samples, emit: taxonomy
    path "*mqc.json", emit: mqc
    path "versions.yml", emit: versions
    path "*abun_parsed*", emit: abun_parsed
    
    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    host_lineage_param = host_lineage ? "--host_lineage $host_lineage" : ""
    """
    parse_sourmash_results_for_qiime.py $gather -n $prefix -l $sketch_log $host_lineage_param
    if [ -f "${prefix}_absabun_parsed_profile.txt" ]; then
        biom convert -i ${prefix}_absabun_parsed_profile.txt -o ${prefix}_absabun_profile.biom --table-type="OTU table" --to-json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biom: \$(biom --version | sed 's/biom, version //')
    END_VERSIONS
    """
}

