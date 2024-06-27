/* 
With the identified COVID samples, generate a read channel with paired end reads
from these COVID samples.
*/ 

process COVID_SAMPLE_EXTRACTION {
    cache false
    tag "$meta.id"
    label 'process_low'

    input:
    path covid_samples_file
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), emit: covid_reads_ch

    script:
    """
    covidSamples=($(cat ${covid_samples_file}))

    samples=()

    for sample_id in "${covidSamples[@]}"; do
        if [[ "${sample_id}" == "${meta.id}" ]]; then
            samples+=($(tuple(meta, reads)))
        fi
    done

    emit(covid_reads_ch, samples)
    """
}