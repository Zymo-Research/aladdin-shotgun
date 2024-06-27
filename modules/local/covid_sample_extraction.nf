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

    exec:
    def covidSamples = file(covid_samples_file).text.readLines().collect { it.trim() as String }.toSet()

    // Check if the sample_id is in covid_samples
    if (covidSamples.contains(meta.id)) {
        println "${meta.id}: ${reads}"
        emit(covid_reads_ch, tuple(meta, reads))
    }
}