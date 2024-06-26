/* 
With the identified COVID samples, generate a read channel with paired end reads
from these COVID samples.
*/ 

process COVID_SAMPLE_EXTRACTION {
    tag "$meta.id"
    label 'process_low'

    input:
    path(covid_samples_file)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads)

    exec:
    // Read COVID sample IDs into a set
    def covid_samples = covid_samples_file.text.splitCsv(sep:'\\n', strip:true).collect { it[0] }.toSet()

    // Filter the reads tuple
    if (covid_samples.contains(meta.id)) {
        emit(tuple(meta, reads))
    }
}