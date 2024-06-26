/* 
With the identified COVID samples, generate a read channel with paired end reads
from these COVID samples.
*/ 

process COVID_SAMPLE_EXTRACTION {
    label 'process_low'

    input:
    path(covid_samples_file)
    channel reads_ch

    output:
    channel covid_reads_ch

    exec:
    // Read COVID sample IDs into a set
    def covid_samples = covid_samples_file.text.split("\n").collect { it.trim() }.toSet()

    // Filter the reads channel
    covid_reads_ch = reads_ch.filter { sample_id, reads -> 
        covid_samples.contains(sample_id)
    }
}