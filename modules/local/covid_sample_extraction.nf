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
    # Read COVID sample IDs into an array
    mapfile -t covidSamples < "${covid_samples_file}"

    # Check if the sample_id is in covidSamples
    if printf '%s\n' "\${covidSamples[@]}" | grep -qx "${meta.id}"; then
        echo "MATCH" > status.txt
    else
        echo "NO_MATCH" > status.txt
    fi
    """
}