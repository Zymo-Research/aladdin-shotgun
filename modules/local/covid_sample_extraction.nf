/* 
With the identified COVID samples, generate a read channel with paired end reads
from these COVID samples.
*/ 

process COVID_SAMPLE_EXTRACTION {
    tag "$meta.id"
    label 'process_low'

    input:
    path covid_samples_file
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), emit: covid_reads_ch

    script:
    """
    # Read COVID sample IDs into a set
    covid_samples=\$(cat ${covid_samples_file} | tr '\\n' ' ')
    
    # Check if the sample_id is in covid_samples
    if [[ "\$covid_samples" =~ "${meta.id}" ]]; then
        echo "${meta.id},${reads[0]},${reads[1]}" > covid_sample_info.csv
    
    fi
    """
}