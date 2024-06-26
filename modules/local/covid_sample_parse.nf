/* 
Parse the qiime filtered counts output and identify COVID positive samples 
as those having minimum read count above a certain threshold
*/ 

process COVID_SAMPLE_PARSE {
    label 'process_low'

    input:
    path(filtered_counts_ch)
    val(covid_threshold_ch)

    output:
    path('covid_samples.txt'), emit:covid_samples_file

    script:
    """
    id_covid_samples.py ${filtered_counts_ch} ${covid_threshold_ch} > covid_samples.txt
    """

}