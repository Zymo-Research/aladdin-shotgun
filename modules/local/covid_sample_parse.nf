/* 
Parse the qiime filtered counts output and identify COVID positive samples 
as those having minimum read count above a certain threshold
*/ 

process COVID_SAMPLE_PARSE {
    label 'process_low'

    input:
    path(taxonomy_sample_counts)
    val(threshold)

    output:
    path('covid_samples.txt')

    script:
    """
    id_covid_samples.py ${taxonomy_sample_counts} ${threshold} >> covid_samples.txt
    """

}