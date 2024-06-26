//
// Run covid variant identification from mixed sample dataset
//

include { COVID_SAMPLE_PARSE                                                } from '../../modules/local/covid_sample_parse'
include { COVID_SAMPLE_EXTRACTION                                           } from '../../modules/local/covid_sample_extraction'
include { COVID_READ_EXTRACTION                                             } from '../../modules/local/covid_read_extraction'
include { COVID_ALIGNMENT                                                   } from '../../modules/local/covid_alignment'
include { DEMIX as COVID_VARID_DEMIX; AGGREGATE as COVID_VARID_AGGREGATE    } from '../../modules/local/covid_varID'


workflow COVID_VAR_ANNOTATION{
    take:
    filtered_counts_ch
    val covid_threshold_ch
    reads_ch
    covid_kraken_ch
    covid_ref_ch

    main:
    COVID_SAMPLE_PARSE(filtered_counts_ch,covid_threshold_ch)
    COVID_SAMPLE_EXTRACTION(COVID_SAMPLE_PARSE.out,reads_ch)
    COVID_READ_EXTRACTION(COVID_SAMPLE_EXTRACTION.out, covid_kraken_ch)
    COVID_ALIGNMENT(COVID_READ_EXTRACTION.out, covid_ref_ch)
    DEMIX(COVID_ALIGNMENT.out, covid_ref_ch)
    AGGREGATE(DEMIX.out.demixed)

    AGGREGATE.out.map{ "${params.outdir}/freyja/" + it.getName() }

}