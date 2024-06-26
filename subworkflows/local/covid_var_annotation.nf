//
// Run covid variant identification from mixed sample dataset
//

include { COVID_SAMPLE_PARSE                              } from '../../modules/local/covid_sample_parse'
include { COVID_SAMPLE_EXTRACTION                         } from '../../modules/local/covid_sample_extraction'
include { COVID_READ_EXTRACTION                           } from '../../modules/local/covid_read_extraction'
include { COVID_ALIGNMENT                                 } from '../../modules/local/covid_alignment'
include { COVID_VARID.DEMIX, COVID_VARID.AGGREGATE        } from '../../modules/local/covid_varID'

workflow COVID_VAR_ANNOTATION{
    take:
    channel taxonomy_sample_counts
    channel covid_threshold
    channel reads
    channel kraken_db
    channel ref_genome

    main:
    COVID_SAMPLE_PARSE(taxonomy_sample_counts,covid_threshold)
    COVID_SAMPLE_EXTRACTION(COVID_SAMPLE_PARSE.out,reads)
    COVID_READ_EXTRACTION(COVID_SAMPLE_EXTRACTION.out, kraken_db)
    COVID_ALIGNMENT(COVID_READ_EXTRACTION.out, ref_genome)
    DEMIX(COVID_ALIGNMENT.out, ref_genome)
    AGGREGATE(DEMIX.out.demixed)

    AGGREGATE.out.map{ "${params.outdir}/freyja/" + it.getName() }

}