//
// Run covid variant identification from mixed sample dataset
//

include { COVID_SAMPLE_PARSE } from './modules/covid_sample_parse.nf'
include { COVID_SAMPLE_EXTRACTION } from './modules/covid_sample_extraction.nf'
include { COVID_READ_EXTRACTION } from './modules/covid_read_extraction.nf'
include { COVID_ALIGNMENT } from './modules/covid_alignment.nf'
include { COVID_VARID } from './modules/covid_varID.nf'

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
    COVID_VARID(COVID_ALIGNMENT.out, ref_genome)

    COVID_VARID.out.aggregated_variants.map{ "${params.outdir}/freyja/" + it.getName() }

}