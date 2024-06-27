//
// Run covid variant identification from mixed sample dataset
//

include { COVID_SAMPLE_PARSE                                                } from '../../modules/local/covid_sample_parse'
include { COVID_SAMPLE_EXTRACTION                                           } from '../../modules/local/covid_sample_extraction'
include { COVID_READ_EXTRACTION                                             } from '../../modules/local/covid_read_extraction'
include { BWA as COVID_ALIGNMENT_BWA; SAMTOOLS as COVID_ALIGNMENT_SAMTOOLS  } from '../../modules/local/covid_alignment'
include { DEMIX as COVID_VARID_DEMIX; AGGREGATE as COVID_VARID_AGGREGATE    } from '../../modules/local/covid_varID'


workflow COVID_VAR_ANNOTATION{
    take:
    filtered_counts_ch
    covid_threshold_ch
    reads_ch
    covid_kraken_ch
    covid_ref_ch

    main:
    COVID_SAMPLE_PARSE(filtered_counts_ch, covid_threshold_ch)

    covid_samples = COVID_SAMPLE_PARSE.out.covid_samples_file.collectFile(name: 'all_covid_samples.txt')

    filtered_reads_ch = reads_ch
        .join(covid_samples)
        .filter { meta, reads, covid_file -> 
            covid_file.text.split('\n').contains(meta.id)
        }
        .map { meta, reads, covid_file -> [meta, reads] }

    COVID_READ_EXTRACTION(filtered_reads_ch, covid_kraken_ch)
    COVID_ALIGNMENT_BWA(COVID_READ_EXTRACTION.out, covid_ref_ch)
    COVID_ALIGNMENT_SAMTOOLS(COVID_ALIGNMENT_BWA.out)
    COVID_VARID_DEMIX(COVID_ALIGNMENT_SAMTOOLS.out, covid_ref_ch)
    COVID_VARID_AGGREGATE(COVID_VARID_DEMIX.out.demixed_dir)

    COVID_VARID_AGGREGATE.out.map{ "${params.outdir}/freyja/" + it.getName() }

}