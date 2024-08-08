//
// Run covid variant identification from mixed sample dataset
//

include { COVID_SAMPLE_PARSE                                                } from '../../modules/local/covid_sample_parse'
include { COVID_READ_EXTRACTION                                             } from '../../modules/local/covid_read_extraction'
include { BWA as COVID_ALIGNMENT_BWA; SAMTOOLS as COVID_ALIGNMENT_SAMTOOLS  } from '../../modules/local/covid_alignment'
include { DEMIX as COVID_VARID_DEMIX; AGGREGATE as COVID_VARID_AGGREGATE    } from '../../modules/local/covid_varID'


workflow COVID_VAR_ANNOTATION {
    take:
    filtered_counts_ch
    covid_threshold_ch
    reads_ch
    covid_kraken_ch
    covid_ref_ch

    main:
    ch_multiqc_files     = Channel.empty()
    ch_output_file_paths = Channel.empty()

    COVID_SAMPLE_PARSE(filtered_counts_ch, covid_threshold_ch)

    // Read the contents of the covid_samples file into a channel
    covid_samples_ch = COVID_SAMPLE_PARSE.out.covid_samples_file
        .splitText()  
        .map { it.trim() }  
        .collect() 
        .map { samples -> [samples] }

    filtered_reads_ch = reads_ch
        .combine(covid_samples_ch) 
        .filter { meta, reads, covid_samples -> 
            covid_samples.contains(meta.id)
        }
        .map { meta, reads, covid_samples -> [meta, reads] }
    
    COVID_READ_EXTRACTION(filtered_reads_ch, covid_kraken_ch.collect())
    COVID_ALIGNMENT_BWA(COVID_READ_EXTRACTION.out, covid_ref_ch.collect())
    COVID_ALIGNMENT_SAMTOOLS(COVID_ALIGNMENT_BWA.out)
    COVID_VARID_DEMIX(COVID_ALIGNMENT_SAMTOOLS.out, covid_ref_ch.collect())
    COVID_VARID_AGGREGATE(COVID_VARID_DEMIX.out.demix_files.collect())

    ch_multiqc_files = ch_multiqc_files.mix(COVID_VARID_AGGREGATE.out.freyja_mqc_json) 
    ch_output_file_paths = ch_output_file_paths.mix(
        COVID_VARID_AGGREGATE.out.aggregated_tsv.map{ "${params.outdir}/freyja/" + it.getName() }
        )

    emit:
    mqc          = ch_multiqc_files
    output_paths = ch_output_file_paths
}