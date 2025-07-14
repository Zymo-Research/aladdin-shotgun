

include { QIIME_METADATAFILTER                          } from '../../modules/nf-core/qiime/metadatafilter/main'
include { QIIME_FILTER_SINGLETON_SAMPLE                 } from '../../modules/nf-core/qiime/filter_singleton_sample/main'
include { QIIME_ALPHARAREFACTION                        } from '../../modules/nf-core/qiime/alpha_rarefaction/main'
include { QIIME_DIVERSITYCORE                           } from '../../modules/nf-core/qiime/diversitycore/main'
include { QIIME_ALPHADIVERSITY                          } from '../../modules/nf-core/qiime/alphadiversity/main'
include { QIIME_BETAGROUPCOMPARE                        } from '../../modules/nf-core/qiime/beta_groupcompare/main'
include { QIIME_PLOT_MULTIQC                            } from '../../modules/nf-core/qiime/plot_multiqc/main'

workflow QIIME2_DIVERSITY {
    take:
    lvlmax_qza
    lvlmax_tsv
    groups // group_metadata.csv

    main:
    ch_versions          = Channel.empty()
    ch_multiqc_files     = Channel.empty()
    ch_output_file_paths = Channel.empty()

    QIIME_METADATAFILTER( groups, lvlmax_tsv )

    QIIME_FILTER_SINGLETON_SAMPLE( lvlmax_qza, QIIME_METADATAFILTER.out.filtered_metadata )
    ch_versions = ch_versions.mix( QIIME_FILTER_SINGLETON_SAMPLE.out.versions )
 
    QIIME_ALPHARAREFACTION( QIIME_METADATAFILTER.out.filtered_metadata, QIIME_FILTER_SINGLETON_SAMPLE.out.abs_qza, QIIME_METADATAFILTER.out.min_total )
    ch_versions = ch_versions.mix( QIIME_ALPHARAREFACTION.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_ALPHARAREFACTION.out.qzv.map{ "${params.outdir}/qiime2/diversity/alpha_rarefaction/" + it.getName() }
        )

    QIIME_DIVERSITYCORE( QIIME_FILTER_SINGLETON_SAMPLE.out.abs_qza, QIIME_METADATAFILTER.out.min_total, QIIME_METADATAFILTER.out.filtered_metadata )
    ch_versions = ch_versions.mix( QIIME_DIVERSITYCORE.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_DIVERSITYCORE.out.qzv.flatten().map{ "${params.outdir}/qiime2/diversity/diversity_core/" + it.getName() }
        )

    QIIME_ALPHADIVERSITY( QIIME_DIVERSITYCORE.out.vector.flatten(), QIIME_METADATAFILTER.out.filtered_metadata.collect() )
    ch_versions = ch_versions.mix( QIIME_ALPHADIVERSITY.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_ALPHADIVERSITY.out.qzv.flatten().map{ "${params.outdir}/qiime2/diversity/alpha_diversity/" + it.getName() }
        )

    QIIME_BETAGROUPCOMPARE ( QIIME_DIVERSITYCORE.out.distance.flatten(), QIIME_METADATAFILTER.out.filtered_metadata.collect() )
    ch_versions = ch_versions.mix( QIIME_BETAGROUPCOMPARE.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_BETAGROUPCOMPARE.out.qzv.flatten().map{ "${params.outdir}/qiime2/diversity/beta_diversity/" + it.getName() }
        )

    QIIME_PLOT_MULTIQC(
        QIIME_METADATAFILTER.out.filtered_metadata,
        QIIME_DIVERSITYCORE.out.pcoa.ifEmpty([]),
        QIIME_ALPHADIVERSITY.out.alphadiversity_tsv.collect().ifEmpty([]),
        QIIME_ALPHARAREFACTION.out.rarefaction_csv.collect().ifEmpty([]),
        false,
        params.skip_alpha_rarefaction )
    ch_multiqc_files = ch_multiqc_files.mix( QIIME_PLOT_MULTIQC.out.mqc_plot.collect() )

    emit:
    versions          = ch_versions          // channel: [ versions.yml ]
    mqc               = ch_multiqc_files
    output_paths      = ch_output_file_paths
    filtered_abs_qza  = QIIME_FILTER_SINGLETON_SAMPLE.out.abs_qza
    filtered_metadata = QIIME_METADATAFILTER.out.filtered_metadata
    ref_metadata      = QIIME_METADATAFILTER.out.ref_comp_metadata
}

