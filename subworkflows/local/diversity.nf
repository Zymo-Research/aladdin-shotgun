//
// Run profiling
//

include { QIIME_IMPORT                                  } from '../../modules/nf-core/qiime/import/main'
include { QIIME_DATAMERGE                               } from '../../modules/nf-core/qiime/datamerge/main'
include { QIIME_METADATAFILTER                          } from '../../modules/nf-core/qiime/metadatafilter/main'
include { QIIME_FILTER_SINGLETON_SAMPLE                 } from '../../modules/nf-core/qiime/filter_singleton_sample/main'
include { QIIME_ALPHARAREFACTION                        } from '../../modules/nf-core/qiime/alpha_rarefaction/main'
include { QIIME_DIVERSITYCORE                           } from '../../modules/nf-core/qiime/diversitycore/main'
include { QIIME_BARPLOT                                 } from '../../modules/nf-core/qiime/barplot/main'
include { GROUP_COMPOSITION                             } from '../../modules/local/group_composition'
include { HEATMAP_INPUT                                 } from '../../modules/local/heatmap_input'
include { QIIME_ALPHADIVERSITY                          } from '../../modules/nf-core/qiime/alphadiversity/main'
include { QIIME_BETAGROUPCOMPARE                        } from '../../modules/nf-core/qiime/beta_groupcompare/main'
include { QIIME_PLOT_MULTIQC                            } from '../../modules/nf-core/qiime/plot_multiqc/main'

workflow DIVERSITY {
    take:
    qiime_profiles // [ [ meta ], absolute counts file ]
    qiime_taxonomy // 
    groups // group_metadata.csv

    main:
    ch_versions          = Channel.empty()
    ch_multiqc_files     = Channel.empty()
    ch_output_file_paths = Channel.empty()
    ch_excel             = Channel.empty()
    ch_warning_message   = Channel.empty()

    QIIME_IMPORT ( qiime_profiles )
    ch_versions = ch_versions.mix( QIIME_IMPORT.out.versions )

    QIIME_DATAMERGE( QIIME_IMPORT.out.absabun_qza.collect(), qiime_taxonomy.collect() )
    ch_versions = ch_versions.mix( QIIME_DATAMERGE.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_DATAMERGE.out.filtered_counts_collapsed_tsv.map{ "${params.outdir}/qiime_mergeddata/" + it.getName() }
        )
    QIIME_DATAMERGE.out.filtered_counts_collapsed_qza
        .ifEmpty('There were no samples or taxa left after filtering! Try lower filtering criteria or examine your data quality.')
        .filter( String )
        .set{ ch_warning_message }
 
    QIIME_BARPLOT( QIIME_DATAMERGE.out.filtered_counts_qza, QIIME_DATAMERGE.out.taxonomy_qza, groups )
    ch_versions = ch_versions.mix( QIIME_BARPLOT.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( QIIME_BARPLOT.out.barplot_composition.collect() )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_BARPLOT.out.qzv.map{ "${params.outdir}/qiime_composition_barplot/" + it.getName() }
        )
    
    if (params.group_of_interest != 'NONE') {
    	ch_excel = Channel.fromPath(params.groupinterest[params.group_of_interest].excel_path, checkIfExists: true)
        GROUP_COMPOSITION( ch_excel, QIIME_BARPLOT.out.barplot_composition.collect() )
        ch_multiqc_files = ch_multiqc_files.mix( GROUP_COMPOSITION.out.compcsv.collect() )
        ch_output_file_paths = ch_output_file_paths.mix(
            GROUP_COMPOSITION.out.search_results.map { "${params.outdir}/groups_of_interest/" + it.getName() }
            )
    }

    QIIME_METADATAFILTER( groups, QIIME_DATAMERGE.out.filtered_counts_collapsed_tsv )

    QIIME_FILTER_SINGLETON_SAMPLE( QIIME_DATAMERGE.out.filtered_counts_collapsed_qza, QIIME_METADATAFILTER.out.filtered_metadata )
    ch_versions = ch_versions.mix( QIIME_FILTER_SINGLETON_SAMPLE.out.versions )
        
    HEATMAP_INPUT( QIIME_BARPLOT.out.barplot_composition.collect(), QIIME_METADATAFILTER.out.filtered_metadata, params.top_taxa )
    ch_multiqc_files = ch_multiqc_files.mix( HEATMAP_INPUT.out.taxo_heatmap.collect()) 

    QIIME_ALPHARAREFACTION( QIIME_METADATAFILTER.out.filtered_metadata, QIIME_FILTER_SINGLETON_SAMPLE.out.abs_qza, QIIME_METADATAFILTER.out.min_total )
    ch_versions = ch_versions.mix( QIIME_ALPHARAREFACTION.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_ALPHARAREFACTION.out.qzv.map{ "${params.outdir}/qiime_alpha_rarefaction/" + it.getName() }
        )

    QIIME_DIVERSITYCORE( QIIME_FILTER_SINGLETON_SAMPLE.out.abs_qza, QIIME_METADATAFILTER.out.min_total, QIIME_METADATAFILTER.out.filtered_metadata )
    ch_versions = ch_versions.mix( QIIME_DIVERSITYCORE.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_DIVERSITYCORE.out.qzv.flatten().map{ "${params.outdir}/qiime_diversity/diversity_core/" + it.getName() }
        )

    QIIME_ALPHADIVERSITY( QIIME_DIVERSITYCORE.out.vector.flatten(), QIIME_METADATAFILTER.out.filtered_metadata.collect() )
    ch_versions = ch_versions.mix( QIIME_ALPHADIVERSITY.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_ALPHADIVERSITY.out.qzv.flatten().map{ "${params.outdir}/qiime_diversity/alpha_diversity/" + it.getName() }
        ) 

    QIIME_BETAGROUPCOMPARE ( QIIME_DIVERSITYCORE.out.distance.flatten(), QIIME_METADATAFILTER.out.filtered_metadata.collect() )
    ch_versions = ch_versions.mix( QIIME_BETAGROUPCOMPARE.out.versions )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_BETAGROUPCOMPARE.out.qzv.flatten().map{ "${params.outdir}/qiime_diversity/beta_diversity/" + it.getName() }
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
    versions     = ch_versions          // channel: [ versions.yml ]
    mqc          = ch_multiqc_files
    output_paths = ch_output_file_paths
    tables       = QIIME_DATAMERGE.out.filtered_counts_qza
    taxonomy     = QIIME_DATAMERGE.out.taxonomy_qza
    metadata     = QIIME_METADATAFILTER.out.ref_comp_metadata
    filtered_counts_ch = Channel.fromPath(QIIME_DATAMERGE.out.filtered_counts_collapsed_tsv)
    warning      = ch_warning_message
}
