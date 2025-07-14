//
// Run profiling
//


/// Subworkflow
include { QIIME2_EXPORT                                } from '../../subworkflows/local/qiime2_export'
include { QIIME2_DIVERSITY                             } from '../../subworkflows/local/qiime2_diversity'
include { QIIME2_ANCOMBC                               } from '../../subworkflows/local/qiime2_ancombc'

// Modules
include { QIIME_IMPORT                                  } from '../../modules/nf-core/qiime/import/main'
include { QIIME2_FILTERSAMPLES                          } from '../../modules/local/qiime2_filtersamples'
include { QIIME2_PREPTAX                                } from '../../modules/local/qiime2_preptax'
include { QIIME_BARPLOT                                 } from '../../modules/nf-core/qiime/barplot/main'
include { KRONA_REPORT                                  } from '../../subworkflows/local/krona_report'
include { GROUP_COMPOSITION                             } from '../../modules/local/group_composition'
include { HEATMAP_INPUT                                 } from '../../modules/local/heatmap_input'

workflow DIVERSITY {
    take:
    qiime_profiles // [ [ meta ], absolute counts file ]
    qiime_taxonomy // 
    groups // group_metadata.csv
    tax_agglom_min
    tax_agglom_max
    ancombc_fdr_cutoff

    main:
    ch_versions          = Channel.empty()
    ch_multiqc_files     = Channel.empty()
    ch_output_file_paths = Channel.empty()
    ch_excel             = Channel.empty()
    ch_warning_message   = Channel.empty()

    QIIME_IMPORT ( qiime_profiles )
    ch_versions = ch_versions.mix( QIIME_IMPORT.out.versions )

    QIIME2_FILTERSAMPLES( QIIME_IMPORT.out.absabun_qza.collect() )
    ch_versions = ch_versions.mix( QIIME2_FILTERSAMPLES.out.versions )
    
    QIIME2_PREPTAX( qiime_taxonomy.collect() )

    QIIME_BARPLOT( QIIME2_FILTERSAMPLES.out.filtered_counts_qza, QIIME2_PREPTAX.out.taxonomy_qza, groups, tax_agglom_max )
    ch_versions = ch_versions.mix( QIIME_BARPLOT.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( QIIME_BARPLOT.out.barplot_composition.collect() )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME_BARPLOT.out.qzv.map{ "${params.outdir}/qiime2/composition_barplot/" + it.getName() }
        )

    if ( params.run_krona ){
        KRONA_REPORT( QIIME_BARPLOT.out.krona_tsv )
        ch_versions = ch_versions.mix( KRONA_REPORT.out.versions )
        ch_output_file_paths = ch_output_file_paths.mix(KRONA_REPORT.out.html.map{ "${params.outdir}/krona/" + it.getName() } )
    }
 
    HEATMAP_INPUT( QIIME_BARPLOT.out.barplot_composition.collect(), groups, params.top_taxa )
    ch_multiqc_files = ch_multiqc_files.mix( HEATMAP_INPUT.out.taxo_heatmap.collect())

    QIIME2_EXPORT( QIIME2_FILTERSAMPLES.out.filtered_counts_qza, QIIME2_PREPTAX.out.taxonomy_qza, QIIME2_PREPTAX.out.taxonomy_tsv, tax_agglom_min, tax_agglom_max )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME2_EXPORT.out.abs_taxa_levels.flatten().map{ "${params.outdir}/qiime2/export/" + it.getName() }
        )

    if (params.group_of_interest != 'NONE') {
    	ch_excel = Channel.fromPath(params.groupinterest[params.group_of_interest].excel_path, checkIfExists: true)
        GROUP_COMPOSITION( ch_excel, QIIME_BARPLOT.out.barplot_composition.collect() )
        ch_multiqc_files = ch_multiqc_files.mix( GROUP_COMPOSITION.out.compcsv.collect() )
        ch_output_file_paths = ch_output_file_paths.mix(
            GROUP_COMPOSITION.out.search_results.map { "${params.outdir}/groups_of_interest/" + it.getName() }
            )
    }

    QIIME2_DIVERSITY( QIIME2_FILTERSAMPLES.out.filtered_counts_qza, QIIME2_FILTERSAMPLES.out.filtered_counts_tsv, groups ) 
    ch_multiqc_files = ch_multiqc_files.mix( QIIME2_DIVERSITY.out.mqc )
    ch_output_file_paths = ch_output_file_paths.mix( QIIME2_DIVERSITY.out.output_paths )
    ch_versions = ch_versions.mix( QIIME2_DIVERSITY.out.versions )
    
    QIIME2_ANCOMBC( QIIME2_DIVERSITY.out.filtered_metadata, QIIME2_EXPORT.out.collapse_qza, QIIME2_PREPTAX.out.taxonomy_qza, tax_agglom_min, tax_agglom_max, ancombc_fdr_cutoff )
    ch_output_file_paths = ch_output_file_paths.mix(
        QIIME2_ANCOMBC.out.ch_output_files.flatten().map{ "${params.outdir}/qiime2/ancombc/visualizations/qzv/" + it.getName() }
        )
    ch_multiqc_files = ch_multiqc_files.mix( QIIME2_ANCOMBC.out.mqc.collect() )

    emit:
    versions     = ch_versions          // channel: [ versions.yml ]
    mqc          = ch_multiqc_files
    output_paths = ch_output_file_paths
    tables       = QIIME2_FILTERSAMPLES.out.filtered_counts_qza
    taxonomy     = QIIME2_PREPTAX.out.taxonomy_qza
    metadata     = QIIME2_DIVERSITY.out.ref_metadata
    warning      = ch_warning_message
}
