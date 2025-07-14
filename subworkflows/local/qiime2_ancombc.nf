//
// Modules
//
include { QIIME2_ANCOMBC_INITIAL             } from '../../modules/local/qiime2_ancombc_initial'
include { QIIME2_ANCOMBC_FILTER              } from '../../modules/local/qiime2_ancombc_filter'
include { QIIME2_ANCOMBC_PARSE               } from '../../modules/local/qiime2_ancombc_parse'

workflow QIIME2_ANCOMBC {
    take:
    ch_metadata
    ch_asv
    ch_tax
    tax_agglom_min
    tax_agglom_max
    fdr_cutoff

    main:
    ch_metadata
        .combine( ch_asv.flatten() )
        .combine( ch_tax )
        .set{ ch_for_ancom_tax }
    QIIME2_ANCOMBC_INITIAL ( ch_for_ancom_tax )

    QIIME2_ANCOMBC_INITIAL.out.failcheck.subscribe{ failfilter, taxlevel -> 
        if (failfilter == "true") {
            log.warn( "WARNING: Summing your data at taxonomic level $taxlevel produced less than three rows (taxa), ANCOM can't proceed.")
        }        
    }

    QIIME2_ANCOMBC_FILTER ( ch_metadata.combine( QIIME2_ANCOMBC_INITIAL.out.ancom ), tax_agglom_max, fdr_cutoff )
    QIIME2_ANCOMBC_PARSE(QIIME2_ANCOMBC_FILTER.out.to_mqc, QIIME2_ANCOMBC_FILTER.out.ref_group )

    emit:
    mqc = QIIME2_ANCOMBC_PARSE.out.ancombc_plot
    ch_output_files = QIIME2_ANCOMBC_FILTER.out.ancombc_vis.collect()
}
