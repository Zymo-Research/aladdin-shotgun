/*
 * Export filtered tables from QIIME2
 */
include { FIND_MAX_AVAILABLE_TAX } from '../../modules/local/find_max_available_tax'
include { QIIME2_EXPORT_ABSOLUTE } from '../../modules/local/qiime2_export_absolute'

workflow QIIME2_EXPORT {
    take:
    ch_asv
    taxonomy_qza
    taxonomy_tsv
    tax_agglom_min
    tax_agglom_max

    main:
    //Find max available taxonomy level, check against specified tax levels
    FIND_MAX_AVAILABLE_TAX ( taxonomy_tsv )
    FIND_MAX_AVAILABLE_TAX.out.max_tax.subscribe { it ->
        if (it.toInteger() < tax_agglom_max) { log.warn "Max available taxonomy is $it, but requested tax_agglom_max=$tax_agglom_max, switching to $it" }
        if (it.toInteger() < tax_agglom_min) { log.warn "Max available taxonomy is $it, but requested tax_agglom_min=$tax_agglom_min, switching to $it" }
    }

    tax_min = FIND_MAX_AVAILABLE_TAX.out.max_tax.toInteger().map{ [it, tax_agglom_min].min() }
    tax_max = FIND_MAX_AVAILABLE_TAX.out.max_tax.toInteger().map{ [it, tax_agglom_max].min() }

    //export_filtered_dada_output (optional)
    QIIME2_EXPORT_ABSOLUTE ( ch_asv, taxonomy_qza, taxonomy_tsv, tax_min, tax_max )

    emit:
    abs_tsv             = QIIME2_EXPORT_ABSOLUTE.out.tsv
    abs_taxa_levels     = QIIME2_EXPORT_ABSOLUTE.out.abundtable
    collapse_qza        = QIIME2_EXPORT_ABSOLUTE.out.collapse_qza
    collapse_tsv        = QIIME2_EXPORT_ABSOLUTE.out.abundtable
}
