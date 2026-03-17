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
    tax_min
    tax_max

    main:
    //export_filtered_dada_output (optional)
    QIIME2_EXPORT_ABSOLUTE ( ch_asv, taxonomy_qza, taxonomy_tsv, tax_min, tax_max )

    emit:
    abs_tsv             = QIIME2_EXPORT_ABSOLUTE.out.tsv
    abs_taxa_levels     = QIIME2_EXPORT_ABSOLUTE.out.abundtable
    collapse_qza        = QIIME2_EXPORT_ABSOLUTE.out.collapse_qza
    collapse_tsv        = QIIME2_EXPORT_ABSOLUTE.out.abundtable
}
