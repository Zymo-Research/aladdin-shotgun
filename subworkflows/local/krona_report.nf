/*
Taxonomic classification with QIIME2
*/

include { KRONA_PARSE                  } from '../../modules/local/krona_parse'
include { KRONA_RUN                    } from '../../modules/local/krona_run'

workflow KRONA_REPORT {
    take:
    ch_reltsv

    main:
    KRONA_PARSE ( ch_reltsv )

    KRONA_RUN ( KRONA_PARSE.out.krona_input.collect() )
    ch_versions = KRONA_RUN.out.versions 

    emit:
    versions  = ch_versions
    html = KRONA_RUN.out.html
}
