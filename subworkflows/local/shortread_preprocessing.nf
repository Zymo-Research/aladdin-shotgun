//
// Perform read trimming and merging
//


include { SHORTREAD_FASTP             } from './shortread_fastp'
include { SHORTREAD_ADAPTERREMOVAL    } from './shortread_adapterremoval'
include { FASTQC as FASTQC_PROCESSED  } from '../../modules/nf-core/fastqc/main'
include { FALCO as FALCO_PROCESSED    } from '../../modules/nf-core/falco/main'

workflow SHORTREAD_PREPROCESSING {
    take:
    reads //  [ [ meta ], [ reads ] ]
    adapterlist // file

    main:
    ch_versions        = Channel.empty()
    ch_multiqc_files   = Channel.empty()
    ch_warning_message = Channel.empty()

    if ( params.shortread_qc_tool == "fastp" ) {
        ch_processed_reads = SHORTREAD_FASTP ( reads, adapterlist ).reads
        ch_versions        = ch_versions.mix( SHORTREAD_FASTP.out.versions )
        ch_multiqc_files   = ch_multiqc_files.mix( SHORTREAD_FASTP.out.mqc )
    } else if ( params.shortread_qc_tool == "adapterremoval" ) {
        ch_processed_reads = SHORTREAD_ADAPTERREMOVAL ( reads, adapterlist ).reads
        ch_versions        = ch_versions.mix( SHORTREAD_ADAPTERREMOVAL.out.versions )
        ch_multiqc_files   = ch_multiqc_files.mix( SHORTREAD_ADAPTERREMOVAL.out.mqc )
    } else {
        ch_processed_reads = reads
    }

    //Filter empty files
    ch_processed_reads
        .branch{
            failed: it[0].single_end ? it[1].size() < 1.KB : it[1][0].size() < 1.KB || it[1][1].size() < 1.KB
            passed: it[0].single_end ? it[1].size() >= 1.KB : it[1][0].size() >= 1.KB && it[1][1].size() >= 1.KB
        }
        .set { ch_filtered_reads }
    ch_filtered_reads.passed.set { ch_processed_reads }
    ch_filtered_reads.failed
        .map { meta, reads -> [ meta.id ] }
        .collect()
        .map {
            msg = "The following samples had too small file size (<1kB) after trimming and read length filtering:\n${it.join("; ")}\nThis could happen when your sequencing quality is poor or your reads are not long enough for k-mer based taxonomy profiling."
            if (params.ignore_failed_samples) {
                log.warn "$msg"
                return msg
            } else {
                return error(msg)
            }
        }
        .set { ch_warning_message }

    if (params.preprocessing_qc_tool == 'fastqc') {
        FASTQC_PROCESSED ( ch_processed_reads )
        ch_versions = ch_versions.mix( FASTQC_PROCESSED.out.versions )
        if (params.posttrim_qc) {
            ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip )
        }
    } else if  (params.preprocessing_qc_tool == 'falco') {
        FALCO_PROCESSED ( ch_processed_reads )
        ch_versions = ch_versions.mix( FALCO_PROCESSED.out.versions )
        if (params.posttrim_qc) {
            ch_multiqc_files = ch_multiqc_files.mix( FALCO_PROCESSED.out.txt )
        }
    }

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
    warning  = ch_warning_message
}

