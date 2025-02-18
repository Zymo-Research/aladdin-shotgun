//
// Run profiling
//

include { MALT_RUN                                      } from '../../modules/nf-core/malt/run/main'
include { MEGAN_RMA2INFO as MEGAN_RMA2INFO_TSV          } from '../../modules/nf-core/megan/rma2info/main'
include { KRAKEN2_KRAKEN2                               } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_STANDARD_REPORT                       } from '../../modules/local/kraken2_standard_report'
include { BRACKEN_BRACKEN                               } from '../../modules/nf-core/bracken/bracken/main'
include { CENTRIFUGE_CENTRIFUGE                         } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT                            } from '../../modules/nf-core/centrifuge/kreport/main'
include { KHMER_TRIM_LOW_ABUND                          } from '../../modules/local/khmer_trim_low_abund'
include { SOURMASH_SKETCH                               } from '../../modules/local/sourmash/sketch/main'
include { SOURMASH_GATHER                               } from '../../modules/local/sourmash/gather/main'
include { METAPHLAN4_METAPHLAN4                         } from '../../modules/nf-core/metaphlan4/metaphlan4/main'
include { METAPHLAN4_QIIMEPREP                          } from '../../modules/nf-core/metaphlan4/qiimeprep/main'
include { METAPHLAN4_UNMAPPED                           } from '../../modules/nf-core/metaphlan4/unmapped/main'
include { SOURMASH_QIIMEPREP                            } from '../../modules/local/sourmash/qiimeprep/main'
include { KAIJU_KAIJU                                   } from '../../modules/nf-core/kaiju/kaiju/main'
include { KAIJU_KAIJU2TABLE as KAIJU_KAIJU2TABLE_SINGLE } from '../../modules/nf-core/kaiju/kaiju2table/main'
include { DIAMOND_BLASTX                                } from '../../modules/nf-core/diamond/blastx/main'
include { MOTUS_PROFILE                                 } from '../../modules/nf-core/motus/profile/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ                } from '../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'

workflow PROFILING {
    take:
    reads // [ [ meta ], [ reads ] ]
    databases // [ [ meta ], path ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_raw_classifications  = Channel.empty()
    ch_raw_profiles         = Channel.empty()
    ch_qiime_profiles       = Channel.empty()
    ch_taxonomy             = Channel.empty()
    
    //COMBINE READS WITH POSSIBLE DATABASES
    // e.g. output [DUMP: reads_plus_db] [['id':'2612', 'run_accession':'combined', 'instrument_platform':'ILLUMINA', 'single_end':1], <reads_path>/2612.merged.fastq.gz, ['tool':'malt', 'db_name':'mal95', 'db_params':'"-id 90"'], <db_path>/malt90]
    ch_input_for_profiling = reads
            .map {
                meta, reads ->
                    [meta, reads]
            }//Not sure if this mapping is needed, will test later
            .combine(databases)

    /*
    PREPARE PROFILER INPUT CHANNELS & RUN PROFILING
    */

    // Each tool as a slightly different input structure and generally separate
    // input channels for reads vs databases. We restructure the channel tuple
    // for each tool and make liberal use of multiMap to keep reads/databases
    // channel element order in sync with each other

    if ( params.profiler == "malt" ) {

        if (!params.shortread_qc_mergepairs) log.warn "[nf-core/taxprofiler] MALT does not accept uncollapsed paired-reads. Pairs will be profiled as separate files."

        // MALT: We groupTuple to have all samples in one channel for MALT as database
        // loading takes a long time, so we only want to run it once per database
        ch_input_for_malt =  ch_input_for_profiling
            .map {
                meta, reads, db_meta, db ->

                    // Reset entire input meta for MALT to just database name,
                    // as we don't run run on a per-sample basis due to huge datbaases
                    // so all samples are in one run and so sample-specific metadata
                    // unnecessary. Set as database name to prevent `null` job ID and prefix.
                    def temp_meta = [ id: meta['db_name'] ]

                    // Extend database parameters to specify whether to save alignments or not
                    def new_db_meta = db_meta.clone()
                    def sam_format = params.malt_save_reads ? ' --alignments ./ -za false' : ""
                    new_db_meta['db_params'] = db_meta['db_params'] + sam_format

                    // Combine reduced sample metadata with updated database parameters metadata,
                    // make sure id is db_name for publishing purposes.
                    def new_meta = temp_meta + new_db_meta
                    new_meta['id'] = new_meta['db_name']

                    [ new_meta, reads, db ]

            }
            .groupTuple(by: [0,2])
            .multiMap {
                meta, reads, db ->
                    reads: [ meta, reads.flatten() ]
                    db: db
            }

        MALT_RUN ( ch_input_for_malt.reads, ch_input_for_malt.db )

        ch_maltrun_for_megan = MALT_RUN.out.rma6
                                .transpose()
                                .map{
                                    meta, rma ->
                                            // re-extract meta from file names, use filename without rma to
                                            // ensure we keep paired-end information in downstream filenames
                                            // when no pair-merging
                                            def meta_new = meta.clone()
                                            meta_new['db_name'] = meta.id
                                            meta_new['id'] = rma.baseName
                                        [ meta_new, rma ]
                                }

        MEGAN_RMA2INFO_TSV (ch_maltrun_for_megan, params.malt_generate_megansummary )
        ch_multiqc_files       = ch_multiqc_files.mix( MALT_RUN.out.log )
        ch_versions            = ch_versions.mix( MALT_RUN.out.versions.first(), MEGAN_RMA2INFO_TSV.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( ch_maltrun_for_megan )
        ch_raw_profiles        = ch_raw_profiles.mix( MEGAN_RMA2INFO_TSV.out.txt )

    }

    if ( params.profiler == "kraken2" ) {
        // Have to pick first element of db_params if using bracken,
        // as db sheet for bracken must have ; sep list to
        // distinguish between kraken and bracken parameters
        ch_input_for_kraken2 = ch_input_for_profiling
                                .map {
                                    meta, reads, db_meta, db ->
                                        def db_meta_new = db_meta.clone()

                                        // Only take second element if one exists
                                        def parsed_params = db_meta_new['db_params'].split(";")
                                        if ( parsed_params.size() == 2 ) {
                                            db_meta_new['db_params'] = parsed_params[0]
                                        } else if ( parsed_params.size() == 0 ) {
                                            db_meta_new['db_params'] = ""
                                        } else {
                                            db_meta_new['db_params'] = parsed_params[0]
                                        }

                                    [ meta, reads, db_meta_new, db ]
                                }
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

        KRAKEN2_KRAKEN2 ( ch_input_for_kraken2.reads, ch_input_for_kraken2.db, params.kraken2_save_reads, params.kraken2_save_readclassification )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKEN2_KRAKEN2.out.classified_reads_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix(
            KRAKEN2_KRAKEN2.out.report
                // Set the tool to be strictly 'kraken2' instead of potentially 'bracken' for downstream use.
                // Will remain distinct from 'pure' Kraken2 results due to distinct database names in file names.
                .map { meta, report -> [meta + [tool: 'kraken2'], report]}
        )

    }
    
    /*
    // Currently cannot support sequential runs of kraken2 and bracken, to be worked on later
    if ( params.run_kraken2 && params.run_bracken ) {
        // Remove files from 'pure' kraken2 runs, so only those aligned against Bracken & kraken2 database are used.
        def ch_kraken2_output = KRAKEN2_KRAKEN2.out.report
            .filter {
                meta, report ->
                    if ( meta['instrument_platform'] == 'OXFORD_NANOPORE' ) log.warn "[nf-core/taxprofiler] Bracken has not been evaluated for Nanopore data. Skipping Bracken for sample ${meta.id}."
                    meta['tool'] == 'bracken' && meta['instrument_platform'] != 'OXFORD_NANOPORE'
            }

        // If necessary, convert the eight column output to six column output.
        if (params.kraken2_save_minimizers) {
            ch_kraken2_output = KRAKEN2_STANDARD_REPORT(ch_kraken2_output).report
        }

        // Extract the database name to combine by.
        ch_bracken_databases = databases
            .filter { meta, db -> meta['tool'] == 'bracken' }
            .map { meta, db -> [meta['db_name'], meta, db] }

        // Combine back with the reads
        ch_input_for_bracken = ch_kraken2_output
            .map { meta, report -> [meta['db_name'], meta, report] }
            .combine(ch_bracken_databases, by: 0)
            .map {

                key, meta, reads, db_meta, db ->
                    def db_meta_new = db_meta.clone()

                    // Have to pick second element if using bracken, as first element
                    // contains kraken parameters
                    if ( db_meta['tool'] == 'bracken' ) {

                        // Only take second element if one exists
                        def parsed_params = db_meta_new['db_params'].split(";")
                        if ( parsed_params.size() == 2 ) {
                            db_meta_new['db_params'] =  parsed_params[1]
                        } else {
                            db_meta_new['db_params'] = ""
                        }

                    } else {
                        db_meta_new['db_params']
                    }

                [ key, meta, reads, db_meta_new, db ]
            }
            .multiMap { key, meta, report, db_meta, db ->
                report: [meta + db_meta, report]
                db: db
            }

        BRACKEN_BRACKEN(ch_input_for_bracken.report, ch_input_for_bracken.db)
        ch_versions     = ch_versions.mix(BRACKEN_BRACKEN.out.versions.first())
        ch_raw_profiles = ch_raw_profiles.mix(BRACKEN_BRACKEN.out.reports)

    }
    */

    if ( params.profiler == "centrifuge" ) {

        ch_input_for_centrifuge =  ch_input_for_profiling
                                .filter{
                                    if (it[0].is_fasta) log.warn "[nf-core/taxprofiler] Centrifuge currently does not accept FASTA files as input. Skipping Centrifuge for sample ${it[0].id}."
                                    !it[0].is_fasta
                                }
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }

        CENTRIFUGE_CENTRIFUGE ( ch_input_for_centrifuge.reads, ch_input_for_centrifuge.db, params.centrifuge_save_reads, params.centrifuge_save_reads, params.centrifuge_save_reads  )
        CENTRIFUGE_KREPORT (CENTRIFUGE_CENTRIFUGE.out.report, ch_input_for_centrifuge.db)
        ch_versions            = ch_versions.mix( CENTRIFUGE_CENTRIFUGE.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( CENTRIFUGE_CENTRIFUGE.out.results )
        ch_raw_profiles        = ch_raw_profiles.mix( CENTRIFUGE_KREPORT.out.kreport )
        ch_multiqc_files       = ch_multiqc_files.mix( CENTRIFUGE_KREPORT.out.kreport )

    }
    
    if ( params.profiler == "metaphlan4" ) {

        ch_input_for_metaphlan4 = ch_input_for_profiling
                            .filter{
                                if (it[0].is_fasta) log.warn "[Zymo-Research/aladdin-shotgun] MetaPhlAn4 currently does not accept FASTA files as input. Skipping MetaPhlAn4 for sample ${it[0].id}."
                                !it[0].is_fasta
                            }
                            .multiMap {
                                it ->
                                    reads: [it[0] + it[2], it[1]]
                                    db: it[3]
                            }

        METAPHLAN4_METAPHLAN4 ( ch_input_for_metaphlan4.reads, ch_input_for_metaphlan4.db )
        ch_versions        = ch_versions.mix( METAPHLAN4_METAPHLAN4.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( METAPHLAN4_METAPHLAN4.out.profile )

        METAPHLAN4_QIIMEPREP ( METAPHLAN4_METAPHLAN4.out.profile )
        ch_versions     = ch_versions.mix( METAPHLAN4_QIIMEPREP.out.versions.first() )
        ch_qiime_profiles = ch_qiime_profiles.mix( METAPHLAN4_QIIMEPREP.out.mpa_biomprofile )
        ch_taxonomy = ch_taxonomy.mix( METAPHLAN4_QIIMEPREP.out.taxonomy )

        // METAPHLAN4_UNMAPPED ( METAPHLAN4_QIIMEPREP.out.mpa_info.collect() )
        // ch_multiqc_files = ch_multiqc_files.mix( METAPHLAN4_UNMAPPED.out.json.ifEmpty([]))

    }

    if ( params.profiler == "sourmash" ) {

        ch_input_for_sourmash =  ch_input_for_profiling
                                .filter{
                                    if (it[0].is_fasta) log.warn "[Zymo-Research/aladdin-shotgun] Sourmash currently does not accept FASTA files as input. Skipping Sourmash for sample ${it[0].id}."
                                    !it[0].is_fasta
                                }
                                .multiMap {
                                    it ->
                                        reads: [ it[0] + it[2], it[1] ]
                                        db: it[3]
                                }
        host_lineage = params.host_lineage ? Channel.fromPath(params.host_lineage) : Channel.empty()
        // Temporary place holder for host lineage file until reconfiguration of database into a config file

        if (params.sourmash_trim_low_abund) {
            KHMER_TRIM_LOW_ABUND ( ch_input_for_sourmash.reads )
            ch_input_for_sourmash_sketch = KHMER_TRIM_LOW_ABUND.out.reads
            ch_versions = ch_versions.mix( KHMER_TRIM_LOW_ABUND.out.versions.first() )
        } else {
            ch_input_for_sourmash_sketch = ch_input_for_sourmash.reads
        }

        SOURMASH_SKETCH ( ch_input_for_sourmash.reads )
        ch_versions = ch_versions.mix( SOURMASH_SKETCH.out.versions.first() )
        SOURMASH_GATHER ( SOURMASH_SKETCH.out.sketch , ch_input_for_sourmash.db )
        ch_versions = ch_versions.mix( SOURMASH_GATHER.out.versions.first() )
        SOURMASH_GATHER.out.gather
            .join( SOURMASH_SKETCH.out.sketch )
            .map { [it[0], it[1], it[3]] }
            .set { qiime2_input }
        SOURMASH_QIIMEPREP ( qiime2_input, host_lineage.collect().ifEmpty([]) )
        ch_versions = ch_versions.mix( SOURMASH_QIIMEPREP.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( SOURMASH_QIIMEPREP.out.mqc.collect().ifEmpty([]) )
        ch_qiime_profiles = ch_qiime_profiles.mix( SOURMASH_QIIMEPREP.out.biom )
        ch_taxonomy = ch_taxonomy.mix( SOURMASH_QIIMEPREP.out.taxonomy )

    }    

    if ( params.profiler == "kaiju" ) {

        ch_input_for_kaiju = ch_input_for_profiling
                            .multiMap {
                                it ->
                                    reads: [it[0] + it[2], it[1]]
                                    db: it[3]
                            }

        KAIJU_KAIJU ( ch_input_for_kaiju.reads, ch_input_for_kaiju.db )
        ch_versions = ch_versions.mix( KAIJU_KAIJU.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KAIJU_KAIJU.out.results )

        KAIJU_KAIJU2TABLE_SINGLE ( KAIJU_KAIJU.out.results, ch_input_for_kaiju.db, params.kaiju_taxon_rank)
        ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE_SINGLE.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE_SINGLE.out.summary )
        ch_raw_profiles    = ch_raw_profiles.mix( KAIJU_KAIJU2TABLE_SINGLE.out.summary )
    }

    if ( params.profiler == "diamond" ) {
        if (params.diamond_save_reads) log.warn "[nf-core/taxprofiler] DIAMOND only allows output of a single format. As --diamond_save_reads supplied, only aligned reads in SAM format will be produced, no taxonomic profiles will be available."

        ch_input_for_diamond = ch_input_for_profiling
                                .multiMap {
                                    it ->
                                        reads: [it[0] + it[2], it[1]]
                                        db: it[3]
                                }

        // diamond only accepts single output file specification, therefore
        // this will replace output file!
        ch_diamond_reads_format = params.diamond_save_reads ? 'sam' : params.diamond_output_format

        DIAMOND_BLASTX ( ch_input_for_diamond.reads, ch_input_for_diamond.db, ch_diamond_reads_format , [] )
        ch_versions        = ch_versions.mix( DIAMOND_BLASTX.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( DIAMOND_BLASTX.out.tsv )
        ch_multiqc_files   = ch_multiqc_files.mix( DIAMOND_BLASTX.out.log )

    }

    if ( params.profiler == "motus" ) {

        ch_input_for_motus = ch_input_for_profiling
                                .filter{
                                    if (it[0].is_fasta) log.warn "[nf-core/taxprofiler] mOTUs currently does not accept FASTA files as input. Skipping mOTUs for sample ${it[0].id}."
                                    !it[0].is_fasta
                                }
                                .multiMap {
                                    it ->
                                        reads: [it[0] + it[2], it[1]]
                                        db: it[3]
                                }

        MOTUS_PROFILE ( ch_input_for_motus.reads, ch_input_for_motus.db )
        ch_versions        = ch_versions.mix( MOTUS_PROFILE.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( MOTUS_PROFILE.out.out )
        ch_multiqc_files   = ch_multiqc_files.mix( MOTUS_PROFILE.out.log )
    }

    if ( params.profiler == "krakenuniq" ) {
        ch_input_for_krakenuniq =  ch_input_for_profiling
                                    .map {
                                        meta, reads, db_meta, db ->
                                            [[id: db_meta.db_name, single_end: meta.single_end], reads, db_meta, db]
                                    }
                                    .groupTuple(by: [0,2,3])
                                    .multiMap {
                                        single_meta, reads, db_meta, db ->
                                            reads: [ single_meta + db_meta, reads.flatten() ]
                                            db: db
                                }
        // Hardcode to _always_ produce the report file (which is our basic output, and goes into)
        KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( ch_input_for_krakenuniq.reads, ch_input_for_krakenuniq.db, params.krakenuniq_ram_chunk_size, params.krakenuniq_save_reads, true, params.krakenuniq_save_readclassifications )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
        ch_versions            = ch_versions.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.classified_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )

    }

    // Find samples that failed profiling and output a warning
    ch_qiime_profiles
        .map { [it[0].id, it[1]] }
        .set { ch_profiling_pass }
    reads
        .map { [it[0].id, it[1]] }
        .join(ch_profiling_pass, remainder:true)
        .filter { !it[2] }
        .map { it[0] }
        .collect()
        .map {
            "The following samples failed taxonomy profiling steps or didn't have any identifiable microbe:\n${it.join("; ")}\nPlease contact us if you want to troubleshoot them."
        }
        .set {ch_warning_message }
    ch_warning_message
        .subscribe {
            log.error "$it"
            params.ignore_failed_samples ? { log.warn "Ignoring failed samples and continue!" } : System.exit(1)
        }

    emit:
    classifications = ch_raw_classifications
    profiles        = ch_raw_profiles    // channel: [ val(meta), [ reads ] ] - should be text files or biom
    qiime_profiles  = ch_qiime_profiles  // channel: [ val(meta), absolute abundance profiles ]
    qiime_taxonomy  = ch_taxonomy
    versions        = ch_versions
    motus_version   = params.profiler == "motus" ? MOTUS_PROFILE.out.versions.first() : Channel.empty()
    mqc             = ch_multiqc_files
    warning         = ch_warning_message
}
