/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTaxprofiler.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.genome,
                            params.outdir, params.longread_hostremoval_index,
                            params.hostremoval_reference, params.shortread_hostremoval_index,
                            params.multiqc_config, params.shortread_qc_adapterlist,
                            params.krona_taxonomy_directory,
                            params.taxpasta_taxonomy_dir,
                            params.multiqc_logo, params.multiqc_methods_description
                        ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if ( params.input ) {
    ch_input              = file(params.input, checkIfExists: true)
} else {
    exit 1, "Input samplesheet not specified"
}

if (params.database) {
    if (!params.databases.containsKey(params.database)) {
        exit 1, "The provided database '${params.database}' is not available. Currently the available databases are ${params.databases.keySet().join(',')}"
    }
    else {
        // Make a db_meta map that is consistent with previous versions of pipeline
        // For others, since we are only allowing one profiler at a time, better to have them in params than channel
        params.db_meta = params.databases[params.database].subMap(['tool','db_param']) + ['db_name': params.database]
        params.profiler = params.databases[params.database].tool
        params.db_path = params.databases[params.database].db_path
        params.host_lineage = params.databases[params.database].host_lineage ?: false
    }
} else {
    exit 1, "Database not specified!"
}

if (params.shortread_qc_includeunmerged && !params.shortread_qc_mergepairs) exit 1, "ERROR: [nf-core/taxprofiler] cannot include unmerged reads when merging is not turned on. Please specify --shortread_qc_mergepairs"

if (params.shortread_complexityfilter_tool == 'fastp' && ( params.perform_shortread_qc == false || params.shortread_qc_tool != 'fastp' ))  exit 1, "ERROR: [nf-core/taxprofiler] cannot use fastp complexity filtering if preprocessing not turned on and/or tool is not fastp. Please specify --perform_shortread_qc and/or --shortread_qc_tool 'fastp'"

if (params.perform_shortread_hostremoval && !params.hostremoval_reference) { exit 1, "ERROR: [nf-core/taxprofiler] --shortread_hostremoval requested but no --hostremoval_reference FASTA supplied. Check input." }
if (!params.hostremoval_reference && params.shortread_hostremoval_index) { exit 1, "ERROR: [nf-core/taxprofiler] --shortread_hostremoval_index provided but no --hostremoval_reference FASTA supplied. Check input." }

if (params.hostremoval_reference           ) { ch_reference = file(params.hostremoval_reference) }
if (params.shortread_hostremoval_index     ) { ch_shortread_reference_index = Channel.fromPath(params.shortread_hostremoval_index).map{[[], it]} } else { ch_shortread_reference_index = [] }
if (params.longread_hostremoval_index      ) { ch_longread_reference_index  = file(params.longread_hostremoval_index     ) } else { ch_longread_reference_index  = [] }

if (params.profiler=='malt' && params.run_krona && !params.krona_taxonomy_directory) log.warn "[nf-core/taxprofiler] Krona can only be run on MALT output if path to Krona taxonomy database supplied to --krona_taxonomy_directory. Krona will not be executed in this run for MALT."
//Not supporting sequential kraken2-bracken yet
//if (params.run_bracken && !params.run_kraken2) exit 1, 'ERROR: [nf-core/taxprofiler] You are attempting to run Bracken without running kraken2. This is not possible! Please set --run_kraken2 as well.'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def mqcPlugins = Channel.fromPath("${baseDir}/assets/mqc_plugins/", checkIfExists: true)

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                   } from '../subworkflows/local/input_check'

include { SHORTREAD_PREPROCESSING       } from '../subworkflows/local/shortread_preprocessing'
include { LONGREAD_PREPROCESSING        } from '../subworkflows/local/longread_preprocessing'
include { SHORTREAD_HOSTREMOVAL         } from '../subworkflows/local/shortread_hostremoval'
include { LONGREAD_HOSTREMOVAL          } from '../subworkflows/local/longread_hostremoval'
include { SHORTREAD_COMPLEXITYFILTERING } from '../subworkflows/local/shortread_complexityfiltering'
include { PROFILING                     } from '../subworkflows/local/profiling'
include { DIVERSITY                     } from '../subworkflows/local/diversity'
include { VISUALIZATION_KRONA           } from '../subworkflows/local/visualization_krona'
include { STANDARDISATION_PROFILES      } from '../subworkflows/local/standardisation_profiles'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR                       } from '../modules/nf-core/untar/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { FALCO                       } from '../modules/nf-core/falco/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TAXPROFILER {

    ch_versions = Channel.empty()
    ch_multiqc_logo= Channel.fromPath("$projectDir/docs/images/nf-core-taxprofiler_logo_custom_light.png")
    ch_multiqc_files = Channel.empty()
    adapterlist = params.shortread_qc_adapterlist ? file(params.shortread_qc_adapterlist) : []

    if ( params.shortread_qc_adapterlist ) {
        if ( params.shortread_qc_tool == 'adapterremoval' && !(adapterlist.extension == 'txt') ) error "[nf-core/taxprofiler] ERROR: AdapterRemoval2 adapter list requires a `.txt` format and extension. Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
        if ( params.shortread_qc_tool == 'fastp' && !adapterlist.extension.matches(".*(fa|fasta|fna|fas)") ) error "[nf-core/taxprofiler] ERROR: fastp adapter list requires a `.fasta` format and extension (or fa, fas, fna). Check input: --shortread_qc_adapterlist ${params.shortread_qc_adapterlist}"
    }

    /*
        SUBWORKFLOW: Read in samplesheet, validate and stage input files
    */
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Untar tar.gz database file
    ch_db = Channel.from([[params.db_meta, file(params.db_path, checkIfExists: true)]])
    if (params.db_path.endsWith(".tar.gz")) {
        UNTAR (ch_db)
        ch_versions = ch_versions.mix(UNTAR.out.versions.first())
        ch_db = UNTAR.out.untar
    }

    /*
        MODULE: Run FastQC
    */
    ch_input_for_fastqc = INPUT_CHECK.out.fastq.mix( INPUT_CHECK.out.nanopore )

    if ( params.preprocessing_qc_tool == 'falco' ) {
        FALCO ( ch_input_for_fastqc )
        ch_versions = ch_versions.mix(FALCO.out.versions.first())
    } else {
        FASTQC ( ch_input_for_fastqc )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }
    /*
        SUBWORKFLOW: PERFORM PREPROCESSING
    */

    if ( params.perform_shortread_qc ) {
        ch_shortreads_preprocessed = SHORTREAD_PREPROCESSING ( INPUT_CHECK.out.fastq, adapterlist ).reads
        ch_versions = ch_versions.mix( SHORTREAD_PREPROCESSING.out.versions )
    } else {
        ch_shortreads_preprocessed = INPUT_CHECK.out.fastq
    }

    if ( params.perform_longread_qc ) {
        ch_longreads_preprocessed = LONGREAD_PREPROCESSING ( INPUT_CHECK.out.nanopore ).reads
                                        .map { it -> [ it[0], [it[1]] ] }
        ch_versions = ch_versions.mix( LONGREAD_PREPROCESSING.out.versions )
    } else {
        ch_longreads_preprocessed = INPUT_CHECK.out.nanopore
    }

    /*
        SUBWORKFLOW: COMPLEXITY FILTERING
    */

    // fastp complexity filtering is activated via modules.conf in shortread_preprocessing
    if ( params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp' ) {
        ch_shortreads_filtered = SHORTREAD_COMPLEXITYFILTERING ( ch_shortreads_preprocessed ).reads
        ch_versions = ch_versions.mix( SHORTREAD_COMPLEXITYFILTERING.out.versions )
    } else {
        ch_shortreads_filtered = ch_shortreads_preprocessed
    }

    /*
        SUBWORKFLOW: HOST REMOVAL
    */

    if ( params.perform_shortread_hostremoval ) {
        if ( !params.host_lineage ) {
            ch_shortreads_hostremoved = SHORTREAD_HOSTREMOVAL ( ch_shortreads_filtered, ch_reference, ch_shortread_reference_index ).reads
            ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_HOSTREMOVAL.out.mqc.collect{it[1]}.ifEmpty([]))
            ch_versions = ch_versions.mix(SHORTREAD_HOSTREMOVAL.out.versions)
        } else {
            log.warn "Host removal requested but profiler database already include host. Host removal not performed!"
            ch_shortreads_hostremoved = ch_shortreads_filtered
        }
    } else {
        ch_shortreads_hostremoved = ch_shortreads_filtered
    }

    if ( params.perform_longread_hostremoval ) {
        ch_longreads_hostremoved = LONGREAD_HOSTREMOVAL ( ch_longreads_preprocessed, ch_reference, ch_longread_reference_index ).reads
        ch_versions = ch_versions.mix(LONGREAD_HOSTREMOVAL.out.versions)
    } else {
        ch_longreads_hostremoved = ch_longreads_preprocessed
    }

    if ( params.perform_runmerging ) {

        ch_reads_for_cat_branch = ch_shortreads_hostremoved
            .mix( ch_longreads_hostremoved )
            .map {
                meta, reads ->
                    def meta_new = meta.clone()
                    meta_new.remove('run_accession')
                    [ meta_new, reads ]
            }
            .groupTuple()
            .map {
                meta, reads ->
                    [ meta, reads.flatten() ]
            }
            .branch {
                meta, reads ->
                // we can't concatenate files if there is not a second run, we branch
                // here to separate them out, and mix back in after for efficiency
                cat: ( meta.single_end && reads.size() > 1 ) || ( !meta.single_end && reads.size() > 2 )
                skip: true
            }

        ch_reads_runmerged = CAT_FASTQ ( ch_reads_for_cat_branch.cat ).reads
            .mix( ch_reads_for_cat_branch.skip )
            .map {
                meta, reads ->
                [ meta, [ reads ].flatten() ]
            }
            .mix( INPUT_CHECK.out.fasta )

        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    } else {
        ch_reads_runmerged = ch_shortreads_hostremoved
            .mix( ch_longreads_hostremoved, INPUT_CHECK.out.fasta )
    }

    /*
        SUBWORKFLOW: PROFILING
    */
    PROFILING ( ch_reads_runmerged, ch_db )
    ch_versions = ch_versions.mix( PROFILING.out.versions )

    /*
        SUBWORKFLOW: DIVERSITY with Qiime2
    */
    DIVERSITY ( PROFILING.out.qiime_profiles, PROFILING.out.qiime_taxonomy, PROFILING.out.qiime_readcount, INPUT_CHECK.out.groups )

    /*
        SUBWORKFLOW: VISUALIZATION_KRONA
    */
    if ( params.run_krona ) {
        VISUALIZATION_KRONA ( PROFILING.out.classifications, PROFILING.out.profiles, ch_db )
        ch_versions = ch_versions.mix( VISUALIZATION_KRONA.out.versions )
    }

    /*
        SUBWORKFLOW: PROFILING STANDARDISATION
    */
    if ( params.run_profile_standardisation ) {
        STANDARDISATION_PROFILES ( PROFILING.out.classifications, PROFILING.out.profiles, ch_db, PROFILING.out.motus_version )
        ch_versions = ch_versions.mix( STANDARDISATION_PROFILES.out.versions )
    }

    /*
        MODULE: MultiQC
    */

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    workflow_summary    = WorkflowTaxprofiler.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowTaxprofiler.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    if ( params.preprocessing_qc_tool == 'falco' ) {
        ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.txt.collect{it[1]}.ifEmpty([]))
    } else {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }


    if (params.perform_shortread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_longread_qc) {
        ch_multiqc_files = ch_multiqc_files.mix( LONGREAD_PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_shortread_complexityfilter && params.shortread_complexityfilter_tool != 'fastp'){
        ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_COMPLEXITYFILTERING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    if (params.perform_longread_hostremoval) {
        ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_HOSTREMOVAL.out.mqc.collect{it[1]}.ifEmpty([]))
    }

    ch_multiqc_files = ch_multiqc_files.mix( PROFILING.out.mqc.collect().ifEmpty([]) )
    ch_multiqc_files = ch_multiqc_files.mix( DIVERSITY.out.mqc.collect().ifEmpty([]) )

    if ( params.run_profile_standardisation ) {
        ch_multiqc_files = ch_multiqc_files.mix( STANDARDISATION_PROFILES.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        mqcPlugins
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
