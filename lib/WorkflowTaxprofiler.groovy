//
// This file holds several functions specific to the workflow/taxprofiler.nf in the Zymo-Research/aladdin-shotgun pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowTaxprofiler {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = "    <dl class=\"dl-horizontal\"; style=\"font-size: 13px\">\n"
        def exclude_groups = ['Core Nextflow options', 
                              'Input/output options', 
                              'Institutional config options', 
                              'Max job request options',
                              'Generic options'] // No need to show these in Aladdin reports
        for (group in summary.keySet()) {
            if (!exclude_groups.contains(group)) {
                def group_params = summary.get(group)  // This gets the parameters of that particular group
                if (group_params) {
                    for (param in group_params.keySet()) {
                        summary_section += "        <dt>$param</dt><dd><samp>${group_params[param]==null?'<span style=\"color:#999999;\">N/A</a>':group_params[param]}</samp></dd>\n"
                    }
                }
            }
        }
        summary_section += "    </dl>\n"

        String yaml_file_text  = "id: 'workflow-summary'\n"
        yaml_file_text        += "description: 'This section summarizes important parameters used in the pipeline. Only parameters that differ from the default are shown.'\n"
        yaml_file_text        += "section_name: 'Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }//
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }
}
