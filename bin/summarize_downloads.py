#!/usr/bin/env python
###
# This program creates a JSON file that lists the locations of files for download
# and how to display those files on aladdin platform.
###

import argparse
import logging
import os
import json
import csv
from itertools import permutations

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def summarize_downloads(locations, design):

    """
    :param locations: a file containing locations of files on S3
    :param design: the design file containing group and sample labels
    """

    file_info = dict()

    # Read the desgin file to collect valid sample and group labels
    logger.info("Reding design file...")
    groups = set()
    samples = set()
    with open(design, 'r') as fh:
        data = csv.DictReader(fh)
        for row in data:
            samples.add(row['sample'])
            if len(row['group']):
                groups.add(row['group'])

    # Define what to do with each type of files
    categories = {
        'multiqc_report.html'                 : ('Report', 'report'),
        'allsamples_compbarplot.qzv'          : ('Barplot Visualization', 'all_samples'),
        'alpha-rarefaction.qzv'               : ('Alpha Rarefaction Visualization', 'all_samples'),
        # diversity core
        'bray_curtis_emperor.qzv' : ('PCoA Visualization by Bray-Curtis Distance', 'all_samples'),
        'jaccard_emperor.qzv'     : ('PCoA Visualization by Jaccard Distance', 'all_samples'),
        # alpha diversity
        'shannon_vector_alpha.qzv'           : ('Alpha Diversity Box Plot by Shannon Diversity Index', 'all_samples'), 
        'evenness_vector_alpha.qzv'          : ('Alpha Diversity Box Plot by Evenness', 'all_samples'), 
        'observed_features_vector_alpha.qzv' : ('Alpha Diversity Box Plot by Observed Features', 'all_samples'), 
        # beta diversity
        'jaccard_distance_matrix-group.qzv'    : ('Beta Diversity Comparison between Groups by Jaccard Distance', 'all_samples'), 
        'bray_curtis_distance_matrix-group.qzv': ('Beta Diversity Comparison between Groups by Bray-Curtis Distance', 'all_samples'),
        #ANCOM-BC
        'ancombc_visualization.qzv'         : ('ANCOM-BC Taxonomy Level Visualizations', 'comparisons'),
        # groups of interest
        'Groups_of_interest.xlsx': ('Detailed abundances among groups of interest', 'all_samples'),
        # AMR plus plus results
        'class_rawcounts_AMR_analytic_matrix.csv':     ('AMR Class Level Raw Counts Results Across Samples', 'all_samples'),
        'mechanism_rawcounts_AMR_analytic_matrix.csv': ('AMR Mechanism Level Raw Counts Results Across Samples', 'all_samples'),
        'genes_rawcounts_AMR_analytic_matrix.csv':     ('AMR Gene Level Raw Counts Results Across Samples', 'all_samples'),
        'genes_SNPconfirmed_analytic_matrix.csv':      ('AMR SNP Verified Gene Level Count Results Across Samples', 'all_samples'),
        'krona_report.html':                           ('Krona Composition Report', 'all_samples'),
        'abs-abund-table-'                    : ('Read Counts Table for All Taxa Filtered Samples', 'all_samples')
    }

    # Read the file locations
    with open(locations, 'r') as fh:
        for line in fh:
            info = dict()
            path = line.strip()
            info['path'] = path
            logger.info("Processing file {}".format(path))
            # Get file name
            fn = os.path.basename(path)
            # Check each file category
            for key, values in categories.items():
                # if key is suffix
                if fn.endswith(key):
                    file_type, scope = values
                    if "/refmerged/" in path:
                        file_type = file_type + " with Reference Dataset"
                        info_key = "Reference-compare-" + fn
                    else:
                        info_key = fn
                    info['file_type'] = file_type
                    info['scope'] = scope
                    if scope in ['samples', 'comparisons']:
                        sname = fn.replace(key, '')
                        if scope == 'samples':
                            # Check if the parsed sample name is in the original design
                            if sname in samples:
                                info['sample'] = sname
                            else:
                                logger.error("Parsed sample name from {} not found in the design file".format(fn))
                        else:
                            # Check if file name match any possible permutations of group names
                            # Sort group labels by length to avoid detecting partial group names first
                            for g1, g2 in permutations(sorted(groups, key=len, reverse=True), 2):
                                comp_name = "{}_vs_{}".format(g1,g2)
                                if comp_name in sname:
                                    info['comparison'] = comp_name
                                    break
                            else:
                                logger.error("Could not find group names in comparison file {}".format(fn))
                    file_info[info_key] = info
                    break
                elif fn.startswith(key):
                    info_key = fn
                    file_type,scope = values
                    suffix = (fn.split('.')[0]).split('-')[-1]
                    info['file_type'] = file_type + suffix
                    info['scope'] = scope
                    file_info[info_key] = info
                    break
            else:
                logger.error("File {} did not match any expected patterns".format(fn))
    
    # sort file collection
    file_info = dict(sorted(file_info.items(),key=lambda x: str(x).lower().replace("-", "_")))
    # Output the dict to JSON
    with open('files_to_download.json', 'w') as fh:
        json.dump(file_info, fh, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Generate a json file for displaying outputs on aladdin platform""")
    parser.add_argument("locations", type=str, help="File with all the locations of files on S3")
    parser.add_argument("-d", "--design", dest="design", required=True, type=str, help="Design file for sanity check purposes")
    args = parser.parse_args()
    summarize_downloads(args.locations, args.design)    
