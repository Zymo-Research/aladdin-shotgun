#!/usr/bin/env python

import argparse
import pandas as pd
import json

def parse_resistome(resistome_results):
    resistome_mqc = {
        "id": "resistome_class_composition",
        "section_name": "Resistome composition (AMRplusplus)",
        "description": ("This plot depicts the composition of reads of antimicrobial resistance gene classes identified (or not) by "
            "<a href='https://github.com/Microbial-Ecology-Group/AMRplusplus/tree/master'>AMRplusplus</a>. "
            "Antimicrobial resistance analysis is currently only supported for paired end samples and single end samples will not appear in the below plot. "
            "Read counts have been normalized to counts per one million reads that passed trimming filters for each sample. "
            "This plot includes counts of anti-microbial resistant gene classes detected in each shotgun sample. "
            "Each class contains multiple genes. The plotted read counts are the sum of read counts of all genes in that class. "
            "For gene level data, please refer to downloadable table gene_AMR_analytic_matrix. "
            "Percentages in this stacked bargraph represent the share of antimicrobial gene reads for a certain class over the " 
            "total number of antimicrobial gene reads detected in the entire sample"),
        "plot_type": "bargraph",
        "anchor": "resistome_amr_composition",
        "pconfig": {
            "id": "resistome_composition_bargraph",
            "title": "AMR Class composition (AMRplusplus)",
            "cpswitch_c_active": True,
            "cpswitch_counts_label": "Read counts",
            "cpswitch_percent_label": "Relative percentage"
            }
    }

    resistome = pd.read_csv(resistome_results, index_col=0)
    resistome["mean"] = resistome.mean(axis=1)
    resistome = resistome.sort_values(by=["mean"], ascending = False)
    resistome = resistome.drop(columns = ['mean'])
    resistome_json = resistome.to_json()
    resistome_parsed = json.loads(resistome_json)

    resistome_mqc["data"] = resistome_parsed
    with open('class_resistomechart_mqc.json', 'w') as ofh:
        json.dump(resistome_mqc, ofh, indent=4)
    
if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="""Summarize AMR resistome data for plotting at the class level""")
    parser.add_argument("resistome_results", type=str, help="Resistome results on the class level")
    args = parser.parse_args()
    parse_resistome(args.resistome_results)
