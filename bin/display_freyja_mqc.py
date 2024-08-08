#!/usr/bin/env python

import sys
import os
import argparse
import re
import pandas as pd
import numpy as np
import json

def parse_args(args=None):
    Description = "Parse freyja aggregated tsv of samples and output json file for MultiQC."
    Epilog = "Example usage: python display_freyja_mqc.py <freyja_aggregated_tsv>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("freyja_aggregated_tsv", help="Path to freyja aggregated tsv of samples")
    return parser.parse_args(args)

def parse_freyja_tsv(freyja_tsv):
    df = pd.read_csv(freyja_tsv, sep='\t',index_col=0) 
    df = df[df['coverage']>70.0] 

    # Function to parse the summarized column
    def parse_summarized(summarized_str):
        data = ast.literal_eval(summarized_str)
        return {item[0]: item[1] for item in data}

    final_df = df['summarized'].apply(parse_summarized).apply(pd.Series)
    final_df = final_df.fillna(0)

    return final_df

def mqc_freyja(freyja_df):
    freyja_mqc = {
        "id": "covid_strain_analysis",
        "section_name": "Covid strain analysis (freyja)",
        "description": ("This plot provides a fractional abundance estimate of different covid strains for all samples."),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "covid_strain_bargraph",
            "title": "Covid Strains bargraph"
            }
    }

    freyja_json = freyja_df.to_json()
    freyja_parsed = json.loads(freyja_json)

    freyja_mqc["data"] = freyja_parsed
    with open('freyja_bargraph_mqc.json', 'w') as fbh:
        json.dump(freyja_mqc, fbh, indent=4)

def main(args=None):
    args = parse_args(args)
    
    freyja_df = parse_freyja_tsv(args.freyja_aggregated_tsv)
    
    mqc_freyja(freyja_df)

if __name__ == "__main__":
    sys.exit(main())

