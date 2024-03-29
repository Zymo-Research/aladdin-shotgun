#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse
import json

def mqc_qiimealpha(alphatable):
    alpha = pd.read_csv(alphatable, sep="\t", index_col=0)
    alpha_json = alpha.to_json(orient = "index")
    alpha_parsed = json.loads(alpha_json)

    description = "The following bargraph provides the alpha diversity Shannon index calculated by QIIME for every sample"

    alpha_mqc = {
        'id' : 'individ_alpha_barchart',
        'section_name' : 'Shannon Alpha Diversity',
        'description': description,
        'plot_type'  : 'bargraph',
        'pconfig' : {
            'id' : 'alpha_diversity_barplot',
            'title' : 'Shannon Alpha Diversity Scores',
            'ylab': 'Shannon Index',
            'yDecimals': 'true',
            'tt_decimals': 2
            },
        'categories' : {
            'shannon_entropy' : {
                'name' : 'Shannon Alpha Diversity Index',
                'color': '#3EE4BF'
                }
            }
        }

    alpha_mqc['data'] = alpha_parsed
    with open('alphachart_mqc.json', 'w') as ofh:
        json.dump(alpha_mqc, ofh, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Alpha diversity table generated by qiime""")
    parser.add_argument("-a", "--alphadiversity", dest="alphatsv", type=str, help="qiime output alpha diversity measurements")
    args = parser.parse_args()
    mqc_qiimealpha(args.alphatsv)
