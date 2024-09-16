#!/usr/bin/env python
"""
Convert biom.tsv into an absolute count table for internal use
"""

import argparse
import pandas as pd


def mod_abund_tab(biom_tsv, output):
    """
    Convert biom.tsv into an absolute count table for internal use
    """
    biom_table = pd.read_csv(biom_tsv, delimiter='\t', skiprows=1)
    tab_transposed = biom_table.transpose()
    tab_transposed.reset_index(inplace=True)
    tab_transposed.columns = tab_transposed.iloc[0]
    tab_transposed = tab_transposed.drop(tab_transposed.index[0])
    tab_transposed.rename(columns={tab_transposed.columns[0]: 'query_name'},
                          inplace=True)
    tab_transposed.to_csv(output, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Convert biom.tsv""")
    parser.add_argument("-i", "--input", dest="biom_tsv",
                        type=str, help="biom generated file")
    parser.add_argument("-o", "--output", dest="output",
                        type=str, default="mod_abund_table.csv",
                        help="output file name")
    args = parser.parse_args()
    mod_abund_tab(args.biom_tsv, args.output)
