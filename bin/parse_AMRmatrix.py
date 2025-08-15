#!/usr/bin/env python

import pandas as pd
import numpy as np
import re
import argparse

def percent_AMRmatrix(AMRmatrix, output):
    AMRmatrix = pd.read_csv(AMRmatrix)
    # Conver to fractions
    num_cols = AMRmatrix.select_dtypes(include=np.number).columns
    AMRmatrix[num_cols] = (AMRmatrix[num_cols].div(AMRmatrix[num_cols].sum()))*100

    AMRmatrix.to_csv(output, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse counts to make top 50 table""")
    parser.add_argument("-i", "--amrinput", dest="AMRmatrix", required=True, help="AMR count matrix")
    parser.add_argument("-o", "--output", dest="output", required=True, help="AMR count matrix top50 output")
    args = parser.parse_args()
    percent_AMRmatrix(args.AMRmatrix, args.output)
