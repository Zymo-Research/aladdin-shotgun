#!/usr/bin/env python

import sys
import argparse
import re
import pandas as pd

def parse_args(args=None):
    Description = "Examine qiime filtered counts file and identify COVID samples"

    Epilog = "Example usage: id_covid_samples.py <FILE_IN> <DETECTION_THRESHOLD>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input qiime filtered counts file.")
    parser.add_argument("DETECTION_THRESHOLD", type=float, help="COVID detection threshold")    
    return parser.parse_args(args)

def filter_samples(file_in, threshold):

    df = pd.read_csv(file_in, sep='\t',skiprows=[0],index_col=False)
    target_row = df[df.iloc[:, 0].str.contains('Severe acute respiratory syndrome-related coronavirus')]

    if target_row.empty:
        print("No samples detected for 'Severe acute respiratory syndrome-related coronavirus'")
        return
    else:
        target_values = target_row.iloc[0, 1:].astype(float)
        filtered_columns = target_values[target_values >= threshold].index
        for column in filtered_columns:
            print(column)
        return 0

def main(args=None):
    args = parse_args(args)
    filter_samples(args.FILE_IN, args.DETECTION_THRESHOLD)

if __name__ == "__main__":
    sys.exit(main())