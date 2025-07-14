#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import re

def qiimebarplot_tomultiqc(table, levelnum):
    f = pd.read_csv(table, sep=",")
    f=f.transpose()
    #Drop first unecessary row
    f.columns = f.iloc[0]
    f = f.iloc[1:,:]
    #remove group row for purposes of generating barplot
    f = f.iloc[:-1,:]
    #set column count to percentage
    for col in f.columns:
        f[col] = (f[col]/f[col].sum())*100
    f.insert(0, "#OTU ID", f.index)

    f.to_csv("./level-"+levelnum+".csv", index=False, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse metaphlan table""")
    parser.add_argument("-d", "--dataframe", dest="table", type=str, help="output csv from qiime barplot export")
    parser.add_argument("-l", "--levelnum", dest="levelnum", type=str, help="level number")
    args = parser.parse_args()
    qiimebarplot_tomultiqc(args.table, args.levelnum)
