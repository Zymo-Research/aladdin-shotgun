#!/usr/bin/env python3

import pandas as pd
import argparse

def add_prefix_to_lineage(lineage):
    # This function is meant to catch scenarios where a databases does not have these prefixes at each taxonomy level
    # This code will add prefixes to taxa lines without for consistency
    prefixes = ['d__','p__','c__','o__','f__','g__','s__']
    levels = lineage.split(';')
    parsed_lineage = []
    for idx,lvl in enumerate(levels):
        if idx<len(prefixes) and not lvl.startswith(prefixes[idx]):
            if lvl == "__":
                parsed_lineage.append(prefixes[idx])
            else:
                parsed_lineage.append(prefixes[idx]+lvl)
        else:
            parsed_lineage.append(lvl)
    return ';'.join(parsed_lineage)


def parse_for_krona(data, name):
    data["lineage"] = data["lineage"].map(add_prefix_to_lineage)
    new_cols = data["lineage"].str.split(";", expand=True, n=6)
    df = pd.concat([data[["percentage"]], new_cols], axis=1)
    df.to_csv("./{}.txt".format(name), sep="\t", index=False, header=False)


def separate_samples_forkrona(reltables):
    #table needs to be separated by sample before krona can accept
    reltable = pd.read_csv(reltables, index_col="#OTU ID")
    for sample in reltable.columns:
        sampletable = pd.DataFrame(reltable.loc[:,sample])
        sampletable["lineage"] = sampletable.index
        sampletable.columns = ["percentage", "lineage"]
        parse_for_krona(sampletable, sample)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse relative taxonomy results for input to krona""")
    parser.add_argument("reltables", type=str, help="relative taxonomy file at species level")
    args = parser.parse_args()
    separate_samples_forkrona(args.reltables)
    
