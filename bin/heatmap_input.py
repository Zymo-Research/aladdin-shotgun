#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import re
import seaborn as sns

def extract_taxa(taxa_string, index):
    unidentified_taxa_patterns = [
        r"^$",
        r"metagenome",
        r"uncultured",
        r"unidentified"
    ]
    levels = taxa_string.split(";")
    levels = [re.sub(r"^\w?__", "", l.strip()) for l in levels]
    if levels[index]:
        result = levels[index]
    else:
        result = "sp."
    while any(re.search(p, levels[index], re.IGNORECASE) for p in unidentified_taxa_patterns) and index>0:
        index -= 1
        if levels[index]:
            result = levels[index] + " " + result
    return result

def make_heatmap_input(reltax_tables, metadata, top_n):
    levels = {
        "level-1.csv": ("Kingdom", 0),
        "level-2.csv": ("Phylum", 1),
        "level-3.csv": ("Class", 2),
        "level-4.csv": ("Order", 3),
        "level-5.csv": ("Family", 4),
        "level-6.csv": ("Genus", 5),
        "level-7.csv": ("Species", 6)
    }

    # Read metadata file
    meta = pd.read_csv(metadata, sep="\t")
    for table in reltax_tables:
        try:
            level, index = levels[table]
        except:
            print("{} not recoginzed as one of the relative abundance table!".format(table))
            continue
        data = pd.read_csv(table)
        # Conver to fractions
        num_cols = data.select_dtypes(include=np.number).columns
        data[num_cols] = data[num_cols].div(data[num_cols].sum())
 
        # Merge with group labels
        data_merge = pd.melt(data, id_vars="#OTU ID", var_name="sampleid", value_name="abundance")
        data_merge = data_merge.merge(meta)

        # Get top n most abundant taxa from each group
        top_taxa = data_merge.groupby(["#OTU ID", "group"])["abundance"].mean().groupby("group").nlargest(top_n)
        top_taxa = top_taxa.droplevel(0).reset_index()["#OTU ID"].unique()

        # Subset abundance data for only top n taxa
        data = data.loc[data["#OTU ID"].isin(top_taxa), ]

        # Extract taxonomy for just this level + uppper levels if necessary
        data["#OTU ID"] = data["#OTU ID"].apply(extract_taxa, index=index)

        # Only keep samples that are in the metadata
        data = data.set_index("#OTU ID")
        data = data.loc[:, data.columns.isin(meta["sampleid"])]

        # Clustering (only if at least 5 taxa)
        if len(data) >= 5:
            # Get the lowest count for pseudocount
            if (data>0).all().all():
                pseudocount = 0
            else:
                pseudocount = data.apply(lambda x:x[x>0].min()).min() / 2
            # Log transform (data+pseudocount)
            transformed = np.log10(data+pseudocount)
            fig = sns.clustermap(transformed)
            data = data.iloc[fig.dendrogram_row.reordered_ind, fig.dendrogram_col.reordered_ind]
            data.to_csv("{}_taxo_heatmap.csv".format(level), sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse and filter relative abundance data for making heatmap""")
    parser.add_argument(dest="reltax_tables", nargs="+", type=str, help="Relative abundance tables")
    parser.add_argument("-m", "--metadata", dest="metadata", required=True, help="Metadata file")
    parser.add_argument("-t", "--top_n", dest="top_n", type=int, default=20, help="")
    args = parser.parse_args()
    make_heatmap_input(args.reltax_tables, args.metadata, args.top_n)