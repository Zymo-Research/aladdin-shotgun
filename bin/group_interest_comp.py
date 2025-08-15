#!/usr/bin/env python
import pandas as pd
import argparse
import re

def group_interest_comp(interest, comp_file):

    valid_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    # Read in the group of interest EXCEL file first
    groups = pd.read_excel(interest, sheet_name=None)
    compositions = pd.read_csv(comp_file, index_col=0)

    with pd.ExcelWriter("Groups_of_interest.xlsx") as fh:
        # Iterate through each category
        for sheet_name, df in groups.items():
            all_results = []
            # Iterate through each item
            for idx, row in df.iterrows():
                lvl = row["type"].lower()
                if lvl in valid_levels:
                    search_term = lvl[0] + '__' + row["ID"]
                    if isinstance(row["Other_name"],str):
                        names = row["Other_name"].split(';')
                        for n in names:
                            search_term += "|" + lvl[0] + '__' + n
                    search_result = compositions[compositions.index.str.contains(search_term, case=False)].copy()
                    if len(search_result):
                        search_result["Matching_ID"] = row["ID"]
                        all_results.append(search_result)
                else:
                    print("{} is not a recognized taxonomy level in {}. Skipped.".format(lvl, sheet_name))
            # Output detail matches to a EXCEL sheet and sum by entries to a CSV file if anything was found
            if len(all_results):
                df = pd.concat(all_results)
                df.to_excel(fh, sheet_name)
                df = df.groupby("Matching_ID").sum()
                sheet_name = re.sub(r"\s+", "_", sheet_name)
                df.to_csv(sheet_name+"_groupinterest_comp.csv")
            else:
                print("No reads were detected for any entry in the category {}".format(sheet_name))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Tally composition within groups of interest""")
    parser.add_argument("-i", "--interest", dest="interest", type=str, required=True, help="EXCEL file containing specified taxa of interests. Can have multiple tabs for multiple groups")
    parser.add_argument("-c", "--composition", dest="comp_file", type=str, required=True, help="Qiime composition file at the species level")
    args = parser.parse_args()
    group_interest_comp(args.interest, args.comp_file)