#!/usr/bin/env python
import argparse
import pandas as pd
import re
import os
import json
import numpy as np

def datapackage_json(group):
    datapackage = {
        "profile": "tabular-data-package",
        "resources": [
            {
                "name": "lfc_slice",
                "path": "lfc_slice.csv",
                "profile": "tabular-data-resource",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "schema": {
                    "fields": [
                        {
                            "name": "id",
                            "type": "string"
                        },
                        {
                            "name": "(Intercept)",
                            "type": "number"
                        },
                        {
                            "name": group,
                            "type": "number"
                        }
                    ]
                }
            },
            {
                "name": "se_slice",
                "path": "se_slice.csv",
                "profile": "tabular-data-resource",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "schema": {
                    "fields": [
                        {
                            "name": "id",
                            "type": "string"
                        },
                        {
                            "name": "(Intercept)",
                            "type": "number"
                        },
                        {
                            "name": group,
                            "type": "number"
                        }
                    ]
                }
            },
            {
                "name": "w_slice",
                "path": "w_slice.csv",
                "profile": "tabular-data-resource",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "schema": {
                    "fields": [
                        {
                            "name": "id",
                            "type": "string"
                        },
                        {
                            "name": "(Intercept)",
                            "type": "number"
                        },
                        {
                            "name": group,
                            "type": "number"
                        }
                    ]
                }
            },
            {
                "name": "p_val_slice",
                "path": "p_val_slice.csv",
                "profile": "tabular-data-resource",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "schema": {
                    "fields": [
                        {
                            "name": "id",
                            "type": "string"
                        },
                        {
                            "name": "(Intercept)",
                            "type": "number"
                        },
                        {
                            "name": group,
                            "type": "number"
                        }
                    ]
                }
            },
            {
                "name": "q_val_slice",
                "path": "q_val_slice.csv",
                "profile": "tabular-data-resource",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "schema": {
                    "fields": [
                        {
                            "name": "id",
                            "type": "string"
                        },
                        {
                            "name": "(Intercept)",
                            "type": "number"
                        },
                        {
                            "name": group,
                            "type": "number"
                        }
                    ]
                }
            }
        ],
        "metadata": {
            "intercept_groups": {}
        }
    }
    with open("./to_multiqc/"+group+"/datapackage.json", "w") as fh:
        json.dump(datapackage, fh, indent=4)

def select_species(table, indexlist, group, cleangroup):
    sliced_table = table.loc[indexlist, :]
    sliced_table = sliced_table[["(Intercept)", group]]
    sliced_table.columns = ["(Intercept)", cleangroup]
    return sliced_table

def make_heatmap_input(lfc, stderr, qval, wslice, pval, metadata, fdr):
    lfc = pd.read_csv(lfc, index_col = "id")
    se_slice = pd.read_csv(stderr, index_col = "id")
    q_val = pd.read_csv(qval, index_col = "id")
    w_slice = pd.read_csv(wslice, index_col = "id")
    p_val   = pd.read_csv(pval, index_col = "id")
    metadata = pd.read_csv(metadata, sep="\t")

    nonref_groups = []
    
    for group in q_val.columns[1:]:
        cleangroup = re.sub("^group", "", group)
     
        grouponly_qval = q_val
        grouponly_qval = grouponly_qval[grouponly_qval[group] < fdr]
        #grouponly_qval = grouponly_qval.sort_values(by=group, ascending=True)
        grouponly_qval = grouponly_qval[["(Intercept)", group]]
        grouponly_qval.columns = ["(Intercept)", cleangroup]
                            
        grouponly_pval = select_species(p_val, grouponly_qval.index, group, cleangroup)
        grouponly_lfc = select_species(lfc, grouponly_qval.index, group, cleangroup)
        grouponly_se_slice = select_species(se_slice, grouponly_qval.index, group, cleangroup)
        grouponly_w_slice = select_species(w_slice, grouponly_qval.index, group, cleangroup)

        os.mkdir("./to_multiqc/"+cleangroup+"/")
        nonref_groups.append(cleangroup)

        datapackage_json(cleangroup)

        grouponly_qval.to_csv("./to_multiqc/"+cleangroup+"/q_val_slice.csv")
        grouponly_pval.to_csv("./to_multiqc/"+cleangroup+"/p_val_slice.csv")
        grouponly_lfc.to_csv("./to_multiqc/"+cleangroup+"/lfc_slice.csv")
        grouponly_se_slice.to_csv("./to_multiqc/"+cleangroup+"/se_slice.csv")
        grouponly_w_slice.to_csv("./to_multiqc/"+cleangroup+"/w_slice.csv")

    unique_groups = metadata["group"].unique()
    refgroup_position = np.where(~pd.Series(unique_groups).isin(nonref_groups))
    if len(refgroup_position) > 1:
        raise ValueError("More than 1 reference group has been found for this ANCOM-BC group set; this is not possible")
    refgroup = unique_groups[refgroup_position][0]
    with open('refgroup.txt', 'w') as fh:
        fh.write(refgroup)

    for group in q_val.columns[1:]:
        cleangroup = re.sub("^group", "", group)
        group_overview = pd.DataFrame(lfc[group]).merge(pd.DataFrame(p_val[group]), left_index=True, right_index=True, how='inner')
        group_overview = group_overview.merge(pd.DataFrame(q_val[group]), left_index=True, right_index=True, how="inner")
        group_overview = group_overview.merge(pd.DataFrame(se_slice[group]), left_index=True, right_index=True, how="inner")
        group_overview.columns = [cleangroup+"_vs_"+refgroup+"_lfc",cleangroup+"_vs_"+refgroup+"_pval",cleangroup+"_vs_"+refgroup+"_qval", cleangroup+"_vs"+refgroup+"_standard_error"]
        group_overview.to_csv(cleangroup+"_vs_"+refgroup+"_ancombc_group_overview.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse and filter relative abundance data for making heatmap""")
    parser.add_argument("-l", "--lfc", dest="lfc", type=str, help="Metadata file")
    parser.add_argument("-se", "--stderr", dest="stderr", type=str, help="")
    parser.add_argument("-q", "--qval", dest="qval", type=str, help="")
    parser.add_argument("-w", "--wslice", dest="wslice", type=str, help="")
    parser.add_argument("-p", "--pval", dest="pval", type=str, help="")
    parser.add_argument("-f", "--filtered_metadata", dest="metadata", type=str, help="")
    parser.add_argument("-c", "--fdr_cutoff", dest="fdr", type=float, help="")
    parser.add_argument("-t", "--taxa_max_display", dest="taxa_max_display", type=int, help="")
    args = parser.parse_args()
    make_heatmap_input(args.lfc, args.stderr, args.qval, args.wslice, args.pval, args.metadata, args.fdr)
