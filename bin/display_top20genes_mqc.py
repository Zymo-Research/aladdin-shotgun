#!/usr/bin/env python
import argparse
import pandas as pd
import json
from collections import defaultdict

def top20genes(resistome_results):
    genelevel_resistome_mqc = {
        "id": "top20genes_amr_composition",
        "section_name": "Top 20 Genes Composition of Each Sample (AMRplusplus)",
        "description": ("This plot depicts the composition of reads of antimicrobial resistance genes identified by "
            "<a href='https://github.com/Microbial-Ecology-Group/AMRplusplus/tree/master'>AMRplusplus</a>. "
            "It highlights the top 20 genes by mean relative abundance among all AMR genes. "
            "Read counts have been normalized to counts per one million reads that passed trimming filters for each sample. "
            "All other genes are labeled as 'Other'."
            ),
        "plot_type": "bargraph",
        "anchor": "top20genes_amr_composition",
        "pconfig": {
            "id": "top20genes_composition_bargraph",
            "title": "AMR Top 20 Gene composition (AMRplusplus)",
            "cpswitch_c_active": True,
            "cpswitch_counts_label": "Read counts",
            "cpswitch_percent_label": "Relative percentage"
            }
    }

    AMRgene = pd.read_csv(resistome_results, index_col=0)
    AMRgene = AMRgene.reindex(sorted(AMRgene.columns), axis=1)
    AMRperc = pd.DataFrame()

    samplecount_columns = AMRgene.columns
    for column in samplecount_columns:
        AMRperc[column] = (AMRgene[column]/sum(AMRgene[column]))*100
                
    AMRperc["mean"] = AMRperc.mean(axis=1)
    AMRperc = AMRperc.sort_values(by=["mean"], ascending = False)
    AMRgene = AMRgene.reindex(AMRperc.index, axis=0)
    others = AMRgene.iloc[20:,:].sum(axis = 'rows')
    others.name = "Others"
    others = pd.DataFrame(others).T
    AMRgene = AMRgene.iloc[:20,:]
    AMRgene = pd.concat([AMRgene,pd.DataFrame(others)], axis = 0)
    AMRgene_json = AMRgene.to_json()
    AMRgene_parsed = json.loads(AMRgene_json)

    geneorder = list(AMRgene.index)
    colors = ['#a1c9f4','#b9f2f0','#ffb482','#cfcfcf','#fbafe4','#8de5a1','#029e73','#fab0e4','#949494','#ca9161','#d55e00','#d0bbff','#debb9b','#56b4e9','#0173b2','#de8f05','#ece133','#cc78bc','#ff9f9b', '#fffea3', '#999999']
    geneorder_dict = defaultdict()
    for gene,color in zip(geneorder,colors):
        geneorder_dict[gene] = {"color":color}

    genelevel_resistome_mqc["categories"] = geneorder_dict
    genelevel_resistome_mqc["data"] = AMRgene_parsed
    with open('genelevel_resistomechart_mqc.json', 'w') as ofh:
        json.dump(genelevel_resistome_mqc, ofh, indent=4)

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="""Parse resistome results to be top 20 AMR genes expressed in each sample""")
    parser.add_argument("resistome_results", type=str, help="Resistome results on the gene level")
    args = parser.parse_args()
    top20genes(args.resistome_results)
