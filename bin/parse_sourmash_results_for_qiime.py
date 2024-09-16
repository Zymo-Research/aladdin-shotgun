#!/usr/bin/env python

import argparse
import pandas as pd
import re
import json

def add_prefix_to_lineage(lineage):
    # Some of the sourmash databases has these prefixes at each taxonomy level
    # but some doesn't. This code will add prefixes to those without for consistency
    prefixes = ['d__','p__','c__','o__','f__','g__','s__']
    levels = lineage.split(';')
    parsed_lineage = []
    for idx,lvl in enumerate(levels):
        if idx<len(prefixes) and not lvl.startswith(prefixes[idx]):
            parsed_lineage.append(prefixes[idx]+lvl)
        else:
            parsed_lineage.append(lvl)
    return ';'.join(parsed_lineage)

def parse_sourmash(sourmash_results, sketch_log, name, filter_fp, host_lineage):
    # Get sample name if not specified
    if not name:
        name = sourmash_results.replace('.with-lineages.csv','')

    # Initialize dict for multiqc output
    mqc_data = {
        "id": "kmer_composition",
        "section_name": "Kmer composition (sourmash)",
        "description": ("This plot depicts the composition of kmers identified (or not) by "
                        "<a href='https://github.com/sourmash-bio/sourmash'>sourmash</a>. "
                        "It includes percentages of kmers that were identified as common host organisms, "
                        "common eukaryotic pathogens/parasites, and other microbes. Unidentified kmers "
                        "are those that do not match with the database or those with a match but fail to reach "
                        "the base-pair threshold set during the sourmash step of the pipeline. Please refer to "
                        "plots in sections below for detailed compositions of other microbes. "
                        "sourmash only outputs percentages of kmers identified, the numbers of reads you see "
                        "here are estimated using percentages and total numbers of reads that are input into sourmash."),
        "plot_type": "bargraph",
        "anchor": "sourmash_kmer_composition",
        "pconfig": {
            "id": "kmer_composition_bargraph",
            "title": "Kmer composition (sourmash)",
            "cpswitch_c_active": False,
            "cpswitch_counts_label": "Read counts (Estimated)",
            "cpswitch_percent_label": "Percent Kmers"
        },
        "data": { name : dict() }
    }

    # Initialize dict for general stats table
    mqc_gs_data = {
        "id": "kmer_composition_gs",
        "plot_type": "generalstats",
        "pconfig": {
            "identified_kmers": {
                "title": "% Kmers w/ Taxonomy",
                "namespace": "sourmash",
                "description": "% Kmers with assigned taxonomy by sourmash",
                "max": 100,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:.2f}%"
            }
        },
        "data": { name: dict() }
    }

    # Read the sketch log to get number of sequences
    pattern = r"calculated \d+ signature for (\d+) sequences taken from \d+ files"
    with open(sketch_log) as fh:
        t = fh.read()
        m = re.search(pattern, t)
        if m:
            readcount = int(m.group(1)) 
        else:
            raise Exception("Failed to parse sourmash sketch log to get number of reads!")
    
    # Read gather results
    profile = pd.read_csv(sourmash_results)
    # Extract accession number from subject names in the gather results
    profile["accession"] = profile["name"].map(lambda x:x.split()[0])

    # Filter false positives if requested
    if filter_fp:
        profile = profile[(profile["match_containment_ani"] >= 0.935) | (profile["unique_intersect_bp"] > 1000000)]

    # Calculate total no. reads of microbes and unidenfied kmers
    mqc_data["data"][name]["Microbes"] = round(profile["f_unique_weighted"].sum()*readcount, 2)
    mqc_data["data"][name]["Unidentified"] = round(readcount-mqc_data["data"][name]["Microbes"], 2)
    mqc_gs_data["data"][name]["identified_kmers"] = profile["f_unique_weighted"].sum()*100

    # Read host lineage file to get a list of host genome accessions and species names
    if host_lineage:
        host_info = pd.read_csv(host_lineage)
        host_info = host_info[["ident","species"]]  
        # Separate host and other matches 
        host = profile.loc[profile["accession"].isin(host_info["ident"])]
        profile = profile.loc[~profile["accession"].isin(host_info["ident"])]
        mqc_data["data"][name]["Microbes"] = round(profile["f_unique_weighted"].sum()*readcount, 2)
        # Add the species name to the host profile
        host = host[["accession","f_unique_weighted"]]
        host = pd.merge(host, host_info, left_on="accession", right_on="ident", how="left")
        # Record the no. reads of each host organism
        mqc_data["data"][name].update(host.set_index("species")["f_unique_weighted"].mul(readcount).round(2).to_dict())

    if len(profile)>0:
        profile = profile[["accession","lineage","f_unique_weighted"]]
        profile["lineage"] = profile["lineage"].map(add_prefix_to_lineage)
        profile = profile.rename(columns={"f_unique_weighted":name, "lineage":"Taxon", "accession":"Feature ID"})
        # Collapsing rows together in case of duplicate accession (unlikely, but just in case)
        profile = profile.groupby("Feature ID").agg({"Taxon":"first", name:"sum"})
        # Calculate absolute abundance (approximate)
        profile[name] = profile[name].mul(readcount)
        # Create a table file
        profile[name].to_csv(name+"_absabun_parsed_profile.txt", sep="\t")
        # Create a taxonomy file
        profile["Taxon"].to_csv(name+"_profile_taxonomy.txt", sep="\t")
    else:
        print("There is no identified microbial genomes after filtering out the host genomes for sample {}!".format(name))

    # Relative abundance file
    rel_profile = profile
    rel_profile[name] = rel_profile[name]/sum(rel_profile[name])
    # Create a table file for relative abundance
    rel_profile[name].to_csv(name+"_relabun_parsed_profile.txt", sep="\t")

    # Output mqc results
    with open(name+"_sourmash_stats_mqc.json", "w") as fh:
        json.dump(mqc_data, fh, indent=4)
    with open(name+"_sourmash_gs_mqc.json", "w") as fh:
        json.dump(mqc_gs_data, fh, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse sourmash gather results to be compatible with Qiime for taxonomy analysis""")
    parser.add_argument("sourmash_results", type=str, help="Sourmash results after gather & tax annotate")
    parser.add_argument("-l", "--sketch_log", dest="sketch_log", type=str, required=True, help="Sourmash sketch log")
    parser.add_argument("-n", "--name", dest="name", type=str, help="Sample name")
    parser.add_argument("-f", "--filter_false_positive", dest="filter_fp", action="store_true", help="Whether to filter potentially false positive results, empirically")
    parser.add_argument("--host_lineage", dest="host_lineage", type=str, help="Host/pathogen lineage file")
    args = parser.parse_args()
    parse_sourmash(args.sourmash_results, args.sketch_log, args.name, args.filter_fp, args.host_lineage)
