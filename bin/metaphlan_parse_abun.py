#!/usr/bin/env python3
import argparse
import pandas as pd

def metaphlan_profileparse(mpa_table, label):
    # 1) Load the full MetaphlAn table
    df = pd.read_csv(
        mpa_table, sep='\t', 
        usecols=['clade_name','clade_taxid','estimated_number_of_reads_from_the_clade']
    )
    df.columns = ['clade_name','clade_taxid','reads']

    # 2) Keep only species-level entries
    is_species = df['clade_name'].str.contains(r'\|s__')
    df = df[is_species].copy()

    # 3) Decide what your "Feature ID" should be:
    #    Option A: entire lineage of taxids, e.g. "2|976|...|816"
    #    Option B: just the final taxid (after last '|')
    # Here we take Option B:
    df['Feature ID'] = df['clade_taxid'].str.split('|').str[-1]

    # 4) Write the abundance table
    out_abun = df[['Feature ID','reads']].copy()
    out_abun.columns = ['Feature ID', label]
    out_abun.to_csv(f"{label}_absabun_parsed_mpaprofile.txt",
                    sep='\t', index=False)

    # 5) Build the taxonomy file for QIIME2
    #    replace '|' → ';' in the clade_name, to match Qiime2 TSVTaxonomyFormat
    tax = df[['Feature ID','clade_name']].copy()
    tax['Taxon'] = tax['clade_name'].str.replace(r'\|', ';', regex=True)
    tax = tax[['Feature ID','Taxon']]
    tax.to_csv(f"{label}_profile_taxonomy.txt",
               sep='\t', index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Parse MetaPhlAn table into feature-counts + taxonomy"
    )
    p.add_argument('-t','--mpa_table', required=True,
                   help="MetaPhlAn output (.txt)")
    p.add_argument('-l','--label',     required=True,
                   help="Sample label prefix")
    args = p.parse_args()
    metaphlan_profileparse(args.mpa_table, args.label)
