#!/usr/bin/env python
import pandas as pd
import argparse
import re
import json
import unicodedata

#the following function takes in all taxonomy files and create a non-redundant taxonomy file
def group_interest_comp(excel, level6, level7):
    genuscsv = pd.read_csv(level6)
    speciescsv = pd.read_csv(level7)

    #For every row in genus csv and species csv, shave OTU ID down to genus
    for place in genuscsv.index:
        genus = re.search(r'g__(.*)', genuscsv["#OTU ID"][place])
        genus = genus.group(1)
        genuscsv.at[place, "#OTU ID"] = genus

    for place in speciescsv.index:
        searchterm = re.search(r'g__(.*);s__(.*)', speciescsv["#OTU ID"][place])
        species = searchterm.group(2)
        speciescsv.at[place, "#OTU ID"] = species

    with pd.ExcelFile(excel) as xls:
        sheetnames = xls.sheet_names

        for sheet,num in zip(sheetnames, range(1,len(sheetnames))):
            df = pd.read_excel(xls, sheet)
            df = df[["ID", "type", "Other_name"]]
            #for every sheet, place genus specific selections and species specific selections into a list
            genuslist = df[df["type"] == "genus"]
            specieslist = df[df["type"] == "species"]
            #specieslist = list(specieslist["ID"]) 

            othernamegenus = []
            for row in genuslist["Other_name"].dropna():
                genusbyrow = row.rsplit(";")
                othernamegenus = othernamegenus+genusbyrow
            genuslist = list(genuslist["ID"])
            genuslist = genuslist+othernamegenus

            othernamespecies = []
            for row in specieslist["Other_name"].dropna():
                speciesbyrow = row.rsplit(";")
                othernamespecies = othernamespecies+speciesbyrow
            specieslist = list(specieslist["ID"])
            specieslist = specieslist + othernamespecies

            genuslist = '|'.join(genuslist)
            specieslist = '|'.join(specieslist)
            genuslist = unicodedata.normalize("NFKD", genuslist)
            specieslist = unicodedata.normalize("NFKD", specieslist)

            if genuslist=="":
                csvhalf1_genus = genuscsv.iloc[:0,:].copy()
            else:
                csvhalf1_genus = genuscsv[genuscsv["#OTU ID"].str.contains(genuslist)]
            if specieslist =="":
                csvhalf2_species = speciescsv.iloc[:0,:].copy()
            else:
                csvhalf2_species = speciescsv[speciescsv["#OTU ID"].str.contains(specieslist)]

            fullcsv = pd.concat([csvhalf1_genus,csvhalf2_species], axis = 0)
            fullcsv.index = fullcsv["#OTU ID"]
            fullcsv = fullcsv.drop(columns = ["#OTU ID"])

            sheet = re.sub(" ", "_", sheet)
            fullcsv.to_csv("./"+sheet+"_groupinterest_comp.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Tally genus and species composition within groups of interest""")
    parser.add_argument("-e", "--excel", dest="excel", type=str, help="excel file containing specified species and genus. Can have multiple tabs for multiple groups")
    parser.add_argument("-g", "--genus", dest="level6", type=str, help="genus percentage count file")
    parser.add_argument("-s", "--species", dest="level7", type=str, help="species percentage count file")
    args = parser.parse_args()
    group_interest_comp(args.excel, args.level6, args.level7)
