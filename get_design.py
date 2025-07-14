"""
This script scans a directory for FASTQ files and generates a CSV file with
sample info.

Steps:
1. Define paths for the input FASTQ directory and the output CSV file.
2. Initialize a dictionary to store sample info.
3. Define regex to match paired-end and single-end FASTQ files.
4. Scan the directory for FASTQ files and populate the dictionary with sample
info.
5. Write the sample info to a CSV file with columns: sample, read_1, read_2,
group, and run_accession.

Variables:
    FASTQ_DIR (str): Path to the directory containing FASTQ files.
    OUTPUT_CSV (str): Path to the output CSV file.
    samples (dict): Dictionary to store sample info.
    PAIRED_PATTERN (re.Pattern): Regex to match paired-end FASTQ files.
    SINGLE_PATTERN (re.Pattern): Regex to match single-end FASTQ files.

Functions:
    None

Usage:
    Run the script to generate a CSV file with sample info from the specified
    FASTQ directory.
"""
import os
import csv
import re

# Define paths
FASTQ_DIR = "/data/namlhs/fastq-gmwi2"
OUTPUT_CSV = "/data/namlhs/fastq-gmwi2/sample_list.csv"

# Dictionary to store sample information
samples = {}

# Regular expression to match read files (handles both paired-end & single-end)
PAIRED_PATTERN = re.compile(r"(.+?)_(1|2)\.fastq\.gz$")
SINGLE_PATTERN = re.compile(r"(.+?)\.fastq\.gz$")  # Single-end files

# Scan directory
for file in os.listdir(FASTQ_DIR):
    file_path = os.path.join(FASTQ_DIR, file)

    paired_match = PAIRED_PATTERN.match(file)
    single_match = SINGLE_PATTERN.match(file)

    if paired_match:
        sample, read_pair = paired_match.groups()
        if sample not in samples:
            samples[sample] = {"read_1": "", "read_2": ""}

        if read_pair == "1":
            samples[sample]["read_1"] = file_path
        elif read_pair == "2":
            samples[sample]["read_2"] = file_path

    elif single_match:
        sample = single_match.group(1)
        if sample not in samples:
            samples[sample] = {"read_1": file_path, "read_2": ""}

# Write to CSV
with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["sample", "read_1", "read_2", "group", "run_accession"])

    for sample, paths in samples.items():
        writer.writerow([sample, paths["read_1"], paths["read_2"], "A", ""])

print(f"CSV file created at {OUTPUT_CSV}")
