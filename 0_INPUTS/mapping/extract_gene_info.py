#!/usr/bin/env python3
# Extracts gene ID and Name from Eggplant FASTA headers to CSV

import re, csv

input_fasta = "Eggplant_V4.1_transcripts.function.fa"
output_csv = "Eggplant_V4.1_transcripts.function.gene_info.csv"

with open(input_fasta, 'r') as fasta, open(output_csv, 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(["Gene_ID", "Name"])
    
    for line in fasta:
        if line.startswith('>'):
            # Extract gene ID (first field after ">")
            gene_id = line.split()[0][1:]
            # Extract Name from Name:"..." pattern
            name_match = re.search(r'Name:"([^"]*)"', line)
            name = name_match.group(1) if name_match else ""
            writer.writerow([gene_id, name])

print(f"Extracted to {output_csv}")
