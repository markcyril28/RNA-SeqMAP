#!/usr/bin/env python3
# Extracts gene ID and Name from Eggplant FASTA headers to CSV

import re
import csv
import sys

# Pre-compiled regex — avoids re-compilation per FASTA header line. O(H) total
# instead of O(H × compile_cost) where H = number of header lines.
_NAME_PATTERN = re.compile(r'Name:"([^"]*)"')

def extract_gene_info(input_fasta, output_csv):
    with open(input_fasta, 'r') as fasta, open(output_csv, 'w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(["Gene_ID", "Name"])

        for line in fasta:
            if line.startswith('>'):
                # Extract gene ID (first field after ">")
                gene_id = line.split()[0][1:]
                # Extract Name from Name:"..." pattern
                name_match = _NAME_PATTERN.search(line)
                name = name_match.group(1) if name_match else ""
                writer.writerow([gene_id, name])

    print(f"Extracted to {output_csv}")

if __name__ == '__main__':
    if len(sys.argv) == 3:
        extract_gene_info(sys.argv[1], sys.argv[2])
    else:
        # Default behavior for backward compatibility
        extract_gene_info(
            "Eggplant_V4.1_transcripts.function.fa",
            "Eggplant_V4.1_transcripts.function.gene_info.csv"
        )
