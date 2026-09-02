#!/usr/bin/env python3
# Extracts gene ID and Name from FASTA headers to CSV.
#
# BANANA TEMPLATE - copied verbatim from the eggplant tree. _NAME_PATTERN below
# matches a Name:"..." field used by the eggplant FASTA headers; Musa headers
# from the Banana Genome Hub do not use that convention, so the Name column will
# be empty until the regex is adapted. Check a header first: head -1 <input.fa>

import re
import csv
import sys
from pathlib import Path

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
        extract_gene_info(str(Path(sys.argv[1]).resolve()), str(Path(sys.argv[2]).resolve()))
    else:
        print(f"Usage: {sys.argv[0]} <input.fasta> <output.csv>", file=sys.stderr)
        print("Example: python3 extract_gene_info.py Musa_acuminata_DH_Pahang_v4_transcripts.fa output.csv", file=sys.stderr)
        sys.exit(1)
