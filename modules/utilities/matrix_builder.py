"""
matrix_builder.py

Description:
This script builds a gene expression matrix from multiple sample files.
- It takes a reference gene name list and a file listing sample file paths.
- For each gene in the reference, it finds the corresponding count in each sample file.
- If a gene is missing in a sample, '0' is inserted.
- Output is a tab-separated matrix: rows are genes, columns are samples.

Usage:
python matrix_builder.py gene_names.txt sample_files_list.txt
"""

import sys

def main(gene_names_file, sample_files_list):
    """
    Build gene expression matrix from sample files.
    
    Args:
        gene_names_file: Path to file containing gene names (one per line)
        sample_files_list: Path to file containing sample file paths (one per line)
    """
    # Read gene names
    with open(gene_names_file) as f:
        gene_names = [line.strip() for line in f]

    # Read sample file paths
    with open(sample_files_list) as f:
        sample_files = [line.strip() for line in f if line.strip()]

    # Build a dict for each sample: {gene_name: count}
    # StringTie uses transcript IDs like SMEL4.1_01g000730.1.01
    # Gene groups CSV uses gene IDs like SMEL4.1_01g000730
    # Need to strip transcript suffix to match
    sample_dicts = []
    for sample_file in sample_files:
        gene_to_count = {}
        with open(sample_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                gene_full = parts[0]
                count = parts[1]
                # Extract base gene ID (remove .1.01 transcript suffix)
                # Pattern: SMEL4.1_XXgYYYYYY.1.01 -> SMEL4.1_XXgYYYYYY
                gene_base = gene_full.rsplit('.', 2)[0] if gene_full.count('.') >= 2 else gene_full
                # Store with base gene ID (may overwrite if multiple transcripts, keep last)
                gene_to_count[gene_base] = count
                # Also keep original for exact matches
                gene_to_count[gene_full] = count
        sample_dicts.append(gene_to_count)

    # Output matrix
    for gene in gene_names:
        row = [gene]
        for sample in sample_dicts:
            count = sample.get(gene, "0")
            row.append(count)
        print("\t".join(row))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python matrix_builder.py gene_names.txt sample_files_list.txt")
        sys.exit(1)
    gene_names_file = sys.argv[1]
    sample_files_list = sys.argv[2]
    main(gene_names_file, sample_files_list)
