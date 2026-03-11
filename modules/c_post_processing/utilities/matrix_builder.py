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
import re

def main(gene_names_file, sample_files_list):
    """
    Build gene expression matrix from sample files.
    
    Args:
        gene_names_file: Path to file containing gene names (one per line)
        sample_files_list: Path to file containing sample file paths (one per line)
    """
    # Read gene names
    with open(gene_names_file) as f:
        gene_names = [line.strip() for line in f if line.strip()]

    # Read sample file paths
    with open(sample_files_list) as f:
        sample_files = [line.strip() for line in f if line.strip()]

    # Build a dict for each sample: {gene_name: count}
    # StringTie uses transcript IDs like SMEL4.1_01g000730.1.01
    # Gene groups CSV uses gene IDs like SMEL4.1_01g000730
    # Need to strip transcript suffix to match
    #
    # DUPLICATE HANDLING: In M2 de novo mode, StringTie can assign multiple
    # "genes" (STRG.*) to the same reference transcript (e.g., different strand
    # calls). For pre-normalized metrics (TPM/FPKM/coverage), we keep the MAX
    # value across duplicates — this represents the dominant isoform's expression.
    # Summing would be incorrect for already-normalized values.
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

                # Collect all gene ID variants (full + stripped suffixes)
                # that this row should be indexed under.
                id_variants = [gene_full]
                # Extract base gene ID by iteratively stripping trailing ".digits"
                # suffixes (max 2 rounds). This safely handles eggplant IDs where
                # the gene ID itself contains a dot (e.g., SMEL4.1_XXgYYYYYY).
                #
                # Example chain for SMEL4.1_06g023900.1.01:
                #   Round 1: strip .01  -> SMEL4.1_06g023900.1  (stored)
                #   Round 2: strip .1   -> SMEL4.1_06g023900    (stored)
                # The chain naturally stops when no trailing .digits remain
                # (e.g., SMEL4.1_06g023900 ends with g023900, not .digits).
                current = gene_full
                for _ in range(2):
                    stripped = re.sub(r'\.\d+$', '', current)
                    if stripped == current:
                        break  # No more trailing .digits to strip
                    id_variants.append(stripped)
                    current = stripped

                # For each variant, keep the MAX value across duplicate entries
                try:
                    count_float = float(count)
                except (ValueError, TypeError):
                    count_float = 0.0
                for gid in id_variants:
                    existing = gene_to_count.get(gid)
                    if existing is None:
                        gene_to_count[gid] = count
                    else:
                        try:
                            existing_float = float(existing)
                        except (ValueError, TypeError):
                            existing_float = 0.0
                        if count_float > existing_float:
                            gene_to_count[gid] = count
        sample_dicts.append(gene_to_count)

    # Output matrix
    for gene in gene_names:
        row = [gene]
        for sample in sample_dicts:
            count = sample.get(gene)
            if count is None:
                # Fallback: strip trailing numeric suffixes from the lookup gene ID
                # e.g. SMEL4.1_XXgYYYYYY.1 -> SMEL4.1_XXgYYYYYY
                gene_base = re.sub(r'\.\d+$', '', gene)
                count = sample.get(gene_base)
            if count is None:
                count = "0"
            row.append(count)
        print("\t".join(row))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python matrix_builder.py gene_names.txt sample_files_list.txt")
        sys.exit(1)
    gene_names_file = sys.argv[1]
    sample_files_list = sys.argv[2]
    main(gene_names_file, sample_files_list)
