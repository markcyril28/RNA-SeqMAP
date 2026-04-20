"""
matrix_builder.py

Description:
This script builds a gene expression matrix from multiple sample files.
- It takes a reference gene name list and a file listing sample file paths.
- For each gene in the reference, it finds the corresponding count in each sample file.
- If a gene is missing in a sample, '0' is inserted.
- Output is a comma-separated matrix: rows are genes, columns are samples.

Usage:
python matrix_builder.py gene_names.txt sample_files_list.txt
"""

import sys
import re
from pathlib import Path

# Pre-compile regex once at module level (avoids re-compilation per call)
_TRAILING_DIGITS_RE = re.compile(r'\.\d+$')
# Track already-warned unparseable values to avoid flooding stderr
_warned_values = set()

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
    # O(S × L × V) where S=samples, L=lines/file, V=ID variants (≤3)
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
                    stripped = _TRAILING_DIGITS_RE.sub('', current)
                    if stripped == current:
                        break  # No more trailing .digits to strip
                    id_variants.append(stripped)
                    current = stripped

                # For each variant, keep the MAX value across duplicate entries
                # Parse float once and store (float, str) tuple to avoid re-parsing
                try:
                    count_float = float(count)
                except (ValueError, TypeError):
                    if count not in _warned_values:
                        _warned_values.add(count)
                        print(f"  Warning: unparseable value '{count}' for gene '{gene_full}' in {sample_file} — defaulting to 0.0",
                              file=sys.stderr)
                    count_float = 0.0
                for gid in id_variants:
                    existing = gene_to_count.get(gid)
                    if existing is None:
                        gene_to_count[gid] = (count_float, count)
                    elif count_float > existing[0]:
                        gene_to_count[gid] = (count_float, count)
        sample_dicts.append(gene_to_count)

    # Build per-sample resolved values using direct lookup.
    # Each sample_dict already has entries keyed by all variant forms (full ID,
    # stripped-once, stripped-twice). We just need to check which gene_names
    # exist as keys. This is O(genes) per sample instead of O(genes × variants).
    #
    # Also handle the reverse case: gene list IDs may have trailing suffixes
    # not present in the abundance data (e.g., CSV has "SMEL5_06g022750.1"
    # but abundance file has "SMEL5_06g022750"). Pre-compute stripped variants
    # of gene_names for fallback lookups.
    gene_name_variants = {}
    for gene in gene_names:
        variants = [gene]
        current = gene
        for _ in range(2):
            stripped = _TRAILING_DIGITS_RE.sub('', current)
            if stripped == current:
                break
            variants.append(stripped)
            current = stripped
        gene_name_variants[gene] = variants

    # O(S × G) per-sample gene resolution with O(1) dict lookups; variant fallback O(V≤3)
    sample_gene_values = []
    for sample in sample_dicts:
        resolved = {}
        for gene in gene_names:
            entry = sample.get(gene)
            if entry is None:
                # Try stripped variants of the gene list ID
                for variant in gene_name_variants[gene][1:]:
                    entry = sample.get(variant)
                    if entry is not None:
                        break
            if entry is not None:
                resolved[gene] = entry[1]
        sample_gene_values.append(resolved)

    # Stream output line-by-line (avoids buffering entire matrix in memory)
    # O(G) write calls instead of O(G×S) — batch each row with join()
    write = sys.stdout.write
    for gene in gene_names:
        write(','.join([gene] + [r.get(gene, '0') for r in sample_gene_values]) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python matrix_builder.py gene_names.txt sample_files_list.txt")
        sys.exit(1)
    # Resolve to absolute paths for portability under Nextflow/Snakemake workDirs
    gene_names_file = str(Path(sys.argv[1]).resolve())
    sample_files_list = str(Path(sys.argv[2]).resolve())
    main(gene_names_file, sample_files_list)
