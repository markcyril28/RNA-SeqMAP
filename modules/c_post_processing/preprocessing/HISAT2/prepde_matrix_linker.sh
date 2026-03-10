#!/bin/bash

# ===============================================
# prepDE Matrix Linker - M1 HISAT2 Ref-Guided
# ===============================================
# Description:
# During M1 alignment, prepDE.py already produced integer count matrices:
#   gene_count_matrix.csv, transcript_count_matrix.csv, sample_metadata.csv
# under: 2_ALIGNMENT_RESULTs/M1_HISAT2_RefGuided/stringtie_WD/<fasta_tag>/deseq2_input/
#
# This script validates those matrices exist and copies them into the canonical
# post-processing location so downstream DESeq2 analysis modules can consume them:
#   3_POST_PROC/M1_HISAT2_RefGuided/count_matrices_from_stringtie/<fasta_tag>/deseq2_input/
# ===============================================

set -euo pipefail

# ===============================================
# CONFIGURATION
# ===============================================

BASE_DIR="${BASE_DIR:-$PWD}"
MASTER_REFERENCE="${MASTER_REFERENCE:-All_Smel_Genes}"

# Source location (produced by prepDE.py during alignment)
SOURCE_DESEQ2_DIR="${BASE_DIR}/2_ALIGNMENT_RESULTs/M1_HISAT2_RefGuided/stringtie_WD/${MASTER_REFERENCE}/deseq2_input"

# Target location (post-processing canonical path)
TARGET_DESEQ2_DIR="${BASE_DIR}/3_POST_PROC/M1_HISAT2_RefGuided/count_matrices_from_stringtie/${MASTER_REFERENCE}/deseq2_input"

ts() { date '+%Y-%m-%d %H:%M:%S'; }

echo "[$(ts)] prepDE Matrix Linker - M1 HISAT2 Ref-Guided"
echo "[$(ts)] MASTER_REFERENCE: $MASTER_REFERENCE"
echo "[$(ts)] Source: $SOURCE_DESEQ2_DIR"
echo "[$(ts)] Target: $TARGET_DESEQ2_DIR"

# ===============================================
# VALIDATE SOURCE
# ===============================================

if [[ ! -d "$SOURCE_DESEQ2_DIR" ]]; then
    echo "[$(ts)] ERROR: Source deseq2_input directory not found: $SOURCE_DESEQ2_DIR" >&2
    echo "[$(ts)] Run the M1 HISAT2 alignment pipeline first (it calls prepDE.py)." >&2
    exit 1
fi

GENE_MATRIX="$SOURCE_DESEQ2_DIR/gene_count_matrix.csv"
if [[ ! -f "$GENE_MATRIX" ]]; then
    echo "[$(ts)] ERROR: gene_count_matrix.csv not found: $GENE_MATRIX" >&2
    echo "[$(ts)] prepDE.py may not have completed successfully during alignment." >&2
    exit 1
fi

# Basic validation: check it has content (header + at least 1 gene row)
local_rows=$(tail -n +2 "$GENE_MATRIX" | wc -l)
if [[ "$local_rows" -lt 1 ]]; then
    echo "[$(ts)] ERROR: gene_count_matrix.csv appears empty (0 gene rows)" >&2
    exit 1
fi
local_samples=$(head -n1 "$GENE_MATRIX" | tr ',' '\n' | tail -n +2 | wc -l)
echo "[$(ts)] Validated: $local_rows genes, $local_samples samples in gene_count_matrix.csv"

# ===============================================
# STAGE TO POST-PROCESSING TARGET
# ===============================================

mkdir -p "$TARGET_DESEQ2_DIR"

for fname in gene_count_matrix.csv transcript_count_matrix.csv sample_metadata.csv; do
    src="$SOURCE_DESEQ2_DIR/$fname"
    dst="$TARGET_DESEQ2_DIR/$fname"
    if [[ ! -f "$src" ]]; then
        echo "[$(ts)] Warning: $fname not found in source, skipping"
    elif [[ -f "$dst" && "${OVERWRITE_EXISTING:-FALSE}" != "TRUE" ]]; then
        echo "[$(ts)] Already staged (skip): $fname"
    else
        cp "$src" "$dst"
        echo "[$(ts)] Staged: $fname"
    fi
done

echo "[$(ts)] M1 DESeq2 matrices ready at: $TARGET_DESEQ2_DIR"

# ===============================================
# BUILD TPM/FPKM/COVERAGE MATRICES FOR HEATMAPS
# ===============================================
# M1 abundance files (_ref_guided_gene_abundances.tsv) carry TPM/FPKM/Coverage
# values needed by heatmap visualization modules.
# Invoke the dedicated M1 matrix builder (separate from M2's stringtie_matrix_builder.sh).

M1_MATRIX_BUILDER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/m1_ref_guided_matrix_builder.sh"
if [[ ! -f "$M1_MATRIX_BUILDER" ]]; then
    echo "[$(ts)] Warning: m1_ref_guided_matrix_builder.sh not found at: $M1_MATRIX_BUILDER" >&2
    echo "[$(ts)] Heatmap matrices will not be built for M1." >&2
    exit 0
fi

echo "[$(ts)] Building M1 TPM/FPKM/Coverage matrices for visualization..."

BASE_DIR="$BASE_DIR" \
MASTER_REFERENCE="$MASTER_REFERENCE" \
GENE_GROUPS_STR="${GENE_GROUPS_STR:-}" \
SRR_COMBINED_LIST_STR="${SRR_COMBINED_LIST_STR:-}" \
SRR_CSV_DIR="${SRR_CSV_DIR:-}" \
CURRENT_DATASET="${CURRENT_DATASET:-}" \
bash "$M1_MATRIX_BUILDER"

echo "[$(ts)] M1 visualization matrices ready at: ${BASE_DIR}/3_POST_PROC/M1_HISAT2_RefGuided/count_matrices_from_stringtie"
