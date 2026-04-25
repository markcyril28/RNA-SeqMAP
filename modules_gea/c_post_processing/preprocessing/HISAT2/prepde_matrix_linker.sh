#!/bin/bash

# ===============================================
# prepDE Matrix Linker - M1 HISAT2 Ref-Guided
# ===============================================
# Description:
# During M1 alignment, prepDE.py already produced integer count matrices:
#   gene_count_matrix.csv, transcript_count_matrix.csv, sample_metadata.csv
# under: II_RESULTS/2_ALIGNMENT_RESULTs/M1_HISAT2_RefGuided/stringtie_WD/<fasta_tag>/deseq2_input/
#
# This script validates those matrices exist and copies them into the canonical
# post-processing location so downstream DESeq2 analysis modules can consume them:
#   II_RESULTS/3_POST_PROC/M1_HISAT2_RefGuided/count_matrices_from_stringtie/<fasta_tag>/deseq2_input/
# ===============================================

set -euo pipefail

# ===============================================
# CONFIGURATION
# ===============================================

# Fail early under orchestrators if BASE_DIR/PROJECT_ROOT are missing
if [[ -n "${WF_MANAGED_ENV:-}" && -z "${BASE_DIR:-}" && -z "${PROJECT_ROOT:-}" ]]; then
    echo "ERROR: WF_MANAGED_ENV is set but neither BASE_DIR nor PROJECT_ROOT is set. Orchestrators must export BASE_DIR." >&2
    exit 1
fi

BASE_DIR="${BASE_DIR:-${PROJECT_ROOT:-$PWD}}"
MASTER_REFERENCE="${MASTER_REFERENCE:-All_Smel_Genes}"

# Source logging utilities for consistent pipeline logging
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" 2>/dev/null && pwd)" || SCRIPT_DIR="."
source "${BASE_DIR}/modules_gea/logging/logging_utils.sh" 2>/dev/null || {
    log_info()  { echo "[INFO] $*"; }
    log_warn()  { echo "[WARN] $*" >&2; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_step()  { echo "=== $* ==="; }
}

# Source location (produced by prepDE.py during alignment)
SOURCE_DESEQ2_DIR="${BASE_DIR}/II_RESULTS/2_ALIGNMENT_RESULTs/M1_HISAT2_RefGuided/stringtie_WD/${MASTER_REFERENCE}/deseq2_input"

# Target location (post-processing canonical path)
TARGET_DESEQ2_DIR="${BASE_DIR}/II_RESULTS/3_POST_PROC/${CURRENT_GENE_GROUP:-_active}/M1_HISAT2_RefGuided/count_matrices_from_stringtie/${MASTER_REFERENCE}/deseq2_input"

log_info "prepDE Matrix Linker - M1 HISAT2 Ref-Guided"
log_info "MASTER_REFERENCE: $MASTER_REFERENCE"
log_info "Source: $SOURCE_DESEQ2_DIR"
log_info "Target: $TARGET_DESEQ2_DIR"

# ===============================================
# VALIDATE SOURCE
# ===============================================

if [[ ! -d "$SOURCE_DESEQ2_DIR" ]]; then
    log_error "Source deseq2_input directory not found: $SOURCE_DESEQ2_DIR"
    log_error "Run the M1 HISAT2 alignment pipeline first (it calls prepDE.py)."
    exit 1
fi

GENE_MATRIX="$SOURCE_DESEQ2_DIR/gene_count_matrix.csv"
if [[ ! -f "$GENE_MATRIX" ]]; then
    log_error "gene_count_matrix.csv not found: $GENE_MATRIX"
    log_error "prepDE.py may not have completed successfully during alignment."
    exit 1
fi

# Validation: single-pass awk extracts all metrics (replaces 4 separate file reads)
read -r local_rows local_samples first_gene first_count < <(awk -F',' '
    NR == 1 { samples = NF - 1 }
    NR == 2 { gene = $1; count = $2 }
    NR > 1  { rows++ }
    END     { print rows+0, samples+0, gene, count }
' "$GENE_MATRIX")

if [[ "$local_rows" -lt 1 ]]; then
    log_error "gene_count_matrix.csv appears empty (0 gene rows)"
    exit 1
fi
if [[ "$local_samples" -lt 1 ]]; then
    log_error "gene_count_matrix.csv has no sample columns (only gene ID column found)"
    exit 1
fi
if [[ -z "$first_gene" || "$first_gene" == "," ]]; then
    log_error "gene_count_matrix.csv has empty gene names in first column"
    exit 1
fi
if ! [[ "$first_count" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    log_error "gene_count_matrix.csv contains non-numeric count values (got: '$first_count')"
    exit 1
fi

log_info "Validated: $local_rows genes, $local_samples samples in gene_count_matrix.csv"

# Warn if the count matrix is missing samples relative to the configured dataset
if [[ -n "${SRR_COMBINED_LIST_STR:-}" ]]; then
    # Pure bash word count — avoids echo|wc subprocess spawn. O(1).
    read -ra _tmp_arr <<< "$SRR_COMBINED_LIST_STR"
    configured_samples=${#_tmp_arr[@]}
    if [[ "$local_samples" -lt "$configured_samples" ]]; then
        log_warn "gene_count_matrix.csv has $local_samples samples but $configured_samples are configured."
        log_warn "DESeq2 Differential_Expression will only cover the $local_samples aligned samples."
        log_warn "Run M1 alignment for all configured samples to include them in the count matrix."
    fi
fi

# ===============================================
# STAGE TO POST-PROCESSING TARGET
# ===============================================

mkdir -p "$TARGET_DESEQ2_DIR" || { log_error "Failed to create target directory: $TARGET_DESEQ2_DIR"; exit 1; }

for fname in gene_count_matrix.csv transcript_count_matrix.csv sample_metadata.csv; do
    src="$SOURCE_DESEQ2_DIR/$fname"
    dst="$TARGET_DESEQ2_DIR/$fname"
    if [[ ! -f "$src" ]]; then
        log_warn "$fname not found in source, skipping"
    elif [[ -f "$dst" && "${OVERWRITE_EXISTING:-FALSE}" != "TRUE" ]]; then
        log_info "Already staged (skip): $fname"
    else
        cp "$src" "$dst"
        log_info "Staged: $fname"
    fi
done

log_info "M1 DESeq2 matrices ready at: $TARGET_DESEQ2_DIR"

# ===============================================
# BUILD TPM/FPKM/COVERAGE MATRICES FOR HEATMAPS
# ===============================================
# M1 abundance files (_ref_guided_gene_abundances.tsv) carry TPM/FPKM/Coverage
# values needed by heatmap visualization modules.
# Invoke the dedicated M1 matrix builder (separate from M2's stringtie_matrix_builder.sh).

M1_MATRIX_BUILDER="${BASH_SOURCE[0]%/*}/m1_ref_guided_matrix_builder.sh"
if [[ ! -f "$M1_MATRIX_BUILDER" ]]; then
    log_warn "m1_ref_guided_matrix_builder.sh not found at: $M1_MATRIX_BUILDER"
    log_warn "Heatmap matrices will not be built for M1."
    exit 0
fi

log_info "Building M1 TPM/FPKM/Coverage matrices for visualization..."

BASE_DIR="$BASE_DIR" \
MASTER_REFERENCE="$MASTER_REFERENCE" \
STRINGTIE_METHOD="M1" \
CURRENT_GENE_GROUP="${CURRENT_GENE_GROUP:-_active}" \
GENE_GROUPS_STR="${GENE_GROUPS_STR:-}" \
GENE_GROUPS_DIR="${GENE_GROUPS_DIR:-}" \
SRR_COMBINED_LIST_STR="${SRR_COMBINED_LIST_STR:-}" \
SRR_CSV_DIR="${SRR_CSV_DIR:-}" \
CURRENT_DATASET="${CURRENT_DATASET:-}" \
bash "$M1_MATRIX_BUILDER"

log_info "M1 visualization matrices ready at: ${BASE_DIR}/II_RESULTS/3_POST_PROC/${CURRENT_GENE_GROUP:-_active}/M1_HISAT2_RefGuided/count_matrices_from_stringtie"
