#!/bin/bash

# ===============================================
# M1 HISAT2 Ref-Guided Matrix Builder
# ===============================================
# Description:
# Generates TPM/FPKM/Coverage matrices from HISAT2 Ref-Guided StringTie
# abundance files (_ref_guided_gene_abundances.tsv) for heatmap visualization.
# Process:
# 1. For each gene group, locate ref-guided abundance files for all samples
# 2. Extract gene names from reference CSV files (centralized in gene_groups_csv/)
# 3. Build matrices (coverage, FPKM, TPM) with genes as rows, samples/organs as columns
# 4. Output matrices to: 3_POST_PROC/M1_HISAT2_RefGuided/count_matrices_from_stringtie/
#
# Called by: prepde_matrix_linker.sh (as part of M1 preprocessing)
# Input:     2_ALIGNMENT_RESULTs/M1_HISAT2_RefGuided/stringtie_WD/
# Abundance: ${SRR}_${MASTER_REFERENCE}_ref_guided_gene_abundances.tsv
# ===============================================

set -euo pipefail

# ===============================================
# CONFIGURATION
# ===============================================

# Gene groups to process (override from environment or use defaults)
# Reads from GENE_GROUPS_STR for parallel compatibility
if [[ -n "${GENE_GROUPS_STR:-}" ]]; then
    IFS=' ' read -ra GENE_GROUPS <<< "$GENE_GROUPS_STR"
elif [[ -z "${GENE_GROUPS:-}" ]]; then
    GENE_GROUPS=(
        "SmelDMPs"
        "SmelGRFs"
        "SmelGIF"
    )
fi

# Master reference
MASTER_REFERENCE="${MASTER_REFERENCE:-All_Smel_Genes}"
MASTER_SUFFIX="_from_${MASTER_REFERENCE}"

# Directories — M1 specific, not shared with M2
BASE_DIR="${BASE_DIR:-$PWD}"
INPUTS_DIR="${INPUTS_DIR:-${BASE_DIR}/2_ALIGNMENT_RESULTs/M1_HISAT2_RefGuided/stringtie_WD}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/3_POST_PROC/M1_HISAT2_RefGuided/count_matrices_from_stringtie}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# M1 ref-guided abundance filename suffix (single-pass quantification with -e)
ABUNDANCE_SUFFIX="_ref_guided_gene_abundances.tsv"

# Utilities directory (contains matrix_builder.py)
UTILITIES_DIR="${UTILITIES_DIR:-$SCRIPT_DIR/../../utilities}"

# Create output directory and logging
mkdir -p "$OUT_DIR/logs"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting M1 Ref-Guided Matrix Builder"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working directory: $BASE_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input directory: $INPUTS_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUT_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Master reference: $MASTER_REFERENCE"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Abundance suffix: $ABUNDANCE_SUFFIX"

# ===============================================
# LOAD SAMPLE IDS FROM CSV
# ===============================================
# SRR_CSV_DIR is exported by run_all_post_processing.sh
# Fallback to inputs/SRR_csv relative to the project root

SRR_CSV_DIR="${SRR_CSV_DIR:-$SCRIPT_DIR/../../../inputs/SRR_csv}"

load_samples_from_csv() {
    local csv_dir="$1"
    local -n sample_ids_ref=$2
    local -n srr_to_organ_ref=$3

    if [[ ! -d "$csv_dir" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: SRR_csv directory not found: $csv_dir"
        return 1
    fi

    for csv_file in "$csv_dir"/*.csv; do
        [[ ! -f "$csv_file" ]] && continue
        while IFS=',' read -r srr_id organ notes || [[ -n "$srr_id" ]]; do
            [[ "$srr_id" =~ ^#.*$ || "$srr_id" == "SRR_ID" || -z "$srr_id" ]] && continue
            srr_id=$(echo "$srr_id" | tr -d '[:space:]')
            organ=$(echo "$organ" | tr -d '[:space:]')
            if [[ ! " ${sample_ids_ref[*]} " =~ " ${srr_id} " ]]; then
                sample_ids_ref+=("$srr_id")
                srr_to_organ_ref["$srr_id"]="$organ"
            fi
        done < "$csv_file"
    done
}

SAMPLE_IDS=()
declare -A SRR_TO_ORGAN=()

load_samples_from_csv "$SRR_CSV_DIR" SAMPLE_IDS SRR_TO_ORGAN || true

if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No samples loaded from CSV files"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expected CSV files in: $SRR_CSV_DIR"
    exit 1
fi

# Filter to configured samples in SRR_COMBINED_LIST_STR order
if [[ -n "${SRR_COMBINED_LIST_STR:-}" ]]; then
    declare -a CONFIGURED_SRRS=()
    for entry in $SRR_COMBINED_LIST_STR; do
        srr_id="${entry%%:*}"
        CONFIGURED_SRRS+=("$srr_id")
    done
    declare -a FILTERED_SAMPLE_IDS=()
    for srr in "${CONFIGURED_SRRS[@]}"; do
        if [[ " ${SAMPLE_IDS[*]} " =~ " ${srr} " ]]; then
            FILTERED_SAMPLE_IDS+=("$srr")
        fi
    done
    SAMPLE_IDS=("${FILTERED_SAMPLE_IDS[@]}")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filtered to ${#SAMPLE_IDS[@]} configured samples"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using ${#SAMPLE_IDS[@]} samples"

# StringTie abundance file column indices (1-based)
# Header: Gene ID | Gene Name | Reference | Strand | Start | End | Coverage | FPKM | TPM
GENENAME_COL=1      # Gene ID column (e.g., SMEL5_01g000100)
COVERAGE_COL=7      # Coverage values
FPKM_COL=8          # FPKM values
TPM_COL=9           # TPM values

# ===============================================
# FUNCTIONS
# ===============================================

get_output_folder_name() {
    local gene_group="$1"
    local dataset="${CURRENT_DATASET:-}"
    if [[ -n "$dataset" ]]; then
        echo "${gene_group}_in_${dataset}"
    else
        echo "$gene_group"
    fi
}

merge_group_counts() {
    local gene_group="$1"
    local ref_csv="$2"
    local group_name
    group_name=$(get_output_folder_name "$gene_group")

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing gene group: $gene_group -> Output: $group_name"

    mkdir -p "$OUT_DIR/$group_name"

    local tmpdir
    tmpdir=$(mktemp -d)
    # Clean up tmpdir on function return. Using a subshell-safe cleanup that also
    # handles early returns (e.g., no abundance files found).
    trap 'rm -rf "$tmpdir"; trap - RETURN' RETURN

    # Count available abundance files (paths are constructed directly per-sample)
    local files_found=0
    for srr in "${SAMPLE_IDS[@]}"; do
        local file_path="$INPUTS_DIR/$MASTER_REFERENCE/$srr/${srr}_${MASTER_REFERENCE}${ABUNDANCE_SUFFIX}"
        if [[ -f "$file_path" ]]; then
            (( files_found++ )) || true
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: File not found: $file_path"
        fi
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found $files_found ref-guided abundance files"

    if [[ $files_found -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No abundance files found for $gene_group"
        return 1
    fi

    # Extract gene names from reference CSV (first column is Gene_ID)
    tail -n +2 "${ref_csv}" | cut -d',' -f1 > "$tmpdir/gene_names.txt"
    echo "Gene names extracted: $(grep -c . "$tmpdir/gene_names.txt") lines."

    for count_type in coverage fpkm tpm; do
        local COUNT_COL_VAR="${count_type^^}_COL"
        local COUNT_COL="${!COUNT_COL_VAR}"

        local -a sample_files=()
        local -a processed_srrs=()

        for srr in "${SAMPLE_IDS[@]}"; do
            # Construct path directly — no need to scan files[] array
            local sample_file="$INPUTS_DIR/$MASTER_REFERENCE/$srr/${srr}_${MASTER_REFERENCE}${ABUNDANCE_SUFFIX}"
            if [[ -f "$sample_file" ]]; then
                tail -n +2 "$sample_file" | cut -f"$GENENAME_COL","$COUNT_COL" > "$tmpdir/${srr}.txt"
                sample_files+=("$tmpdir/${srr}.txt")
                processed_srrs+=("$srr")
            fi
        done

        # NOTE: Filename uses "geneName" (camelCase) while the TSV header column is "GeneName" (PascalCase).
        # build_input_path() in 0_shared_config.R maps gene_type=="Shortened_Name" -> "geneName" to match this convention.
        local output_geneName_SRR_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_SRR${MASTER_SUFFIX}.tsv"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating SRR matrix: $(basename "$output_geneName_SRR_tsv")"

        if [[ ${#sample_files[@]} -eq 0 ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: No sample files for $count_type in $gene_group, skipping matrix"
            continue
        fi

        printf "%s\n" "${sample_files[@]}" > "$tmpdir/sample_files_list.txt"

        {
            printf "GeneName"
            for srr in "${processed_srrs[@]}"; do
                printf "\t%s" "$srr"
            done
            printf "\n"
            python3 "$UTILITIES_DIR/matrix_builder.py" "$tmpdir/gene_names.txt" "$tmpdir/sample_files_list.txt"
        } > "$output_geneName_SRR_tsv" || { echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: matrix_builder.py failed for SRR matrix: $group_name"; return 1; }

        local output_geneName_Organ_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_Organ${MASTER_SUFFIX}.tsv"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating Organ matrix: $(basename "$output_geneName_Organ_tsv")"

        {
            printf "GeneName"
            for srr in "${processed_srrs[@]}"; do
                local organ="${SRR_TO_ORGAN[$srr]:-Unknown}"
                printf "\t%s" "$organ"
            done
            printf "\n"
            python3 "$UTILITIES_DIR/matrix_builder.py" "$tmpdir/gene_names.txt" "$tmpdir/sample_files_list.txt"
        } > "$output_geneName_Organ_tsv" || { echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: matrix_builder.py failed for Organ matrix: $group_name"; return 1; }

        rm -f "${sample_files[@]}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed $count_type matrix generation"
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed processing for $group_name"
}

# ===============================================
# FULL TRANSCRIPTOME MATRIX
# ===============================================
# Build a full-transcriptome matrix (all genes from abundance files) so that
# WGCNA and genome-wide analyses can consume M1 data.
# Uses MASTER_REFERENCE as the gene group name to match build_input_path() conventions.

build_full_transcriptome_matrix() {
    local group_name
    group_name=$(get_output_folder_name "$MASTER_REFERENCE")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building full-transcriptome matrix: $group_name"

    # Collect the UNION of gene IDs from ALL sample abundance files.
    # In ref-guided mode all samples share the same gene set, but taking the
    # union keeps this robust if any file is truncated or filtered.
    local tmp_csv
    tmp_csv=$(mktemp --suffix=.csv)
    trap 'rm -f "$tmp_csv"' RETURN
    echo "Gene_ID" > "$tmp_csv"

    local files_found=0
    for srr in "${SAMPLE_IDS[@]}"; do
        local file_path="$INPUTS_DIR/$MASTER_REFERENCE/$srr/${srr}_${MASTER_REFERENCE}${ABUNDANCE_SUFFIX}"
        if [[ -f "$file_path" ]]; then
            tail -n +2 "$file_path" | cut -f"$GENENAME_COL" >> "$tmp_csv"
            files_found=$((files_found + 1))
        fi
    done

    if [[ "$files_found" -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: No abundance files found - skipping full-transcriptome matrix"
        return 1
    fi

    # De-duplicate the collected gene IDs (header line is already first)
    local tmp_dedup
    tmp_dedup=$(mktemp --suffix=.csv)
    head -1 "$tmp_csv" > "$tmp_dedup"
    tail -n +2 "$tmp_csv" | sort -u >> "$tmp_dedup"
    mv "$tmp_dedup" "$tmp_csv"

    local gene_count
    gene_count=$(tail -n +2 "$tmp_csv" | grep -c .)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Full transcriptome: $gene_count genes (union from $files_found samples)"

    # Reuse merge_group_counts with MASTER_REFERENCE as the gene group name
    if merge_group_counts "$MASTER_REFERENCE" "$tmp_csv"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Full-transcriptome matrix complete"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed to build full-transcriptome matrix"
    fi
}

# ===============================================
# MAIN EXECUTION
# ===============================================

GENE_GROUPS_CSV_DIR="${GENE_GROUPS_DIR:-${GENE_GROUPS_CSV_DIR:-$SCRIPT_DIR/../../../inputs/gene_groups_csv}}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Gene groups CSV directory: $GENE_GROUPS_CSV_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting count matrix generation for ${#GENE_GROUPS[@]} gene groups"

# Build full-transcriptome matrix first (enables WGCNA and genome-wide analyses)
build_full_transcriptome_matrix

for gene_group in "${GENE_GROUPS[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing gene group: $gene_group"

    REF_CSV="${GENE_GROUPS_CSV_DIR}/${gene_group}.csv"

    # Search subdirectories if not found at top level
    if [[ ! -f "$REF_CSV" ]]; then
        REF_CSV=$(find "$GENE_GROUPS_CSV_DIR" -maxdepth 3 -name "${gene_group}.csv" -type f -print -quit 2>/dev/null)
    fi

    if [[ -z "$REF_CSV" || ! -f "$REF_CSV" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: Reference CSV not found: ${GENE_GROUPS_CSV_DIR}/${gene_group}.csv, skipping $gene_group"
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found reference CSV with $(tail -n +2 "$REF_CSV" | grep -c .) genes"

    if merge_group_counts "$gene_group" "$REF_CSV"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Successfully processed $gene_group"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed to process $gene_group"
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] M1 Ref-Guided matrix generation completed"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUT_DIR"
