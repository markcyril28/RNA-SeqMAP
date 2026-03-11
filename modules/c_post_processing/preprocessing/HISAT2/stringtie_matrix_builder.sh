#!/bin/bash

# ===============================================
# StringTie Matrix Builder - M2 HISAT2 De Novo
# ===============================================
# Description:
# Generates TPM/FPKM/Coverage matrices from M2 HISAT2 De Novo StringTie
# abundance files (_gene_abundances_de_novo.tsv) for heatmap visualization.
# Process:
# 1. For each gene group, locate de novo abundance files for all samples
# 2. Extract gene names from reference CSV files (centralized in gene_groups_csv/)
# 3. Build matrices (coverage, FPKM, TPM) with genes as rows, samples/organs as columns
# 4. Output matrices to: 3_POST_PROC/M2_HISAT2_DeNovo/count_matrices_from_stringtie/
#
# M1 HISAT2 Ref-Guided: use m1_ref_guided_matrix_builder.sh instead.
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

# Directories
BASE_DIR="${BASE_DIR:-$PWD}"
INPUTS_DIR="${INPUTS_DIR:-${BASE_DIR}/2_ALIGNMENT_RESULTs/M2_HISAT2_DeNovo/stringtie_WD}"

# M2 de novo abundance filename suffix — hardcoded, not configurable
# For M1 ref-guided, use m1_ref_guided_matrix_builder.sh
ABUNDANCE_SUFFIX="_gene_abundances_de_novo.tsv"

OUT_DIR="${OUT_DIR:-${BASE_DIR}/3_POST_PROC/M2_HISAT2_DeNovo/count_matrices_from_stringtie}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Utilities directory (contains matrix_builder.py)
UTILITIES_DIR="${UTILITIES_DIR:-$SCRIPT_DIR/../../utilities}"

# Create output directory and logging
mkdir -p "$OUT_DIR/logs"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting StringTie Matrix Builder"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working directory: $BASE_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input directory: $INPUTS_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUT_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Master reference: $MASTER_REFERENCE"

# ===============================================
# LOAD SAMPLE IDS FROM CSV (DRY principle)
# ===============================================
# Load from SRR_csv directory instead of hardcoding
# SRR_CSV_DIR is exported by run_all_post_processing.sh
# Fallback to inputs/SRR_csv relative to the project root

SRR_CSV_DIR="${SRR_CSV_DIR:-$SCRIPT_DIR/../../../inputs/SRR_csv}"

# Parse CSV and build arrays (reuse function from pipeline_utils.sh if available)
load_samples_from_csv() {
    local csv_dir="$1"
    local -n sample_ids_ref=$2
    local -n srr_to_organ_ref=$3

    if [[ ! -d "$csv_dir" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: SRR_csv directory not found: $csv_dir"
        return 1
    fi

    # O(1) dedup via associative array (replaces O(n) string-scan per entry)
    local -A _seen=()

    # Sort CSV files for deterministic processing order across filesystems
    while IFS= read -r csv_file; do
        [[ ! -f "$csv_file" ]] && continue
        while IFS=',' read -r srr_id organ notes || [[ -n "$srr_id" ]]; do
            # Skip header and comments
            [[ "$srr_id" =~ ^#.*$ || "$srr_id" == "SRR_ID" || -z "$srr_id" ]] && continue
            srr_id="${srr_id//[[:space:]]/}"
            organ="${organ//[[:space:]]/}"

            # Add if not already present (O(1) hash lookup)
            if [[ -z "${_seen[$srr_id]+x}" ]]; then
                _seen["$srr_id"]=1
                sample_ids_ref+=("$srr_id")
                srr_to_organ_ref["$srr_id"]="$organ"
            fi
        done < "$csv_file"
    done < <(find "$csv_dir" -maxdepth 1 -name "*.csv" -type f | sort)
}

# Initialize arrays
SAMPLE_IDS=()
declare -A SRR_TO_ORGAN=()

# Load from CSV files (DRY - single source of truth)
load_samples_from_csv "$SRR_CSV_DIR" SAMPLE_IDS SRR_TO_ORGAN || true

if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No samples loaded from CSV files"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expected CSV files in: $SRR_CSV_DIR"
    exit 1
fi

# Filter to only configured samples (from SRR_COMBINED_LIST_STR environment variable)
# IMPORTANT: Use the order from SRR_COMBINED_LIST_STR to preserve CSV file order
if [[ -n "${SRR_COMBINED_LIST_STR:-}" ]]; then
    # Parse SRR_COMBINED_LIST_STR (space-separated "SRR_ID:Organ" pairs)
    declare -a CONFIGURED_SRRS=()
    for entry in $SRR_COMBINED_LIST_STR; do
        srr_id="${entry%%:*}"  # Extract SRR ID before colon
        CONFIGURED_SRRS+=("$srr_id")
    done
    
    # Use CONFIGURED_SRRS order (preserves CSV file order)
    # Only include samples that exist in SAMPLE_IDS (have data files)
    # Build O(1) lookup set from SAMPLE_IDS
    declare -A _sample_set=()
    for srr in "${SAMPLE_IDS[@]}"; do _sample_set["$srr"]=1; done
    declare -a FILTERED_SAMPLE_IDS=()
    for srr in "${CONFIGURED_SRRS[@]}"; do
        if [[ -n "${_sample_set[$srr]+x}" ]]; then
            FILTERED_SAMPLE_IDS+=("$srr")
        fi
    done
    SAMPLE_IDS=("${FILTERED_SAMPLE_IDS[@]}")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filtered to ${#SAMPLE_IDS[@]} configured samples"

    if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No configured samples found in CSV data"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] SRR_COMBINED_LIST_STR entries did not match any SRR_IDs in: $SRR_CSV_DIR"
        exit 1
    fi
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using ${#SAMPLE_IDS[@]} samples"

# StringTie -A abundance file column indices (1-based)
# Columns: 1=Gene_ID, 2=Gene_Name, 3=Reference, 4=Strand, 5=Start, 6=End, 7=Coverage, 8=FPKM, 9=TPM
# For M2 de novo (transcript FASTA alignment), column 3 (Reference) contains the
# transcript ID from the FASTA header, which we use to match against gene group CSVs.
GENENAME_COL=3      # Reference column (transcript/contig ID from FASTA)
COVERAGE_COL=7      # Coverage values (per-base depth, NOT raw fragment counts)
FPKM_COL=8          # FPKM values
TPM_COL=9           # TPM values

# ===============================================
# FUNCTIONS
# ===============================================

# Generate combined output folder name: GeneGroup_in_Dataset
get_output_folder_name() {
    local gene_group="$1"
    local dataset="${CURRENT_DATASET:-}"
    if [[ -n "$dataset" ]]; then
        echo "${gene_group}_in_${dataset}"
    else
        echo "$gene_group"
    fi
}

# Function: merge_group_counts
# Purpose: Process abundance files for a gene group and create count matrices
merge_group_counts() {
    local gene_group="$1"
    local ref_csv="$2"
    local group_name
    group_name=$(get_output_folder_name "$gene_group")

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing gene group: $gene_group -> Output: $group_name"
    
    mkdir -p "$OUT_DIR/$group_name"

    local tmpdir
    tmpdir=$(mktemp -d)
    # NOTE: Do NOT use 'trap ... RETURN' here. This function is called from
    # build_full_transcriptome_matrix(), and in bash nested RETURN traps
    # replace each other — the inner trap would clobber the outer, causing
    # unbound-variable errors under set -u. Use explicit rm at each exit.

    # Collect abundance files and build SRR→file map (O(n) instead of O(n²) nested loops)
    local -A srr_to_file=()
    local -a processed_srrs=()
    local files_found=0
    for srr in "${SAMPLE_IDS[@]}"; do
        local file_path="$INPUTS_DIR/$MASTER_REFERENCE/$srr/${srr}_${MASTER_REFERENCE}${ABUNDANCE_SUFFIX}"
        if [[ -f "$file_path" ]]; then
            srr_to_file["$srr"]="$file_path"
            processed_srrs+=("$srr")
            files_found=$((files_found + 1))
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: File not found: $file_path"
        fi
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found $files_found abundance files"

    if [[ $files_found -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No abundance files found for gene group '$gene_group'"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Hint: MASTER_REFERENCE='$MASTER_REFERENCE' must match the fasta_tag used during M2 alignment"
        rm -rf "$tmpdir"
        return 1
    fi

    # Extract gene names from reference CSV (first column is Gene_ID)
    tail -n +2 "${ref_csv}" | cut -d',' -f1 > "$tmpdir/gene_names.txt" \
        || { echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: Failed to extract gene names from $ref_csv"; rm -rf "$tmpdir"; return 1; }
    local gene_name_count
    gene_name_count=$(grep -c . "$tmpdir/gene_names.txt" || true)
    if [[ "$gene_name_count" -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No genes found in reference CSV: $ref_csv"
        rm -rf "$tmpdir"
        return 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Gene names extracted: $gene_name_count lines."

    for count_type in coverage fpkm tpm; do
        local COUNT_COL_VAR="${count_type^^}_COL"
        local COUNT_COL="${!COUNT_COL_VAR}"

        # Extract count data from each sample file (O(1) lookup via srr_to_file map)
        local sample_files=()

        for srr in "${processed_srrs[@]}"; do
            tail -n +2 "${srr_to_file[$srr]}" | cut -f"$GENENAME_COL","$COUNT_COL" > "$tmpdir/${srr}.txt"
            sample_files+=("$tmpdir/${srr}.txt")
        done

        if [[ ${#sample_files[@]} -eq 0 ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: No sample files for $count_type in $gene_group, skipping matrix"
            continue
        fi

        local output_geneName_SRR_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_SRR${MASTER_SUFFIX}.tsv"
        local output_geneName_Organ_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_Organ${MASTER_SUFFIX}.tsv"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating SRR + Organ matrices: $(basename "$output_geneName_SRR_tsv")"

        printf "%s\n" "${sample_files[@]}" > "$tmpdir/sample_files_list.txt"

        # Call matrix_builder.py ONCE; reuse body for both SRR and Organ header variants
        local matrix_body="$tmpdir/matrix_body_${count_type}.txt"
        python3 "$UTILITIES_DIR/matrix_builder.py" "$tmpdir/gene_names.txt" "$tmpdir/sample_files_list.txt" \
            > "$matrix_body" \
            || { echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: matrix_builder.py failed for $group_name"; rm -rf "$tmpdir"; return 1; }

        # SRR header + body
        {
            printf "GeneName"
            for srr in "${processed_srrs[@]}"; do printf "\t%s" "$srr"; done
            printf "\n"
            cat "$matrix_body"
        } > "$output_geneName_SRR_tsv"

        # Organ header + body
        {
            printf "GeneName"
            for srr in "${processed_srrs[@]}"; do printf "\t%s" "${SRR_TO_ORGAN[$srr]:-Unknown}"; done
            printf "\n"
            cat "$matrix_body"
        } > "$output_geneName_Organ_tsv"

        rm -f "$matrix_body"

        rm -f "${sample_files[@]}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed $count_type matrix generation"
    done
    
    rm -rf "$tmpdir"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed processing for $group_name"
}

# ===============================================
# FULL TRANSCRIPTOME MATRIX
# ===============================================
# Build a full-transcriptome matrix (all genes from abundance files) so that
# WGCNA and genome-wide analyses can consume M2 data.
# Uses MASTER_REFERENCE as the gene group name to match build_input_path() conventions.

build_full_transcriptome_matrix() {
    local group_name
    group_name=$(get_output_folder_name "$MASTER_REFERENCE")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building full-transcriptome matrix: $group_name"

    # Collect the UNION of gene IDs from ALL sample abundance files.
    # In de novo mode, StringTie omits zero-coverage transcripts, so any single
    # file may be missing genes that are expressed in other samples.
    local tmp_csv
    tmp_csv=$(mktemp --suffix=.csv)
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
        rm -f "$tmp_csv"
        return 1
    fi

    # De-duplicate the collected gene IDs (header line is already first)
    local tmp_dedup
    tmp_dedup=$(mktemp --suffix=.csv)
    head -1 "$tmp_csv" > "$tmp_dedup"
    tail -n +2 "$tmp_csv" | sort -u >> "$tmp_dedup"
    mv "$tmp_dedup" "$tmp_csv"

    local gene_count
    gene_count=$(tail -n +2 "$tmp_csv" | grep -c . || true)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Full transcriptome: $gene_count genes (union from $files_found samples)"

    # Reuse merge_group_counts with MASTER_REFERENCE as the gene group name
    if merge_group_counts "$MASTER_REFERENCE" "$tmp_csv"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Full-transcriptome matrix complete"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed to build full-transcriptome matrix"
    fi
    rm -f "$tmp_csv"
}

# ===============================================
# MAIN EXECUTION
# ===============================================

# Centralized gene groups CSV directory
# Use GENE_GROUPS_DIR from environment (set by run_all_post_processing.sh)
# Fallback to inputs/gene_groups_csv relative to the project root
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
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found reference CSV with $(tail -n +2 "$REF_CSV" | grep -c . || true) genes"
    
    if merge_group_counts "$gene_group" "$REF_CSV"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Successfully processed $gene_group"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed to process $gene_group"
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Matrix generation completed"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUT_DIR"
