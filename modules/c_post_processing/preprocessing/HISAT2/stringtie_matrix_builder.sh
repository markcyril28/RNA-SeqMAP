#!/bin/bash

# ===============================================
# StringTie Matrix Builder - Generic Version
# ===============================================
# Description:
# Generates gene expression count matrices from StringTie abundance files.
# Process:
# 1. For each gene group, locate abundance files for all samples
# 2. Extract gene names from reference CSV files (centralized in gene_groups/)
# 3. Build matrices (coverage, FPKM, TPM) with genes as rows, samples/organs as columns
# 4. Output matrices to results directory
# ===============================================

set -uo pipefail

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
        "SmelGIFs"
    )
fi

# Master reference
MASTER_REFERENCE="${MASTER_REFERENCE:-All_Smel_Genes}"
MASTER_SUFFIX="_from_${MASTER_REFERENCE}"

# Directories
BASE_DIR="${BASE_DIR:-$PWD}"
INPUTS_DIR="${INPUTS_DIR:-stringtie_WD/a_Method_2_RAW_RESULTs}"
OUT_DIR="${OUT_DIR:-stringtie_WD/b_Method_2_COUNT_MATRICES}"
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
# Fallback to 0_INPUTs/SRR_csv relative to the project root

SRR_CSV_DIR="${SRR_CSV_DIR:-$SCRIPT_DIR/../../../0_INPUTs/SRR_csv}"

# Parse CSV and build arrays (reuse function from pipeline_utils.sh if available)
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
            # Skip header and comments
            [[ "$srr_id" =~ ^#.*$ || "$srr_id" == "SRR_ID" || -z "$srr_id" ]] && continue
            srr_id=$(echo "$srr_id" | tr -d '[:space:]')
            organ=$(echo "$organ" | tr -d '[:space:]')
            
            # Add if not already present
            if [[ ! " ${sample_ids_ref[*]} " =~ " ${srr_id} " ]]; then
                sample_ids_ref+=("$srr_id")
                srr_to_organ_ref["$srr_id"]="$organ"
            fi
        done < "$csv_file"
    done
}

# Initialize arrays
SAMPLE_IDS=()
declare -A SRR_TO_ORGAN=()

# Load from CSV files (DRY - single source of truth)
load_samples_from_csv "$SRR_CSV_DIR" SAMPLE_IDS SRR_TO_ORGAN

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
GENENAME_COL=3      # Gene name column
COVERAGE_COL=7      # Coverage values
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
    local gene_group_path="$2"
    local ref_csv="$3"
    local group_name=$(get_output_folder_name "$gene_group")

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing gene group: $gene_group -> Output: $group_name"
    
    mkdir -p "$OUT_DIR/$group_name"

    local tmpdir
    tmpdir=$(mktemp -d)

    # Collect abundance files
    files=()
    local files_found=0
    for srr in "${SAMPLE_IDS[@]}"; do
        file_path="$INPUTS_DIR/$MASTER_REFERENCE/$srr/${srr}_${MASTER_REFERENCE}_gene_abundances_de_novo.tsv"
        if [[ -f "$file_path" ]]; then
            files+=("$file_path")
            ((files_found++))
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: File not found: $file_path"
        fi
    done
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found $files_found abundance files"
    
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No abundance files found"
        rm -r "$tmpdir"
        return 1
    fi

    # Extract gene names from reference CSV (first column is Gene_ID)
    tail -n +2 "${ref_csv}" | cut -d',' -f1 > "$tmpdir/gene_names.txt"
    echo "Gene names extracted: $(wc -l < "$tmpdir/gene_names.txt") lines."

    for count_type in coverage fpkm tpm; do
        local COUNT_COL_VAR="${count_type^^}_COL"
        local COUNT_COL="${!COUNT_COL_VAR}"

        # Extract count data from each sample file
        sample_files=()
        
        for srr in "${SAMPLE_IDS[@]}"; do
            sample_file=""
            for f in "${files[@]}"; do
                if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                    sample_file="$f"
                    break
                fi
            done
            
            if [[ -n "$sample_file" ]]; then
                tail -n +2 "$sample_file" | cut -f"$GENENAME_COL","$COUNT_COL" > "$tmpdir/${srr}.txt"
                sample_files+=("$tmpdir/${srr}.txt")
            fi
        done

        local output_geneName_SRR_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_SRR${MASTER_SUFFIX}.tsv"
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating SRR matrix: $(basename "$output_geneName_SRR_tsv")"
        
        printf "%s\n" "${sample_files[@]}" > "$tmpdir/sample_files_list.txt"
        
        # Create matrix with SRR headers
        {
            printf "GeneName"
            for srr in "${SAMPLE_IDS[@]}"; do
                for f in "${files[@]}"; do
                    if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                        printf "\t%s" "$srr"
                        break
                    fi
                done
            done
            printf "\n"

            python3 "$UTILITIES_DIR/matrix_builder.py" "$tmpdir/gene_names.txt" "$tmpdir/sample_files_list.txt"
        } > "$output_geneName_SRR_tsv"
        
        # Create matrix with Organ headers
        local output_geneName_Organ_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_Organ${MASTER_SUFFIX}.tsv"
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating Organ matrix: $(basename "$output_geneName_Organ_tsv")"
        
        {
            printf "GeneName"
            for srr in "${SAMPLE_IDS[@]}"; do
                for f in "${files[@]}"; do
                    if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                        local organ="${SRR_TO_ORGAN[$srr]:-Unknown}"
                        printf "\t%s" "$organ"
                        break
                    fi
                done
            done
            printf "\n"

            python3 "$UTILITIES_DIR/matrix_builder.py" "$tmpdir/gene_names.txt" "$tmpdir/sample_files_list.txt"
        } > "$output_geneName_Organ_tsv"

        rm -f "${sample_files[@]}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed $count_type matrix generation"
    done
    
    rm -r "$tmpdir"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed processing for $group_name"
}

# ===============================================
# MAIN EXECUTION
# ===============================================
 
# Centralized gene groups CSV directory
# Use GENE_GROUPS_DIR from environment (set by run_all_post_processing.sh)
# Fallback to 0_INPUTs/gene_groups relative to the project root
GENE_GROUPS_CSV_DIR="${GENE_GROUPS_DIR:-${GENE_GROUPS_CSV_DIR:-$SCRIPT_DIR/../../../0_INPUTs/gene_groups}}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Gene groups CSV directory: $GENE_GROUPS_CSV_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting count matrix generation for ${#GENE_GROUPS[@]} gene groups"

for gene_group in "${GENE_GROUPS[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing gene group: $gene_group"
    
    REF_CSV="${GENE_GROUPS_CSV_DIR}/${gene_group}.csv"
    
    if [[ ! -f "$REF_CSV" ]]; then 
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: Reference CSV not found: $REF_CSV, skipping $gene_group"
        continue
    fi
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found reference CSV with $(tail -n +2 "$REF_CSV" | wc -l) genes"
    
    gene_group_path="$INPUTS_DIR/$gene_group"
    
    if merge_group_counts "$gene_group" "$gene_group_path" "$REF_CSV"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Successfully processed $gene_group"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed to process $gene_group"
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Matrix generation completed"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUT_DIR"
