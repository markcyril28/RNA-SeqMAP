#!/bin/bash

# ===============================================
# Unified StringTie Matrix Builder (M1 + M2)
# ===============================================
# Description:
# Generates TPM/FPKM/Coverage matrices from HISAT2 StringTie abundance files
# for heatmap visualization. Handles both M1 (ref-guided) and M2 (de novo).
#
# Process:
# 1. For each gene group, locate abundance files for all samples
# 2. Extract gene names from reference CSV files (centralized in gene_groups_csv/)
# 3. Call matrix_builder.py ONCE per count_type; reuse body for both SRR + Organ headers
# 4. Output matrices to: 3_POST_PROC/{method}/count_matrices_from_stringtie/
#
# Usage:
#   STRINGTIE_METHOD=M1  bash stringtie_matrix_builder.sh   # ref-guided
#   STRINGTIE_METHOD=M2  bash stringtie_matrix_builder.sh   # de novo (default)
#
# Called by: prepde_matrix_linker.sh (M1) or pipeline_utils.sh (M2)
# ===============================================

set -euo pipefail

# Source logging utilities for consistent pipeline logging (minimal fallback if unavailable)
source "${BASE_DIR:-$PWD}/modules/logging/logging_utils.sh" 2>/dev/null || {
    log_info()  { echo "[INFO] $*"; }
    log_warn()  { echo "[WARN] $*" >&2; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_step()  { echo "=== $* ==="; }
}

# ===============================================
# METHOD CONFIGURATION
# ===============================================
# STRINGTIE_METHOD selects M1 (ref-guided) vs M2 (de novo) defaults.
# All derived values can still be overridden individually via environment.

STRINGTIE_METHOD="${STRINGTIE_METHOD:-M2}"

case "$STRINGTIE_METHOD" in
    M1)
        _DEFAULT_INPUTS_SUBDIR="M1_HISAT2_RefGuided/stringtie_WD"
        _DEFAULT_OUT_SUBDIR="M1_HISAT2_RefGuided/count_matrices_from_stringtie"
        _DEFAULT_ABUNDANCE_SUFFIX="_ref_guided_gene_abundances.tsv"
        _DEFAULT_GENENAME_COL=1    # Gene ID column (e.g., SMEL5_01g000100)
        _METHOD_LABEL="M1 Ref-Guided"
        ;;
    M2)
        _DEFAULT_INPUTS_SUBDIR="M2_HISAT2_DeNovo/stringtie_WD"
        _DEFAULT_OUT_SUBDIR="M2_HISAT2_DeNovo/count_matrices_from_stringtie"
        _DEFAULT_ABUNDANCE_SUFFIX="_gene_abundances_de_novo.tsv"
        _DEFAULT_GENENAME_COL=3    # Reference column (transcript/contig ID from FASTA)
        _METHOD_LABEL="M2 De Novo"
        ;;
    *)
        log_error "Unknown STRINGTIE_METHOD='$STRINGTIE_METHOD' (expected M1 or M2)"
        exit 1
        ;;
esac

# ===============================================
# CONFIGURATION (all overridable via environment)
# ===============================================

# Gene groups to process
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
INPUTS_DIR="${INPUTS_DIR:-${BASE_DIR}/2_ALIGNMENT_RESULTs/${_DEFAULT_INPUTS_SUBDIR}}"
OUT_DIR="${OUT_DIR:-${BASE_DIR}/3_POST_PROC/${_DEFAULT_OUT_SUBDIR}}"
# Resolve SCRIPT_DIR without nested dirname subshell
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."

# Abundance filename suffix
ABUNDANCE_SUFFIX="${ABUNDANCE_SUFFIX:-$_DEFAULT_ABUNDANCE_SUFFIX}"

# Gene name column index (1-based) in StringTie abundance file
# Columns: 1=Gene_ID, 2=Gene_Name, 3=Reference, 4=Strand, 5=Start, 6=End, 7=Coverage, 8=FPKM, 9=TPM
GENENAME_COL="${GENENAME_COL:-$_DEFAULT_GENENAME_COL}"

# Utilities directory (contains matrix_builder.py)
UTILITIES_DIR="${UTILITIES_DIR:-$SCRIPT_DIR/../../utilities}"

# Create output directory and logging
mkdir -p "$OUT_DIR/logs"

log_step "Starting StringTie Matrix Builder ($_METHOD_LABEL)"
log_info "Working directory: $BASE_DIR"
log_info "Input directory: $INPUTS_DIR"
log_info "Output directory: $OUT_DIR"
log_info "Master reference: $MASTER_REFERENCE"
log_info "Abundance suffix: $ABUNDANCE_SUFFIX"
log_info "Gene name column: $GENENAME_COL"

# Fixed column indices for count types
COVERAGE_COL=7
FPKM_COL=8
TPM_COL=9

# ===============================================
# LOAD SAMPLE IDS FROM CSV
# ===============================================
# SRR_CSV_DIR is exported by run_post_processing.sh
# Fallback to inputs/3_post_proc_inputs/SRR_csv relative to the project root

SRR_CSV_DIR="${SRR_CSV_DIR:-$SCRIPT_DIR/../../../../inputs/3_post_proc_inputs/SRR_csv}"

load_samples_from_csv() {
    local csv_dir="$1"
    local -n sample_ids_ref=$2
    local -n srr_to_organ_ref=$3

    if [[ ! -d "$csv_dir" ]]; then
        log_warn "SRR_csv directory not found: $csv_dir"
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
    log_error "No samples loaded from CSV files"
    log_info "Expected CSV files in: $SRR_CSV_DIR"
    exit 1
fi

# Filter to only configured samples (from SRR_COMBINED_LIST_STR environment variable)
# IMPORTANT: Use the order from SRR_COMBINED_LIST_STR to preserve CSV file order
if [[ -n "${SRR_COMBINED_LIST_STR:-}" ]]; then
    declare -a CONFIGURED_SRRS=()
    for entry in $SRR_COMBINED_LIST_STR; do
        srr_id="${entry%%:*}"
        CONFIGURED_SRRS+=("$srr_id")
    done
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
    log_info "Filtered to ${#SAMPLE_IDS[@]} configured samples"

    if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
        log_error "No configured samples found in CSV data"
        log_info "SRR_COMBINED_LIST_STR entries did not match any SRR_IDs in: $SRR_CSV_DIR"
        exit 1
    fi
fi

log_info "Using ${#SAMPLE_IDS[@]} samples"

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
# Calls matrix_builder.py ONCE per count_type; reuses body for both SRR and Organ headers
merge_group_counts() {
    local gene_group="$1"
    local ref_csv="$2"
    local group_name
    group_name=$(get_output_folder_name "$gene_group")

    log_info "Processing gene group: $gene_group -> Output: $group_name"

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
            log_warn "File not found: $file_path"
        fi
    done

    log_info "Found $files_found/${#SAMPLE_IDS[@]} abundance files"

    if [[ $files_found -eq 0 ]]; then
        log_error "No abundance files found for gene group '$gene_group'"
        log_info "Hint: MASTER_REFERENCE='$MASTER_REFERENCE' must match the fasta_tag used during alignment"
        rm -rf "$tmpdir"
        return 1
    fi

    if [[ $files_found -lt ${#SAMPLE_IDS[@]} ]]; then
        local missing=$(( ${#SAMPLE_IDS[@]} - files_found ))
        log_warn "$missing/${#SAMPLE_IDS[@]} samples missing abundance files for '$gene_group' — matrices will be incomplete"
    fi

    # Extract gene names from reference CSV (first column is Gene_ID)
    # Single awk pass replaces tail|cut pipeline (1 process instead of 2)
    awk -F',' 'NR>1 && NF>0 {print $1}' "${ref_csv}" > "$tmpdir/gene_names.txt" \
        || { log_error "Failed to extract gene names from $ref_csv"; rm -rf "$tmpdir"; return 1; }
    # O(1) fork via wc vs O(L) bash read loop
    local gene_name_count
    gene_name_count=$(wc -l < "$tmpdir/gene_names.txt")
    if [[ "$gene_name_count" -eq 0 ]]; then
        log_error "No genes found in reference CSV: $ref_csv"
        rm -rf "$tmpdir"
        return 1
    fi
    log_info "Gene names extracted: $gene_name_count lines."

    # Extract ALL 3 count types (coverage, fpkm, tpm) in a SINGLE awk pass per sample.
    # Replaces 3 separate tail|cut pipelines per sample (was 3×N process spawns, now 1×N).
    # Write to $tmpdir (not alongside input) so cleanup is guaranteed on error.
    #
    # Parallel extraction: launch up to THREADS background awk jobs concurrently.
    # Each awk invocation is I/O-bound (reads one TSV, writes 3 small files), so
    # parallelism yields near-linear speedup until disk bandwidth saturates.
    local _max_jobs="${THREADS:-4}"
    local _running=0
    local _awk_failed=0
    # wait -n requires bash >= 4.3; detect once and fall back to bare wait.
    # Fallback limitation: bare `wait` waits for ALL children and always returns 0,
    # so throttling becomes bursty and _awk_failed undercounts. Acceptable since
    # bash < 4.3 is rare and awk failures here are non-critical (data still merges).
    local _has_wait_n=false
    if (( BASH_VERSINFO[0] > 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] >= 3) )); then
        _has_wait_n=true
    fi
    for srr in "${processed_srrs[@]}"; do
        awk -F'\t' -v gc="$GENENAME_COL" -v cov="$COVERAGE_COL" -v fpkm="$FPKM_COL" -v tpm="$TPM_COL" \
            -v outdir="$tmpdir" -v srr="$srr" 'NR > 1 && $gc != "" && $gc != "." && $gc != "-" {
            print $gc "\t" $cov > outdir "/" srr ".cov"
            print $gc "\t" $fpkm > outdir "/" srr ".fpkm"
            print $gc "\t" $tpm > outdir "/" srr ".tpm"
        }' "${srr_to_file[$srr]}" &
        _running=$((_running + 1))
        if (( _running >= _max_jobs )); then
            if $_has_wait_n; then
                wait -n 2>/dev/null || _awk_failed=$((_awk_failed + 1))
            else
                wait || _awk_failed=$((_awk_failed + 1))
            fi
            _running=$((_running - 1))
        fi
    done
    # Drain remaining background jobs (wait -n catches individual exit statuses;
    # bare wait always returns 0 in bash, so _awk_failed undercounts on bash < 4.3)
    while (( _running > 0 )); do
        if $_has_wait_n; then
            wait -n 2>/dev/null || _awk_failed=$((_awk_failed + 1))
        else
            wait || _awk_failed=$((_awk_failed + 1))
        fi
        _running=$((_running - 1))
    done
    if (( _awk_failed > 0 )); then
        log_error "awk extraction failed for $_awk_failed sample(s)"
        rm -rf "$tmpdir"
        return 1
    fi

    for count_type in coverage fpkm tpm; do
        local -a sample_files=()
        local ext
        case "$count_type" in
            coverage) ext="cov" ;;
            fpkm)     ext="fpkm" ;;
            tpm)      ext="tpm" ;;
        esac

        # Build sample_files and matched_srrs in a single pass (was two identical loops)
        local -a matched_srrs=()
        for srr in "${processed_srrs[@]}"; do
            local extracted="$tmpdir/${srr}.${ext}"
            if [[ -f "$extracted" ]]; then
                sample_files+=("$extracted")
                matched_srrs+=("$srr")
            fi
        done

        if [[ ${#sample_files[@]} -eq 0 ]]; then
            log_warn "No sample files for $count_type in $gene_group, skipping matrix"
            continue
        fi

        # NOTE: Filename uses "geneName" (camelCase) while the CSV header column is "GeneName" (PascalCase).
        # build_input_path() in 0_shared_config.R maps gene_type=="Shortened_Name" -> "geneName" to match this convention.
        local output_geneName_SRR_csv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_SRR${MASTER_SUFFIX}.csv"
        local output_geneName_Organ_csv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_Organ${MASTER_SUFFIX}.csv"

        log_info "Creating SRR + Organ matrices: ${output_geneName_SRR_csv##*/}"

        printf "%s\n" "${sample_files[@]}" > "$tmpdir/sample_files_list.txt"

        # Call matrix_builder.py ONCE; reuse body for both SRR and Organ header variants
        local matrix_body="$tmpdir/matrix_body_${count_type}.txt"
        python3 "$UTILITIES_DIR/matrix_builder.py" "$tmpdir/gene_names.txt" "$tmpdir/sample_files_list.txt" \
            > "$matrix_body" \
            || { log_error "matrix_builder.py failed for $group_name"; rm -rf "$tmpdir"; return 1; }

        # Read matrix body once into variable, write to both CSVs (avoids 2 cat forks).
        local _body
        _body=$(<"$matrix_body")
        rm -f "$matrix_body"

        # SRR header + body (use matched_srrs to align with matrix body columns)
        {
            printf "GeneName"
            for srr in "${matched_srrs[@]}"; do printf ",%s" "$srr"; done
            printf "\n%s\n" "$_body"
        } > "$output_geneName_SRR_csv"

        # Organ header + body
        {
            printf "GeneName"
            for srr in "${matched_srrs[@]}"; do printf ",%s" "${SRR_TO_ORGAN[$srr]:-Unknown}"; done
            printf "\n%s\n" "$_body"
        } > "$output_geneName_Organ_csv"
        # sample_files are in $tmpdir; cleaned by rm -rf "$tmpdir" at function exit
        log_info "Completed $count_type matrix generation"
    done

    rm -rf "$tmpdir"
    log_info "Completed processing for $group_name"
}

# ===============================================
# FULL TRANSCRIPTOME MATRIX
# ===============================================
# Build a full-transcriptome matrix (all genes from abundance files) so that
# WGCNA and genome-wide analyses can consume StringTie data.
# Uses MASTER_REFERENCE as the gene group name to match build_input_path() conventions.

build_full_transcriptome_matrix() {
    local group_name
    group_name=$(get_output_folder_name "$MASTER_REFERENCE")
    log_info "========================================"
    log_info "Building full-transcriptome matrix: $group_name"

    # Collect the UNION of gene IDs from ALL sample abundance files.
    # In de novo mode, StringTie omits zero-coverage transcripts, so any single
    # file may be missing genes that are expressed in other samples.
    # In ref-guided mode all samples share the same gene set, but taking the
    # union keeps this robust if any file is truncated or filtered.
    local tmp_csv
    tmp_csv=$(mktemp --suffix=.csv)
    # NOTE: Do NOT use 'trap ... RETURN' here — merge_group_counts() is called
    # below, and nested RETURN traps clobber each other in bash. Use explicit rm.
    echo "Gene_ID" > "$tmp_csv"

    # Collect existing abundance file paths (avoids spawning awk on missing files)
    local -a _abund_files=()
    for srr in "${SAMPLE_IDS[@]}"; do
        local file_path="$INPUTS_DIR/$MASTER_REFERENCE/$srr/${srr}_${MASTER_REFERENCE}${ABUNDANCE_SUFFIX}"
        [[ -f "$file_path" ]] && _abund_files+=("$file_path")
    done
    local files_found=${#_abund_files[@]}

    if (( files_found > 0 )); then
        # Single awk invocation across ALL files (replaces N separate awk spawns)
        awk -F'\t' -v c="$GENENAME_COL" 'NR>1 && FNR>1 && $c!="" && $c!="." && $c!="-" {print $c}' \
            "${_abund_files[@]}" >> "$tmp_csv" \
            || { log_error "awk gene extraction failed"; rm -f "$tmp_csv"; return 1; }
    fi

    if [[ "$files_found" -eq 0 ]]; then
        log_warn "No abundance files found - skipping full-transcriptome matrix"
        rm -f "$tmp_csv"
        return 1
    fi

    if [[ $files_found -lt ${#SAMPLE_IDS[@]} ]]; then
        local missing=$(( ${#SAMPLE_IDS[@]} - files_found ))
        log_warn "Full-transcriptome: $missing/${#SAMPLE_IDS[@]} samples missing abundance files"
    fi

    # De-duplicate gene IDs AND count in single awk pass — eliminates separate wc -l fork.
    # awk writes dedup output to tmp_dedup, prints only the count to stdout.
    local tmp_dedup gene_count
    tmp_dedup=$(mktemp --suffix=.csv)
    gene_count=$(awk -v out="$tmp_dedup" '
        NR==1 { print > out; next }
        !seen[$0]++ { n++; print > out }
        END { print n+0 }
    ' "$tmp_csv") \
        || { log_error "awk dedup failed"; rm -f "$tmp_csv" "$tmp_dedup"; return 1; }
    mv "$tmp_dedup" "$tmp_csv" \
        || { log_error "mv dedup failed"; rm -f "$tmp_csv" "$tmp_dedup"; return 1; }
    log_info "Full transcriptome: $gene_count genes (union from $files_found samples)"

    # Reuse merge_group_counts with MASTER_REFERENCE as the gene group name
    if merge_group_counts "$MASTER_REFERENCE" "$tmp_csv"; then
        log_info "Full-transcriptome matrix complete"
    else
        log_warn "Failed to build full-transcriptome matrix"
    fi
    rm -f "$tmp_csv"
}

# ===============================================
# MAIN EXECUTION
# ===============================================

# Centralized gene groups CSV directory
GENE_GROUPS_CSV_DIR="${GENE_GROUPS_DIR:-${GENE_GROUPS_CSV_DIR:-$SCRIPT_DIR/../../../../inputs/3_post_proc_inputs/gene_groups_csv}}"

log_info "Gene groups CSV directory: $GENE_GROUPS_CSV_DIR"
log_step "Starting count matrix generation for ${#GENE_GROUPS[@]} gene groups"

# Pre-build CSV lookup map: single find call replaces N per-group find spawns
declare -A _CSV_LOOKUP=()
if [[ -d "$GENE_GROUPS_CSV_DIR" ]]; then
    while IFS= read -r _csv_path; do
        # Pure bash: strip directory + .csv suffix (avoids basename subshell per CSV file)
        _csv_base="${_csv_path##*/}"; _csv_base="${_csv_base%.csv}"
        # First match wins (skip duplicates)
        [[ -z "${_CSV_LOOKUP[$_csv_base]+x}" ]] && _CSV_LOOKUP["$_csv_base"]="$_csv_path"
    done < <(find "$GENE_GROUPS_CSV_DIR" -maxdepth 3 -name "*.csv" -type f 2>/dev/null)
fi

# Build full-transcriptome matrix first (enables WGCNA and genome-wide analyses)
build_full_transcriptome_matrix

for gene_group in "${GENE_GROUPS[@]}"; do
    log_info "========================================"
    log_info "Processing gene group: $gene_group"

    # O(1) lookup from pre-built map (replaces per-group find subprocess)
    REF_CSV="${_CSV_LOOKUP[$gene_group]:-${GENE_GROUPS_CSV_DIR}/${gene_group}.csv}"

    if [[ -z "$REF_CSV" || ! -f "$REF_CSV" ]]; then
        log_error "Reference CSV not found: ${GENE_GROUPS_CSV_DIR}/${gene_group}.csv, skipping $gene_group"
        continue
    fi

    # O(1) fork via wc vs O(L) bash read loop
    # NOTE: no 'local' here — this runs at top-level (outside any function)
    _ref_lines=$(wc -l < "$REF_CSV")
    log_info "Found reference CSV with $((_ref_lines - 1)) genes"

    if merge_group_counts "$gene_group" "$REF_CSV"; then
        log_info "Successfully processed $gene_group"
    else
        log_warn "Failed to process $gene_group"
    fi
done

log_info "========================================"
log_info "$_METHOD_LABEL matrix generation completed"
log_info "Output directory: $OUT_DIR"
