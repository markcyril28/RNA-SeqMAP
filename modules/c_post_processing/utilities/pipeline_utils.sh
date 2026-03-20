#!/bin/bash
#===============================================================================
# PIPELINE UTILITIES - RNA-SEQ ANALYSIS
# Shared functions for post-processing scripts
# Note: Logging functions are provided by logging_utils.sh
#===============================================================================

# Guard against double-sourcing
[[ "${PIPELINE_UTILS_SOURCED:-}" == "true" ]] && return 0
export PIPELINE_UTILS_SOURCED="true"

# Pre-compute utility directory paths at module load (avoids ~15-30 local var assignments
# per run_single_analysis call × N methods × N analyses). O(1) lookup thereafter.
_PIPELINE_UTIL_DIR="${UTILITIES_DIR:-${BASE_DIR:-.}/modules/c_post_processing/utilities}"
_PIPELINE_MODS_DIR="${ANALYSIS_MODULES_DIR:-${BASE_DIR:-.}/modules/c_post_processing/analysis_modules}"

#===============================================================================
# CSV PARSING
#===============================================================================

# Parse CSV and extract SRR entries as SRR_ID:Organ format. O(n) single awk pass.
# Replaces while-read loop (avoids per-line bash overhead: fork+IFS parse per row).
# Usage: mapfile -t SRR_LIST < <(parse_srr_csv "path/to/file.csv")
parse_srr_csv() {
    local csv_file="$1"
    [[ ! -f "$csv_file" ]] && { log_warn "CSV not found: $csv_file"; return 1; }

    awk -F',' '
        NR == 1 { next }
        /^#/ || /^[[:space:]]*$/ { next }
        $1 == "SRR_ID" { next }
        { gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
          if ($1 != "") print $1 ":" $2 }
    ' "$csv_file"
}

#===============================================================================
# ARRAY REBUILD HELPER (bash arrays cannot be exported to subshells)
#===============================================================================

# Rebuild GENE_GROUPS and ANALYSES arrays from exported string representations.
# GNU Parallel workers arrive with arrays empty — call this in any parallel entry point.
# Guard: skip if already rebuilt in this process (saves ~100 redundant splits in large runs).
_rebuild_exported_arrays() {
    [[ "${_ARRAYS_REBUILT:-}" == "$$" ]] && return 0
    # Save/restore IFS so inline IFS=' ' does not leak into caller's environment.
    local _saved_ifs="${IFS:-}"
    [[ -n "${GENE_GROUPS_STR:-}" ]] && { IFS=' ' read -ra GENE_GROUPS <<< "$GENE_GROUPS_STR"; }
    [[ -n "${ANALYSES_STR:-}" ]]    && { IFS=' ' read -ra ANALYSES    <<< "$ANALYSES_STR"; }
    IFS="$_saved_ifs"
    _ARRAYS_REBUILT="$$"
}

#===============================================================================
# ANALYSIS SCRIPT MAPPING
#===============================================================================

# Map analysis name to R script filename. O(1) via case pattern match.
# Uses case statement (not associative array) so the function works correctly
# when exported to GNU Parallel subshells (bash cannot export associative arrays).
# NOTE: Preprocessing scripts (Tximport_Salmon, Tximport_RSEM, Stringtie_Matrix)
#       are now handled via get_preprocessing_script() in preprocessing/ folder
# NOTE: Matrix_Creation uses method-specific scripts for M3/M4/M5 (see
#       get_matrix_creation_script()); this provides the fallback.
get_analysis_script() {
    case "$1" in
        "Matrix_Creation")                echo "3_Matrix_Creation.R" ;;
        "Basic_Heatmap")                  echo "4_Basic_Heatmap.R" ;;
        "Heatmap_with_CV")               echo "5_Heatmap_with_CV.R" ;;
        *)                                echo "" ;;
    esac
}

# Return the method-specific Matrix_Creation script (M3/M4/M5) or fallback.
# M1/M2 (HISAT2+StringTie) do not use Matrix_Creation — preprocessing is handled
# by stringtie_matrix_builder.sh (M2) or prepde_matrix_linker.sh (M1).
# Usage: get_matrix_creation_script "method_name"
get_matrix_creation_script() {
    case "$1" in
        "M5_RSEM_Bowtie2")       echo "3_Matrix_Creation_RSEM.R" ;;
        "M4_Salmon_Saf")         echo "3_Matrix_Creation_Salmon.R" ;;
        "M3_STAR_Align")         echo "3_Matrix_Creation_STAR.R" ;;
        "M1_HISAT2_RefGuided"|"M2_HISAT2_DeNovo") echo "" ;;  # Not applicable — use preprocessing script instead
        *)                       echo "3_Matrix_Creation.R" ;;  # fallback for unknown/legacy methods
    esac
}

# Get method-specific preprocessing script
# Usage: get_preprocessing_script "method_name"
# Returns the full path to the script that converts raw quant output to count matrices
get_preprocessing_script() {
    local method=$1
    local PREPROCESSING_DIR="${BASE_DIR}/modules/c_post_processing/preprocessing"
    
    case "$method" in
        "M1_HISAT2_RefGuided")
            # M1: prepDE.py already ran during alignment; linker validates and stages the integer count matrices
            echo "${PREPROCESSING_DIR}/HISAT2/prepde_matrix_linker.sh"
            ;;
        "M2_HISAT2_DeNovo")
            # M2: builds FPKM/TPM/Coverage matrices from StringTie abundance files
            echo "${PREPROCESSING_DIR}/HISAT2/stringtie_matrix_builder.sh"
            ;;
        "M3_STAR_Align")
            # STAR+Salmon: runs tximport on 6_salmon/quant outputs to create standardized matrices
            echo "${PREPROCESSING_DIR}/STAR/tximport_star_to_matrices.R"
            ;;
        "M4_Salmon_Saf")
            echo "${PREPROCESSING_DIR}/Salmon/tximport_salmon_to_matrices.R"
            ;;
        "M5_RSEM_Bowtie2")
            echo "${PREPROCESSING_DIR}/RSEM/tximport_rsem_to_matrices.R"
            ;;
        *)
            echo ""
            ;;
    esac
}

#===============================================================================
# ANALYSIS CLASSIFICATION
#===============================================================================

# Check whether an analysis belongs to the parallelisable (figure) set.
# Matrix_Creation remains a thread-intensive (non-figure) analysis.
# Note: list is inlined because bash cannot export arrays to GNU Parallel subshells.
# Uses case for O(1) pattern match instead of O(n) loop.
# Usage: is_figure_analysis "Basic_Heatmap" && echo yes
is_figure_analysis() {
    case "$1" in
        Basic_Heatmap|Heatmap_with_CV) return 0 ;;
        *) return 1 ;;
    esac
}

#===============================================================================
# METHOD ANALYSIS RUNNER
#===============================================================================

# Prepare the method directory, temp config files, and run preprocessing.
# Called once per method before any analysis phase.
# Usage: setup_method_env "method_name" "master_reference"
setup_method_env() {
    local method=$1 master_ref=$2
    local method_dir="$BASE_DIR/3_POST_PROC/$method"

    # Create the method output directory if it doesn't exist yet (first post-processing run).
    mkdir -p "$method_dir" || { log_error "Failed to create directory: $method_dir"; return 1; }

    pushd "$method_dir" > /dev/null || { log_error "Cannot cd to $method_dir"; return 1; }

    export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"
    export METHOD_BASE_DIR="$method_dir"

    # Export method-specific quant directory so R preprocessing scripts use the exact path
    if [[ "$method" == "M4_Salmon_Saf" ]]; then
        export SALMON_QUANT_ROOT="$BASE_DIR/2_ALIGNMENT_RESULTs/M4_Salmon_Saf/Salmon_Quant/$master_ref"
    elif [[ "$method" == "M5_RSEM_Bowtie2" ]]; then
        export RSEM_QUANT_ROOT="$BASE_DIR/2_ALIGNMENT_RESULTs/M5_RSEM_Bowtie2/RSEM_Quant_WD/$master_ref"
    fi

    export GENE_GROUPS_DIR="$BASE_DIR/inputs/3_post_proc_inputs/gene_groups_csv"

    # Rebuild arrays from exported strings (bash arrays are not exported to subshells)
    _rebuild_exported_arrays

    # Setup temp config files for R scripts (written to the method's post-proc dir)
    # Consolidate into a single write operation to reduce disk I/O (3 files → 1 atomic write each)
    {
        printf '%s\n' "${GENE_GROUPS[@]}"
    } > ".gene_groups_temp.txt" || { log_error "Failed to write .gene_groups_temp.txt in $method_dir"; popd > /dev/null; return 1; }
    printf '%s\n%s\n' "$master_ref" "${OVERWRITE_EXISTING:-FALSE}" > ".method_config_temp.txt" || { log_error "Failed to write .method_config_temp.txt in $method_dir"; popd > /dev/null; return 1; }
    # Legacy compat: still write individual files for any R scripts that read them directly
    echo "$master_ref"                > ".master_reference_temp.txt"
    echo "${OVERWRITE_EXISTING:-FALSE}" > ".overwrite_temp.txt"

    export SRR_COMBINED_LIST_STR="${SRR_COMBINED_LIST_STR:-}"

    popd > /dev/null || log_warn "popd failed in setup_method_env (was in $method_dir)"
}

# Run preprocessing (tximport / prepDE / StringTie) for a single method.
# Usage: run_method_preprocessing "method_name" "master_reference"
run_method_preprocessing() {
    local method=$1 master_ref=$2
    local method_dir="$BASE_DIR/3_POST_PROC/$method"

    pushd "$method_dir" > /dev/null || { log_error "Cannot cd to $method_dir"; return 1; }

    # Re-export env vars (needed in GNU Parallel subshells where setup_method_env
    # ran in parent). Skip _rebuild_exported_arrays if already in same process.
    # ASSUMPTION: setup_method_env() was called before this function, so MASTER_REFERENCE
    # is already correct. If the same method is called with a different master_ref without
    # re-calling setup_method_env(), MASTER_REFERENCE would be stale.
    if [[ "${CURRENT_METHOD:-}" != "$method" ]]; then
        export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"
        export METHOD_BASE_DIR="$method_dir"
        _rebuild_exported_arrays
    fi

    # For M3/M4/M5: skip tximport preprocessing when Matrix_Creation is also enabled —
    # the method-specific 3_Matrix_Creation_*.R script supersedes the tximport step.
    local skip_preprocess=false
    if [[ "$method" =~ ^(M3_STAR_Align|M4_Salmon_Saf|M5_RSEM_Bowtie2)$ ]]; then
        for _a in "${ANALYSES[@]}"; do
            [[ "$_a" == "Matrix_Creation" ]] && { skip_preprocess=true; break; }
        done
    fi

    local preprocess_path
    preprocess_path=$(get_preprocessing_script "$method")
    if [[ "$skip_preprocess" == "true" ]]; then
        log_info "Skipping preprocessing for $method — Matrix_Creation will handle import"
    elif [[ -n "$preprocess_path" && -f "$preprocess_path" ]]; then
        log_info "Running preprocessing: $(basename "$preprocess_path")"
        if [[ "$preprocess_path" == *.R ]]; then
            run_with_error_capture Rscript "$preprocess_path" || { log_error "Failed: preprocessing ($method)"; popd > /dev/null; return 1; }
        else
            run_with_error_capture bash "$preprocess_path" || { log_error "Failed: preprocessing ($method)"; popd > /dev/null; return 1; }
        fi
    elif [[ -n "$preprocess_path" ]]; then
        log_warn "Preprocessing script not found: $preprocess_path"
        popd > /dev/null || true
        return 1
    fi

    popd > /dev/null || log_warn "popd failed in run_method_preprocessing (was in $method_dir)"
}

# Resolve and run a single analysis script inside the method directory.
# Usage: run_single_analysis "method_name" "master_reference" "analysis_name"
run_single_analysis() {
    local method=$1 master_ref=$2 analysis=$3
    local method_dir="$BASE_DIR/3_POST_PROC/$method"

    pushd "$method_dir" > /dev/null || { log_error "Cannot cd to $method_dir"; return 1; }

    # Re-export env vars only when in a GNU Parallel subshell (CURRENT_METHOD differs).
    # In same-process sequential runs, setup_method_env already set these — skip for speed.
    # ASSUMPTION: setup_method_env() was called before this function, so MASTER_REFERENCE
    # is already correct. If the same method is called with a different master_ref without
    # re-calling setup_method_env(), MASTER_REFERENCE would be stale.
    if [[ "${CURRENT_METHOD:-}" != "$method" ]]; then
        export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"
        export METHOD_BASE_DIR="$method_dir"
        _rebuild_exported_arrays
    fi

    # Skip legacy preprocessing analysis names
    if [[ "$analysis" =~ ^(Tximport_Salmon|Tximport_RSEM|Tximport_STAR|Stringtie_Matrix)$ ]]; then
        popd > /dev/null || true
        return 0
    fi

    # For Matrix_Creation, use a method-specific script when available
    local script
    if [[ "$analysis" == "Matrix_Creation" ]]; then
        script=$(get_matrix_creation_script "$method")
        if [[ -z "$script" ]]; then
            log_info "Matrix_Creation not applicable for $method — skipping"
            popd > /dev/null || true
            return 0
        fi
    else
        script=$(get_analysis_script "$analysis")
    fi

    # Resolve script path: check utilities dir, then analysis_modules dir
    # Uses module-level pre-computed paths (_PIPELINE_UTIL_DIR, _PIPELINE_MODS_DIR)
    # to avoid repeated local variable allocation. O(1) path resolution.
    local script_path=""
    if [[ -f "${_PIPELINE_UTIL_DIR}/$script" ]]; then
        script_path="${_PIPELINE_UTIL_DIR}/$script"
    elif [[ -f "${_PIPELINE_MODS_DIR}/$script" ]]; then
        script_path="${_PIPELINE_MODS_DIR}/$script"
    fi

    if [[ -n "$script_path" && -f "$script_path" ]]; then
        log_info "Running: $analysis ($method)"
        if [[ "$script_path" == *.sh ]]; then
            run_with_error_capture bash "$script_path" || log_error "Failed: $analysis ($method)"
        else
            run_with_error_capture Rscript "$script_path" || log_error "Failed: $analysis ($method)"
        fi
    else
        log_warn "Script not found for: $analysis ($method)"
    fi

    popd > /dev/null || log_warn "popd failed in run_single_analysis (was in $method_dir)"
}

# Legacy wrapper — runs ALL analyses for a single method sequentially.
# Kept for backward compatibility; the main script now uses the phased approach.
# Usage: run_method_analysis "method_name" "master_reference"
run_method_analysis() {
    local method=$1 master_ref=$2

    log_step "Processing Method: $method"
    setup_method_env "$method" "$master_ref"
    run_method_preprocessing "$method" "$master_ref"

    _rebuild_exported_arrays

    for analysis in "${ANALYSES[@]}"; do
        [[ -z "$analysis" ]] && continue
        run_single_analysis "$method" "$master_ref" "$analysis"
    done

    log_info "Method $method complete"
}

#===============================================================================
# EXPORT FUNCTIONS FOR PARALLEL
#===============================================================================

export_utils_for_parallel() {
    # Guard: skip if already exported this session (saves ~20 export -f calls per dataset iteration)
    [[ "${_PARALLEL_UTILS_EXPORTED:-}" == "true" ]] && return 0
    # Export logging functions (from logging_utils.sh)
    export -f log log_info log_warn log_error log_step timestamp
    # Export error capture (from logging_utils.sh)
    export -f run_with_error_capture capture_stderr_errors strip_ansi_stream 2>/dev/null || true
    # Export error/warning regex patterns used by capture_stderr_errors
    export _ERROR_PATTERN _WARN_PATTERN 2>/dev/null || true
    # Export pipeline functions
    export -f _rebuild_exported_arrays
    export -f run_method_analysis run_single_analysis setup_method_env run_method_preprocessing
    export -f is_figure_analysis get_analysis_script get_matrix_creation_script get_preprocessing_script parse_srr_csv
    _PARALLEL_UTILS_EXPORTED="true"
}
