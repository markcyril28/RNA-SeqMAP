#!/bin/bash
#===============================================================================
# PIPELINE UTILITIES - RNA-SEQ ANALYSIS
# Shared functions for post-processing scripts
# Note: Logging functions are provided by logging_utils.sh
#===============================================================================

#===============================================================================
# CSV PARSING
#===============================================================================

# Parse CSV and extract SRR entries as SRR_ID:Organ format
# Usage: mapfile -t SRR_LIST < <(parse_srr_csv "path/to/file.csv")
parse_srr_csv() {
    local csv_file="$1"
    [[ ! -f "$csv_file" ]] && { echo "Warning: CSV not found: $csv_file" >&2; return; }
    
    while IFS=',' read -r srr_id organ notes || [[ -n "$srr_id" ]]; do
        [[ "$srr_id" =~ ^#.*$ || "$srr_id" == "SRR_ID" || -z "$srr_id" ]] && continue
        echo "${srr_id}:${organ}"
    done < "$csv_file"
}

#===============================================================================
# ANALYSIS SCRIPT MAPPING
#===============================================================================

# Map analysis name to R script filename
# NOTE: Preprocessing scripts (Tximport_Salmon, Tximport_RSEM, Stringtie_Matrix)
#       are now handled via get_preprocessing_script() in preprocessing/ folder
# NOTE: Matrix_Creation uses method-specific scripts for M3/M4/M5 (see
#       get_matrix_creation_script()); this table provides the fallback.
get_analysis_script() {
    local -A scripts=(
        ["Matrix_Creation"]="3_Matrix_Creation.R"
        ["Basic_Heatmap"]="4_Basic_Heatmap.R"
        ["Heatmap_with_CV"]="5_Heatmap_with_CV.R"
        ["BarGraph"]="6_BarGraph.R"
        ["Coexpression_using_WGCNA"]="7_Coexpression_WGCNA.R"
        ["Differential_Expression"]="8_Differential_Expression.R"
        ["Gene_Set_Enrichment"]="9_Gene_Set_Enrichment.R"
        ["PCA_Dimensionality_Reduction"]="10_PCA_Dimensionality_Reduction.R"
        ["Sample_Correlation_Clustering"]="11_Sample_Correlation_Clustering.R"
        ["Tissue_Specificity"]="12_Tissue_Specificity.R"
    )
    echo "${scripts[$1]:-}"
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

# Analyses that are thread-intensive and must NOT be parallelised across methods.
# These get full THREADS when running sequentially.
_SEQUENTIAL_ANALYSES=(
    "Matrix_Creation"
    "Differential_Expression"
    "Coexpression_using_WGCNA"
    "Gene_Set_Enrichment"
    "PCA_Dimensionality_Reduction"
)

# Analyses that are lightweight figure-generation tasks — safe to parallelise.
_PARALLEL_ANALYSES=(
    "Basic_Heatmap"
    "Heatmap_with_CV"
    "BarGraph"
    "Sample_Correlation_Clustering"
    "Tissue_Specificity"
)

# Check whether an analysis belongs to the parallelisable (figure) set.
# Note: list is inlined because bash cannot export arrays to GNU Parallel subshells.
# Usage: is_figure_analysis "Basic_Heatmap" && echo yes
is_figure_analysis() {
    local a
    for a in Basic_Heatmap Heatmap_with_CV BarGraph Sample_Correlation_Clustering Tissue_Specificity; do
        [[ "$1" == "$a" ]] && return 0
    done
    return 1
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
    mkdir -p "$method_dir"

    pushd "$method_dir" > /dev/null

    export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"
    export METHOD_BASE_DIR="$method_dir"

    # Export method-specific quant directory so R preprocessing scripts use the exact path
    if [[ "$method" == "M4_Salmon_Saf" ]]; then
        export SALMON_QUANT_ROOT="$BASE_DIR/2_ALIGNMENT_RESULTs/M4_Salmon_Saf/Salmon_Quant/$master_ref"
    fi

    export GENE_GROUPS_DIR="$BASE_DIR/inputs/gene_groups_csv"

    # Rebuild arrays from exported strings — bash arrays are not exported to subshells,
    # so GNU Parallel workers arrive with GENE_GROUPS/ANALYSES empty.
    [[ -n "${GENE_GROUPS_STR:-}" ]] && IFS=' ' read -ra GENE_GROUPS <<< "$GENE_GROUPS_STR"
    [[ -n "${ANALYSES_STR:-}" ]]    && IFS=' ' read -ra ANALYSES    <<< "$ANALYSES_STR"

    # Setup temp config files for R scripts (written to the method's post-proc dir)
    printf '%s\n' "${GENE_GROUPS[@]}" > ".gene_groups_temp.txt"
    echo "$master_ref"                > ".master_reference_temp.txt"
    echo "${OVERWRITE_EXISTING:-FALSE}" > ".overwrite_temp.txt"

    export SRR_COMBINED_LIST_STR="${SRR_COMBINED_LIST_STR:-}"

    popd > /dev/null
}

# Run preprocessing (tximport / prepDE / StringTie) for a single method.
# Usage: run_method_preprocessing "method_name" "master_reference"
run_method_preprocessing() {
    local method=$1 master_ref=$2
    local method_dir="$BASE_DIR/3_POST_PROC/$method"

    pushd "$method_dir" > /dev/null

    export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"
    export METHOD_BASE_DIR="$method_dir"

    [[ -n "${ANALYSES_STR:-}" ]] && IFS=' ' read -ra ANALYSES <<< "$ANALYSES_STR"

    # For M3/M4/M5: skip tximport preprocessing when Matrix_Creation is also enabled —
    # the method-specific 3_Matrix_Creation_*.R script supersedes the tximport step.
    local skip_preprocess=false
    if [[ "$method" =~ ^(M3_STAR_Align|M4_Salmon_Saf|M5_RSEM_Bowtie2)$ ]]; then
        printf '%s\n' "${ANALYSES[@]}" | grep -q "^Matrix_Creation$" && skip_preprocess=true
    fi

    local preprocess_path
    preprocess_path=$(get_preprocessing_script "$method")
    if [[ "$skip_preprocess" == "true" ]]; then
        log_info "Skipping preprocessing for $method — Matrix_Creation will handle import"
    elif [[ -n "$preprocess_path" && -f "$preprocess_path" ]]; then
        log_info "Running preprocessing: $(basename "$preprocess_path")"
        if [[ "$preprocess_path" == *.R ]]; then
            run_with_error_capture Rscript "$preprocess_path" || log_error "Failed: preprocessing ($method)"
        else
            run_with_error_capture bash "$preprocess_path" || log_error "Failed: preprocessing ($method)"
        fi
    elif [[ -n "$preprocess_path" ]]; then
        log_warn "Preprocessing script not found: $preprocess_path"
    fi

    popd > /dev/null
}

# Resolve and run a single analysis script inside the method directory.
# Usage: run_single_analysis "method_name" "master_reference" "analysis_name"
run_single_analysis() {
    local method=$1 master_ref=$2 analysis=$3
    local method_dir="$BASE_DIR/3_POST_PROC/$method"

    pushd "$method_dir" > /dev/null

    export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"
    export METHOD_BASE_DIR="$method_dir"

    # Rebuild arrays in case we are inside a GNU Parallel subshell
    [[ -n "${GENE_GROUPS_STR:-}" ]] && IFS=' ' read -ra GENE_GROUPS <<< "$GENE_GROUPS_STR"

    # Skip legacy preprocessing analysis names
    if [[ "$analysis" =~ ^(Tximport_Salmon|Tximport_RSEM|Tximport_STAR|Stringtie_Matrix)$ ]]; then
        popd > /dev/null
        return 0
    fi

    # For Matrix_Creation, use a method-specific script when available
    local script
    if [[ "$analysis" == "Matrix_Creation" ]]; then
        script=$(get_matrix_creation_script "$method")
        if [[ -z "$script" ]]; then
            log_info "Matrix_Creation not applicable for $method — skipping"
            popd > /dev/null
            return 0
        fi
    else
        script=$(get_analysis_script "$analysis")
    fi

    # Resolve script path: check utilities dir, then analysis_modules dir
    local script_path=""
    local _util_dir="${UTILITIES_DIR:-$BASE_DIR/modules/c_post_processing/utilities}"
    local _mods_dir="${ANALYSIS_MODULES_DIR:-$BASE_DIR/modules/c_post_processing/analysis_modules}"
    if [[ -f "$_util_dir/$script" ]]; then
        script_path="$_util_dir/$script"
    elif [[ -f "$_mods_dir/$script" ]]; then
        script_path="$_mods_dir/$script"
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

    popd > /dev/null
}

# Legacy wrapper — runs ALL analyses for a single method sequentially.
# Kept for backward compatibility; the main script now uses the phased approach.
# Usage: run_method_analysis "method_name" "master_reference"
run_method_analysis() {
    local method=$1 master_ref=$2

    log_step "Processing Method: $method"
    setup_method_env "$method" "$master_ref"
    run_method_preprocessing "$method" "$master_ref"

    [[ -n "${ANALYSES_STR:-}" ]] && IFS=' ' read -ra ANALYSES <<< "$ANALYSES_STR"

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
    # Export logging functions (from logging_utils.sh)
    export -f log log_info log_warn log_error log_step timestamp
    # Export error capture (from logging_utils.sh)
    export -f run_with_error_capture capture_stderr_errors strip_ansi_stream 2>/dev/null || true
    # Export error/warning regex patterns used by capture_stderr_errors
    export _ERROR_PATTERN _WARN_PATTERN 2>/dev/null || true
    # Export pipeline functions
    export -f run_method_analysis run_single_analysis setup_method_env run_method_preprocessing
    export -f is_figure_analysis get_analysis_script get_matrix_creation_script get_preprocessing_script parse_srr_csv
}
