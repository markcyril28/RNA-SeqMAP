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
# Usage: get_matrix_creation_script "method_name"
get_matrix_creation_script() {
    case "$1" in
        "M5_RSEM_Bowtie2") echo "3_Matrix_Creation_RSEM.R" ;;
        "M4_Salmon_Saf")   echo "3_Matrix_Creation_Salmon.R" ;;
        "M3_STAR_Align")   echo "3_Matrix_Creation_STAR.R" ;;
        *)                 echo "3_Matrix_Creation.R" ;;  # fallback for unknown/legacy methods
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
# METHOD ANALYSIS RUNNER
#===============================================================================

# Run analysis for a single method
# Usage: run_method_analysis "method_name" "master_reference"
run_method_analysis() {
    local method=$1 master_ref=$2
    local method_dir="$BASE_DIR/3_POST_PROC/$method"
    
    # Create the method output directory if it doesn't exist yet (first post-processing run).
    # The alignment pipeline normally creates this; creating it here allows post-processing
    # to run independently (e.g., when tximport is the first step that creates content).
    mkdir -p "$method_dir"

    log_step "Processing Method: $method"
    pushd "$method_dir" > /dev/null
    
    export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"

    # Export method-specific quant directory so R preprocessing scripts use the exact path
    # rather than reconstructing it from MASTER_REFERENCE (which may differ from fasta_tag).
    if [[ "$method" == "M4_Salmon_Saf" ]]; then
        export SALMON_QUANT_ROOT="$BASE_DIR/2_ALIGNMENT_RESULTs/M4_Salmon_Saf/Salmon_Quant/$master_ref"
    fi

    # Export GENE_GROUPS_DIR as absolute path for R scripts
    export GENE_GROUPS_DIR="$BASE_DIR/inputs/gene_groups"
    
    # Rebuild arrays from exported strings first — bash arrays are not exported to subshells,
    # so GNU Parallel workers arrive with GENE_GROUPS/ANALYSES empty.  This must happen
    # BEFORE writing .gene_groups_temp.txt so R scripts receive the correct gene groups.
    [[ -n "${GENE_GROUPS_STR:-}" ]] && IFS=' ' read -ra GENE_GROUPS <<< "$GENE_GROUPS_STR"
    [[ -n "${ANALYSES_STR:-}" ]]    && IFS=' ' read -ra ANALYSES  <<< "$ANALYSES_STR"

    # Setup temp config files for R scripts (written to the method's post-proc dir)
    local modules_dir="."

    printf '%s\n' "${GENE_GROUPS[@]}" > "$modules_dir/.gene_groups_temp.txt"
    echo "$master_ref" > "$modules_dir/.master_reference_temp.txt"
    echo "${OVERWRITE_EXISTING:-FALSE}" > "$modules_dir/.overwrite_temp.txt"
    
    # Re-export SRR_COMBINED_LIST_STR for R scripts (ensures it's available in subprocesses)
    export SRR_COMBINED_LIST_STR="${SRR_COMBINED_LIST_STR:-}"
    
    # Run method-specific preprocessing if needed.
    # For M3/M4/M5: skip tximport preprocessing when Matrix_Creation is also enabled —
    # the method-specific 3_Matrix_Creation_*.R script performs the same import and
    # supersedes the tximport preprocessing step.
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
            run_with_error_capture Rscript "$preprocess_path" || log_error "Failed: preprocessing"
        else
            run_with_error_capture bash "$preprocess_path" || log_error "Failed: preprocessing"
        fi
    elif [[ -n "$preprocess_path" ]]; then
        log_warn "Preprocessing script not found: $preprocess_path"
    fi
    
    # Run each enabled analysis
    for analysis in "${ANALYSES[@]}"; do
        [[ -z "$analysis" ]] && continue
        
        # Skip preprocessing analyses entirely - they run via get_preprocessing_script()
        if [[ "$analysis" == "Tximport_Salmon" || "$analysis" == "Tximport_RSEM" || \
              "$analysis" == "Tximport_STAR" || "$analysis" == "Stringtie_Matrix" ]]; then
            continue
        fi
        
        # For Matrix_Creation, use a method-specific script when available
        local script
        if [[ "$analysis" == "Matrix_Creation" ]]; then
            script=$(get_matrix_creation_script "$method")
        else
            script=$(get_analysis_script "$analysis")
        fi

        # Check for script in utilities directory first (for shell scripts like stringtie_matrix_builder.sh)
        local script_path=""
        if [[ -f "$UTILITIES_DIR/$script" ]]; then
            script_path="$UTILITIES_DIR/$script"
        elif [[ -f "$ANALYSIS_MODULES_DIR/$script" ]]; then
            script_path="$ANALYSIS_MODULES_DIR/$script"
        elif [[ -f "$modules_dir/${script%.R}.R" ]]; then
            script_path="$modules_dir/${script%.R}.R"
        fi
        
        if [[ -n "$script_path" && -f "$script_path" ]]; then
            log_info "Running: $analysis"
            if [[ "$script_path" == *.sh ]]; then
                run_with_error_capture bash "$script_path" || log_error "Failed: $analysis"
            else
                run_with_error_capture Rscript "$script_path" || log_error "Failed: $analysis"
            fi
        else
            log_warn "Script not found for: $analysis"
        fi
    done
    
    popd > /dev/null
    log_info "Method $method complete"
}

#===============================================================================
# EXPORT FUNCTIONS FOR PARALLEL
#===============================================================================

export_utils_for_parallel() {
    # Export logging functions (from logging_utils.sh)
    export -f log log_info log_warn log_error log_step timestamp
    # Export error capture (from logging_utils.sh)
    export -f run_with_error_capture capture_stderr_errors 2>/dev/null || true
    # Export pipeline functions
    export -f run_method_analysis get_analysis_script get_matrix_creation_script get_preprocessing_script parse_srr_csv
}
