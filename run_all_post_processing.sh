#!/bin/bash
#===============================================================================
# MASTER POST-PROCESSING SCRIPT - RNA-SEQ ANALYSIS PIPELINE
# Runs selected analyses across all configured methods and gene groups
#===============================================================================

#set -euo pipefail

#===============================================================================
# DIRECTORY PATHS and SOURCE UTILITIES
#===============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$SCRIPT_DIR"  # HeatSeq is the project root
ANALYSIS_MODULES_DIR="$BASE_DIR/modules/c_post_processing/analysis_modules"
GENE_GROUPS_DIR="$BASE_DIR/0_INPUTs/gene_groups"
SRR_CSV_DIR="$BASE_DIR/0_INPUTs/SRR_csv"
UTILITIES_DIR="$BASE_DIR/modules/c_post_processing/utilities"
LOGGING_UTILS="$BASE_DIR/modules/logging/logging_utils.sh"

# Source logging utilities first (provides log_info, log_step, log_error, etc.)
source "$LOGGING_UTILS"
# Source pipeline utilities (provides parse_srr_csv, run_method_analysis, etc.)
source "$UTILITIES_DIR/pipeline_utils.sh"

#===============================================================================
# GLOBAL SETTINGS
#===============================================================================

# System resources (32GB RAM available - prioritize accuracy over memory)
THREADS=12
ENABLE_GPU="false"              # To be fixed 
ENABLE_GNU_PARALLEL="FALSE"      # Enable parallel processing for faster execution
DESIRED_CPU_PER_JOB=4
AVAILABLE_RAM_GB=32
GPU_VRAM_GB=8

# Logging and output
CLEAR_LOGS="FALSE"
CLEAR_OUTPUT_FOLDER="FALSE"

# Skip existing outputs when CLEAR_OUTPUT_FOLDER is FALSE
# TRUE = regenerate all outputs, FALSE = skip if output already exists
if [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]]; then
    OVERWRITE_EXISTING="TRUE"
else
    OVERWRITE_EXISTING="FALSE"
fi

# Export RAM info for R scripts
export AVAILABLE_RAM_GB GPU_VRAM_GB OVERWRITE_EXISTING

#===============================================================================
# DATASETS (comment/uncomment to enable/disable)
#===============================================================================

# Master reference (uncomment ONE)
MASTER_REFERENCES=(
    #"All_Smel_Genes"
    "Eggplant_V4.1_transcripts.function"
)

# Extract first enabled reference
for ref in "${MASTER_REFERENCES[@]}"; do
    MASTER_REFERENCE="$ref"
    break
done

# Gene groups to analyze
# NOTE: For HISAT2 methods, ensure these groups have pre-built matrices in
#       count_matrices_from_stringtie/{gene_group}/
GENE_GROUPS=(
    #"SmelDMPs"
    "SmelDMPs_with_1_18s_rRNA"
    #"SmelDMPs_with_2_18s_rRNA"
    #"SmelDMPs_with_SmelCyclo"
    #"SmelGRF-GIFs"              # Needs stringtie_matrix_builder.sh to be run first
    #"SmelGRF-GIFs_with_1_18s_rRNA"
    #"SmelGRFs"                 # Alternative: already has StringTie matrices
    #"SmelGIFs"                 # Alternative: already has StringTie matrices
    #"Selected_GRF_GIF_Genes_vAll_GIFs"
)

# Each name corresponds to a CSV file in 0_INPUTs/SRR_csv/
# CSV format: SRR_ID,Organ,Notes
SRR_DATASETS=(
    #"PRJNA328564"      # Main Dataset - Eggplant tissue atlas (PRJNA328564)
    "PRJNA328564_selected"
    #"SAMN28540077"     # Chinese Dataset 1 - replicability
    #"SAMN28540068"     # Chinese Dataset 2 - replicability
    "PRJNA865018"     # SET_1: Good Dataset for SmelDMP GEA
    "PRJNA941250"     # SET_2: Good Dataset for SmelDMP GEA
    "PRJNA865018_and_PRJNA941250"
    #"OTHER_SRR_LIST"  # Other miscellaneous samples
)

#===============================================================================
# CONFIGURATION (comment/uncomment to enable/disable)
#===============================================================================

# Alignment Methods
METHODS=(
    #"M1_HISAT2_RefGuided"
    #"M2_HISAT2_DeNovo"
    #"M3_STAR_Align"
    #"M4_Salmon_Saf"
    "M5_RSEM_Bowtie2"
)

# Analysis modules (comment/uncomment to enable/disable)
# 
# METHOD-SPECIFIC NOTES:
# ----------------------
# HISAT2 (M1/M2): Matrices are built by stringtie_matrix_builder.sh during preprocessing.
#                 Use Stringtie_Matrix for matrix creation. Do NOT use Tximport_Salmon or Tximport_RSEM.
#                 Heatmaps read from: count_matrices_from_stringtie/
#
# Salmon (M4):    Requires Tximport_Salmon for matrix creation (auto-run as preprocessing).
#                 Heatmaps read from: count_matrices_from_Salmon_Quant/
#
# RSEM (M5):      Requires Tximport_RSEM for matrix creation (auto-run as preprocessing).
#                 Heatmaps read from: count_matrices_from_RSEM_Quant/

ANALYSES=(
    "Stringtie_Matrix"         # For HISAT2-based methods (M1/M2) - builds matrices from StringTie output

    "Tximport_Salmon"          # For Salmon-based methods (M4) - auto-run based on method
    "Tximport_RSEM"            # For RSEM-based methods (M5) - auto-run based on method
    "Matrix_Creation"          # For tximport methods only (M4/M5) - NOT for HISAT2!

    "Basic_Heatmap"             # Works with all methods (auto-detects input path)
    #"Heatmap_with_CV"           # Works with all methods (auto-detects input path)
    #"BarGraph"
    #"Coexpression_using_WGCNA"  # Works with all methods (auto-detects input path)
    #"Differential_Expression"   # DESeq2-based differential expression
    #"Gene_Set_Enrichment"
    #"PCA_Dimensionality_Reduction"  # PCA, t-SNE, UMAP
    #"Sample_Correlation_Clustering" # Sample QC and clustering
    #"Tissue_Specificity"
)


#===============================================================================
# INITIALIZATION
#===============================================================================

# Build combined SRR list from selected datasets
SRR_COMBINED_LIST=()
for dataset in "${SRR_DATASETS[@]}"; do
    csv_file="$SRR_CSV_DIR/${dataset}.csv"
    if [[ -f "$csv_file" ]]; then
        mapfile -t -O "${#SRR_COMBINED_LIST[@]}" SRR_COMBINED_LIST < <(parse_srr_csv "$csv_file")
    else
        echo "Warning: SRR CSV not found: $csv_file"
    fi
done

# Validate arrays (after SRR list is built)
[[ -z "$MASTER_REFERENCE" ]] && { echo "ERROR: No master reference configured"; exit 1; }
[[ ${#METHODS[@]} -eq 0 ]] && { echo "ERROR: No methods configured"; exit 1; }
[[ ${#GENE_GROUPS[@]} -eq 0 ]] && { echo "ERROR: No gene groups configured"; exit 1; }
[[ ${#SRR_COMBINED_LIST[@]} -eq 0 ]] && { echo "ERROR: No SRR samples loaded from datasets"; exit 1; }

# Calculate parallel jobs
JOBS=1
[[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && JOBS=$((THREADS / DESIRED_CPU_PER_JOB)) && [[ $JOBS -lt 1 ]] && JOBS=1

# Initialize conda environment
eval "$(conda shell.bash hook)"
conda activate gea 2>/dev/null || echo "Warning: conda env 'gea' not found, using current env"

# Source optional shared utilities
[[ -f "$BASE_DIR/modules/config.sh" ]] && source "$BASE_DIR/modules/config.sh"
[[ -f "$BASE_DIR/modules/shared_utils.sh" ]] && source "$BASE_DIR/modules/shared_utils.sh"

# Setup logging directories with ABSOLUTE paths (critical for subprocesses that change directories)
LOG_DIR="$BASE_DIR/3_POST_PROC/logs/log_files"
TIME_DIR="$BASE_DIR/3_POST_PROC/logs/time_logs"
SPACE_DIR="$BASE_DIR/3_POST_PROC/logs/space_logs"
SPACE_TIME_DIR="$BASE_DIR/3_POST_PROC/logs/space_time_logs"
ERROR_WARN_DIR="$BASE_DIR/3_POST_PROC/logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="$BASE_DIR/3_POST_PROC/logs/software_catalogs"
GPU_LOG_DIR="$BASE_DIR/3_POST_PROC/logs/gpu_log"
export LOG_DIR TIME_DIR SPACE_DIR SPACE_TIME_DIR ERROR_WARN_DIR SOFTWARE_CATALOG_DIR GPU_LOG_DIR

# Initialize logging system
setup_logging "$CLEAR_LOGS"

# Export actual log FILE paths after setup_logging has set them (for subprocesses)
export LOG_FILE TIME_FILE SPACE_FILE SPACE_TIME_FILE ERROR_WARN_FILE SOFTWARE_FILE GPU_LOG_FILE

# Map analysis name to output folder name
get_output_folder_name() {
    local analysis="$1"
    case "$analysis" in
        "Matrix_Creation")                  echo "0_Matrix_Creation" ;;
        "Basic_Heatmap")                    echo "I_Basic_Heatmap" ;;
        "Heatmap_with_CV")                  echo "II_Heatmap_with_CV" ;;
        "BarGraph")                         echo "III_Bar_Graphs" ;;
        "Coexpression_using_WGCNA")         echo "III_Coexpression_WGCNA" ;;
        "Differential_Expression")          echo "V_Differential_Expression" ;;
        "Gene_Set_Enrichment")              echo "VI_Gene_Set_Enrichment" ;;
        "PCA_Dimensionality_Reduction")     echo "VII_Dimensionality_Reduction" ;;
        "Sample_Correlation_Clustering")    echo "VIII_Sample_Correlation" ;;
        "Tissue_Specificity")               echo "IX_Tissue_Specificity" ;;
        *)                                  echo "" ;;  # No output folder for preprocessing analyses
    esac
}

# Clear output folders if requested (only for enabled analyses, not entire output dir)
if [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]]; then
    log_info "Clearing output folders for enabled analyses..."
    for method in "${METHODS[@]}"; do
        output_base="$BASE_DIR/3_POST_PROC/$method/Figure_Outputs"
        for analysis in "${ANALYSES[@]}"; do
            folder_name=$(get_output_folder_name "$analysis")
            if [[ -n "$folder_name" && -d "$output_base/$folder_name/$MASTER_REFERENCE" ]]; then
                log_info "  Clearing: $method/$folder_name/$MASTER_REFERENCE"
                rm -rf "$output_base/$folder_name/$MASTER_REFERENCE"/* 2>/dev/null || true
            fi
        done
    done
fi

# Export environment variables for R scripts and subprocesses
export BASE_DIR THREADS ENABLE_GPU ANALYSIS_MODULES_DIR GENE_GROUPS_DIR UTILITIES_DIR SRR_CSV_DIR

# Export arrays as strings for subprocesses (needed by stringtie_matrix_builder.sh)
export GENE_GROUPS_STR="${GENE_GROUPS[*]}"
export ANALYSES_STR="${ANALYSES[*]}"
export SRR_DATASETS_STR="${SRR_DATASETS[*]}"

#===============================================================================
# MAIN EXECUTION
#===============================================================================

log_step "Starting Post-Processing Pipeline"
log_info "Threads: $THREADS | GPU: $ENABLE_GPU (${GPU_VRAM_GB}GB VRAM) | Parallel: $ENABLE_GNU_PARALLEL (Jobs: $JOBS)"
log_info "Master Reference: $MASTER_REFERENCE"
log_info "Methods: ${METHODS[*]}"
log_info "Gene Groups: ${GENE_GROUPS[*]}"
log_info "SRR Datasets: ${SRR_DATASETS[*]} (${#SRR_COMBINED_LIST[@]} total samples)"
log_info "Analyses: ${ANALYSES[*]}"
log_info "Note: Processing each dataset separately with combined output naming (e.g., GeneGroup_in_Dataset)"

# Process each dataset separately (outer loop)
for dataset in "${SRR_DATASETS[@]}"; do
    csv_file="$SRR_CSV_DIR/${dataset}.csv"
    if [[ ! -f "$csv_file" ]]; then
        log_warn "SRR CSV not found: $csv_file, skipping dataset $dataset"
        continue
    fi
    
    # Build SRR list for this dataset only
    mapfile -t CURRENT_SRR_LIST < <(parse_srr_csv "$csv_file")
    if [[ ${#CURRENT_SRR_LIST[@]} -eq 0 ]]; then
        log_warn "No samples found in $dataset, skipping"
        continue
    fi
    
    # Export current dataset info for subprocesses
    export CURRENT_DATASET="$dataset"
    export SRR_COMBINED_LIST_STR="${CURRENT_SRR_LIST[*]}"
    
    log_step "Processing Dataset: $dataset (${#CURRENT_SRR_LIST[@]} samples)"
    
    # Run methods for this dataset (parallel or sequential)
    if [[ "$ENABLE_GNU_PARALLEL" == "TRUE" && $JOBS -gt 1 ]] && command -v parallel &>/dev/null; then
        log_info "Using GNU Parallel with $JOBS jobs"
        
        # Export for parallel subshells
        export_utils_for_parallel
        export SCRIPT_DIR LOG_FILE ERROR_WARN_FILE RUN_ID ANALYSIS_MODULES_DIR GENE_GROUPS_DIR
        export MASTER_REFERENCE THREADS ENABLE_GPU CURRENT_DATASET
        export GENE_GROUPS_STR="${GENE_GROUPS[*]}" ANALYSES_STR="${ANALYSES[*]}"
        
        printf '%s\n' "${METHODS[@]}" | parallel -j "$JOBS" run_method_analysis {} "$MASTER_REFERENCE"
    else
        log_info "Running methods sequentially for $dataset"
        for method in "${METHODS[@]}"; do
            run_method_analysis "$method" "$MASTER_REFERENCE"
        done
    fi
done

#===============================================================================
# SUMMARY
#===============================================================================

log_step "Post-Processing Complete"
log_info "Datasets processed: ${#SRR_DATASETS[@]} (${SRR_DATASETS[*]})"
log_info "Methods processed: ${#METHODS[@]} per dataset"
log_info "Output naming: GeneGroup_in_Dataset (e.g., ${GENE_GROUPS[0]}_in_${SRR_DATASETS[0]})"
log_info "Log file: $LOG_FILE"
log_info "Time metrics: $TIME_FILE"

if [[ -f "$ERROR_WARN_FILE" && -s "$ERROR_WARN_FILE" ]]; then
    log_info "Errors/Warnings: $(wc -l < "$ERROR_WARN_FILE") (see $ERROR_WARN_FILE)"
else
    log_info "No errors encountered"
fi

echo -e "\n========================================"
echo "Pipeline completed at $(date)"
echo "========================================"