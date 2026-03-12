#!/bin/bash
#===============================================================================
# MASTER POST-PROCESSING SCRIPT - RNA-SEQ ANALYSIS PIPELINE
#===============================================================================

set -o pipefail   # -e/-u omitted intentionally (sourced functions use boolean returns)

# ==============================================================================
# SYSTEM RESOURCES
# ==============================================================================

THREADS=$(nproc)
ENABLE_GPU="FALSE"
ENABLE_GNU_PARALLEL="TRUE"
DESIRED_CPU_PER_JOB=1
AVAILABLE_RAM_GB=24
GPU_VRAM_GB=8

# ==============================================================================
# LOGGING AND OUTPUT
# ==============================================================================

CLEAR_LOGS="TRUE"
CLEAR_OUTPUT_FOLDER="TRUE"

PIPELINE_CONFIGS=(
    # ── Full — Eggplant_V4.1 ──
    "config/3_post_proc_configs/HPC_full_M1_Eggplant_V4.1.sh"      # M1 HISAT2 RefGuided   Eggplant_V4.1 genome
    "config/3_post_proc_configs/HPC_full_M2_Eggplant_V4.1.sh"      # M2 HISAT2 DeNovo      Eggplant_V4.1 transcript
    "config/3_post_proc_configs/HPC_full_M3_Eggplant_V4.1.sh"      # M3 STAR Align         Eggplant_V4.1 genome
    "config/3_post_proc_configs/HPC_full_M4_Eggplant_V4.1.sh"      # M4 Salmon SAF         Eggplant_V4.1 transcript
    "config/3_post_proc_configs/HPC_full_M5_Eggplant_V4.1.sh"      # M5 RSEM Bowtie2       Eggplant_V4.1 transcript

    # ── Full — GPE001970 ──
    "config/3_post_proc_configs/HPC_full_M1_GPE001970.sh"           # M1 HISAT2 RefGuided   GPE001970 genome
    "config/3_post_proc_configs/HPC_full_M2_GPE001970.sh"           # M2 HISAT2 DeNovo      GPE001970 transcript
    "config/3_post_proc_configs/HPC_full_M3_GPE001970.sh"           # M3 STAR Align         GPE001970 genome
    "config/3_post_proc_configs/HPC_full_M4_GPE001970.sh"           # M4 Salmon SAF         GPE001970 transcript
    "config/3_post_proc_configs/HPC_full_M5_GPE001970.sh"           # M5 RSEM Bowtie2       GPE001970 transcript

    # ── Legacy (original unsplit full configs) ──
    #"config/3_post_proc_configs/HPC_full_ref_guided.sh"             # HPC full,  genome (M1+M3)
    #"config/3_post_proc_configs/HPC_full_non_ref_guided.sh"         # HPC full,  transcript (M2+M4+M5)
)

#===============================================================================
# PATHS AND UTILITIES
#===============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$SCRIPT_DIR"
ANALYSIS_MODULES_DIR="$BASE_DIR/modules/c_post_processing/analysis_modules"
GENE_GROUPS_DIR="$BASE_DIR/inputs/gene_groups_csv"
SRR_CSV_DIR="$BASE_DIR/inputs/SRR_csv"
UTILITIES_DIR="$BASE_DIR/modules/c_post_processing/utilities"

source "$BASE_DIR/modules/logging/logging_utils.sh"
source "$UTILITIES_DIR/pipeline_utils.sh"


#===============================================================================
# FUNCTIONS
#===============================================================================

get_output_folder_name() {
    case "$1" in
        "Matrix_Creation")               echo "0_Matrix_Creation" ;;
        "Basic_Heatmap")                 echo "I_Basic_Heatmap" ;;
        "Heatmap_with_CV")               echo "II_Heatmap_with_CV" ;;
        "BarGraph")                      echo "III_Bar_Graphs" ;;
        "Coexpression_using_WGCNA")      echo "IV_Coexpression_WGCNA" ;;
        "Differential_Expression")       echo "V_Differential_Expression" ;;
        "Gene_Set_Enrichment")           echo "VI_Gene_Set_Enrichment" ;;
        "PCA_Dimensionality_Reduction")  echo "VII_PCA" ;;
        "Sample_Correlation_Clustering") echo "VIII_Sample_Clustering" ;;
        "Tissue_Specificity")            echo "IX_Tissue_Specificity" ;;
        *)                               echo "" ;;
    esac
}

#===============================================================================
# INITIALIZATION
#===============================================================================

[[ ${#PIPELINE_CONFIGS[@]} -eq 0 ]] && { log_error "No configs enabled in PIPELINE_CONFIGS"; exit 1; }

eval "$(conda shell.bash hook)"
conda activate gea 2>/dev/null || log_warn "conda env 'gea' not found, using current env"

# Log dirs use absolute paths so subprocesses that change directories still resolve correctly
LOG_DIR="$BASE_DIR/3_POST_PROC/logs/log_files"
TIME_DIR="$BASE_DIR/3_POST_PROC/logs/time_logs"
SPACE_DIR="$BASE_DIR/3_POST_PROC/logs/space_logs"
SPACE_TIME_DIR="$BASE_DIR/3_POST_PROC/logs/space_time_logs"
ERROR_WARN_DIR="$BASE_DIR/3_POST_PROC/logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="$BASE_DIR/3_POST_PROC/logs/software_catalogs"
GPU_LOG_DIR="$BASE_DIR/3_POST_PROC/logs/gpu_log"
export LOG_DIR TIME_DIR SPACE_DIR SPACE_TIME_DIR ERROR_WARN_DIR SOFTWARE_CATALOG_DIR GPU_LOG_DIR

setup_logging "$CLEAR_LOGS"
export LOG_FILE TIME_FILE SPACE_FILE SPACE_TIME_FILE ERROR_WARN_FILE SOFTWARE_FILE GPU_LOG_FILE

# Skip software catalog if already generated this session (saves ~2s)
if [[ ! -f "${SOFTWARE_FILE:-}" ]] || [[ ! -s "${SOFTWARE_FILE:-}" ]]; then
	catalog_all_software
else
	log_info "Software catalog already exists, skipping: $SOFTWARE_FILE"
fi

log_step "Starting Post-Processing Pipeline (${#PIPELINE_CONFIGS[@]} config(s) enabled)"

# Cache parallel availability once (avoids 3 PATH lookups per dataset iteration)
_HAS_PARALLEL=false
command -v parallel &>/dev/null && _HAS_PARALLEL=true

#===============================================================================
# CONFIG LOOP
#===============================================================================

for CONFIG_FILE in "${PIPELINE_CONFIGS[@]}"; do

    [[ "$CONFIG_FILE" != /* ]] && CONFIG_FILE="$SCRIPT_DIR/$CONFIG_FILE"
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_warn "Config file not found, skipping: $CONFIG_FILE"
        continue
    fi

    # Load config (sets THREADS, METHODS, ANALYSES, GENE_GROUPS, SRR_DATASETS, etc.)
    source "$CONFIG_FILE"
    log_step "Config: $(basename "$CONFIG_FILE")"

    # Snapshot error/warning line count so we can report per-config delta
    _err_baseline=0
    [[ -f "$ERROR_WARN_FILE" ]] && _err_baseline=$(wc -l < "$ERROR_WARN_FILE")

    MASTER_REFERENCE="${MASTER_REFERENCES[0]}"
    [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]] && OVERWRITE_EXISTING="TRUE" || OVERWRITE_EXISTING="FALSE"
    export AVAILABLE_RAM_GB GPU_VRAM_GB OVERWRITE_EXISTING

    # Build combined SRR list and cache per-dataset results (avoids re-parsing later)
    SRR_COMBINED_LIST=()
    declare -A _CACHED_SRR_LISTS=()
    for dataset in "${SRR_DATASETS[@]}"; do
        csv_file="$SRR_CSV_DIR/${dataset}.csv"
        if [[ -f "$csv_file" ]]; then
            _cached=$(parse_srr_csv "$csv_file")
            _CACHED_SRR_LISTS["$dataset"]="$_cached"
            mapfile -t -O "${#SRR_COMBINED_LIST[@]}" SRR_COMBINED_LIST <<< "$_cached"
        else
            log_warn "SRR CSV not found: $csv_file"
        fi
    done

    # Validate
    [[ -z "$MASTER_REFERENCE" ]]          && { log_error "No master reference in $CONFIG_FILE"; continue; }
    [[ ${#METHODS[@]} -eq 0 ]]            && { log_error "No methods in $CONFIG_FILE"; continue; }
    [[ ${#GENE_GROUPS[@]} -eq 0 ]]        && { log_error "No gene groups in $CONFIG_FILE"; continue; }
    [[ ${#SRR_COMBINED_LIST[@]} -eq 0 ]]  && { log_error "No SRR samples loaded from $CONFIG_FILE"; continue; }

    # Parallel job count
    JOBS=1
    if [[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]]; then
        JOBS=$((THREADS / DESIRED_CPU_PER_JOB))
        (( JOBS < 1 )) && JOBS=1
    fi

    # Clear output folders if requested (batch find+delete instead of per-folder rm)
    if [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]]; then
        log_info "Clearing output folders for $MASTER_REFERENCE..."
        _cleared=0
        for method in "${METHODS[@]}"; do
            output_base="$BASE_DIR/3_POST_PROC/$method/Figure_Outputs"
            [[ -d "$output_base" ]] || continue
            for analysis in "${ANALYSES[@]}"; do
                folder_name=$(get_output_folder_name "$analysis")
                target="$output_base/$folder_name/$MASTER_REFERENCE"
                if [[ -n "$folder_name" && -d "$target" ]]; then
                    find "$target" -mindepth 1 -delete 2>/dev/null || true
                    ((_cleared++)) || true
                fi
            done
        done
        [[ $_cleared -gt 0 ]] && log_info "  Cleared $_cleared output directories"
    fi

    # Export for R scripts and subprocesses
    export BASE_DIR THREADS ENABLE_GPU ANALYSIS_MODULES_DIR GENE_GROUPS_DIR UTILITIES_DIR SRR_CSV_DIR MASTER_REFERENCE
    export GENE_GROUPS_STR="${GENE_GROUPS[*]}" ANALYSES_STR="${ANALYSES[*]}" SRR_DATASETS_STR="${SRR_DATASETS[*]}"

    # Run summary
    log_info "Threads: $THREADS | GPU: $ENABLE_GPU (${GPU_VRAM_GB}GB VRAM) | Parallel: $ENABLE_GNU_PARALLEL (Jobs: $JOBS)"
    log_info "Reference: $MASTER_REFERENCE | Methods: ${METHODS[*]}"
    log_info "Gene Groups: ${GENE_GROUPS[*]} | Datasets: ${SRR_DATASETS[*]} (${#SRR_COMBINED_LIST[@]} samples)"
    log_info "Analyses: ${ANALYSES[*]}"

    # Process each dataset (reuse cached CSV parse — avoids redundant file I/O)
    for dataset in "${SRR_DATASETS[@]}"; do
        if [[ -z "${_CACHED_SRR_LISTS[$dataset]+x}" ]]; then
            log_warn "No cached SRR data for $dataset, skipping"
            continue
        fi

        mapfile -t CURRENT_SRR_LIST <<< "${_CACHED_SRR_LISTS[$dataset]}"
        if [[ ${#CURRENT_SRR_LIST[@]} -eq 0 ]]; then
            log_warn "No samples found in $dataset, skipping"
            continue
        fi

        export CURRENT_DATASET="$dataset" SRR_COMBINED_LIST_STR="${CURRENT_SRR_LIST[*]}"
        log_step "Dataset: $dataset (${#CURRENT_SRR_LIST[@]} samples)"

        # ── Phase 0: Setup method environments ──
        for method in "${METHODS[@]}"; do
            setup_method_env "$method" "$MASTER_REFERENCE" || {
                log_error "Failed to set up environment for $method — skipping dataset $dataset"
                continue 2
            }
        done

        # Export utils once per dataset (avoids repeated export -f in each phase)
        if $_HAS_PARALLEL && [[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]]; then
            export_utils_for_parallel
            export SCRIPT_DIR LOG_FILE ERROR_WARN_FILE RUN_ID
        fi

        # ── Phase 1: Preprocessing (parallel across methods when possible) ──
        if [[ ${#METHODS[@]} -gt 1 && "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && $_HAS_PARALLEL; then
            log_info "Phase 1: Preprocessing (parallel across ${#METHODS[@]} methods)"
            printf '%s\n' "${METHODS[@]}" | parallel \
                -j "${#METHODS[@]}" \
                --halt soon,fail=1 \
                --joblog "$LOG_DIR/parallel_preproc_${dataset}.log" \
                run_method_preprocessing {} "$MASTER_REFERENCE"
        else
            log_info "Phase 1: Preprocessing (sequential, all threads)"
            for method in "${METHODS[@]}"; do
                run_method_preprocessing "$method" "$MASTER_REFERENCE"
            done
        fi

        # ── Phase 2: Thread-heavy analyses (sequential — needs full threads) ──
        # Includes: Matrix_Creation, Differential_Expression, WGCNA, GSEA, PCA
        HEAVY_ANALYSES=()
        FIGURE_ANALYSES=()
        for analysis in "${ANALYSES[@]}"; do
            [[ -z "$analysis" ]] && continue
            if is_figure_analysis "$analysis"; then
                FIGURE_ANALYSES+=("$analysis")
            else
                HEAVY_ANALYSES+=("$analysis")
            fi
        done

        if [[ ${#HEAVY_ANALYSES[@]} -gt 0 ]]; then
            if [[ ${#METHODS[@]} -gt 1 && "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && $_HAS_PARALLEL; then
                # Each method's heavy analyses are independent of other methods.
                # Run methods in parallel, each getting THREADS/N_METHODS cores.
                _n_methods=${#METHODS[@]}
                _threads_per_method=$(( THREADS / _n_methods ))
                (( _threads_per_method < 1 )) && _threads_per_method=1
                log_info "Phase 2: Heavy analyses (parallel across ${_n_methods} methods, ${_threads_per_method} threads each): ${HEAVY_ANALYSES[*]}"
                # Worker function: run all heavy analyses for one method with reduced thread count
                _run_heavy_for_method() {
                    local _method="$1" _ref="$2" _orig_threads="$THREADS"
                    shift 2
                    export THREADS="$_threads_per_method"
                    for _analysis in "$@"; do
                        run_single_analysis "$_method" "$_ref" "$_analysis"
                    done
                    export THREADS="$_orig_threads"
                }
                export -f _run_heavy_for_method
                export _threads_per_method
                printf '%s\n' "${METHODS[@]}" | parallel \
                    -j "$_n_methods" \
                    --halt soon,fail=1 \
                    --joblog "$LOG_DIR/parallel_heavy_${dataset}.log" \
                    _run_heavy_for_method {} "$MASTER_REFERENCE" "${HEAVY_ANALYSES[@]}"
            else
                log_info "Phase 2: Heavy analyses (sequential): ${HEAVY_ANALYSES[*]}"
                for method in "${METHODS[@]}"; do
                    log_step "Processing Method (heavy): $method"
                    for analysis in "${HEAVY_ANALYSES[@]}"; do
                        run_single_analysis "$method" "$MASTER_REFERENCE" "$analysis"
                    done
                done
            fi
        fi

        # ── Phase 3: Figure generation (parallelisable — lightweight per job) ──
        if [[ ${#FIGURE_ANALYSES[@]} -gt 0 ]]; then
            if [[ "$ENABLE_GNU_PARALLEL" == "TRUE" && $JOBS -gt 1 ]] && $_HAS_PARALLEL; then
                log_info "Phase 3: Figure generation (GNU Parallel, $JOBS jobs): ${FIGURE_ANALYSES[*]}"

                # Build method×analysis pairs (tab-separated) and parallelise
                PARALLEL_TASKS=()
                for method in "${METHODS[@]}"; do
                    for analysis in "${FIGURE_ANALYSES[@]}"; do
                        PARALLEL_TASKS+=("${method}"$'\t'"${analysis}")
                    done
                done

                printf '%s\n' "${PARALLEL_TASKS[@]}" | parallel \
                    -j "$JOBS" \
                    --colsep '\t' \
                    --halt soon,fail=1 \
                    --joblog "$LOG_DIR/parallel_figures_${dataset}.log" \
                    run_single_analysis {1} "$MASTER_REFERENCE" {2}
            else
                log_info "Phase 3: Figure generation (sequential): ${FIGURE_ANALYSES[*]}"
                for method in "${METHODS[@]}"; do
                    for analysis in "${FIGURE_ANALYSES[@]}"; do
                        run_single_analysis "$method" "$MASTER_REFERENCE" "$analysis"
                    done
                done
            fi
        fi
    done

    # Config summary
    log_step "Config Complete: $(basename "$CONFIG_FILE")"
    log_info "Datasets: ${SRR_DATASETS[*]} | Methods: ${#METHODS[@]} per dataset"
    _err_total=0
    [[ -f "$ERROR_WARN_FILE" ]] && _err_total=$(wc -l < "$ERROR_WARN_FILE")
    _err_delta=$(( _err_total - _err_baseline ))
    if [[ $_err_delta -gt 0 ]]; then
        log_info "Errors/Warnings this config: $_err_delta (total: $_err_total, see $ERROR_WARN_FILE)"
    else
        log_info "No errors encountered for this config"
    fi

done

#===============================================================================
# FINAL SUMMARY
#===============================================================================

log_step "All Configs Complete"
log_info "Configs run: ${#PIPELINE_CONFIGS[@]} | Log: $LOG_FILE | Time: $TIME_FILE"
log_step "Pipeline completed at $(date)"
