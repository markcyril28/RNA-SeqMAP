#!/bin/bash
#===============================================================================
# MASTER POST-PROCESSING SCRIPT - RNA-SEQ ANALYSIS PIPELINE
#===============================================================================

set -o pipefail   # -e/-u omitted intentionally (sourced functions use boolean returns)

# ==============================================================================
# SYSTEM RESOURCES
# ==============================================================================

# Respect pre-set THREADS from environment; only probe nproc if unset
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 12)}"
ENABLE_GPU="TRUE"
ENABLE_GNU_PARALLEL="TRUE"
DESIRED_CPU_PER_JOB=1

# Number of pipeline configs to process concurrently.
# "auto" = RAM-aware auto-detection (recommended for HPC); 1 = sequential (safest).
# TRADE-OFF: parallel configs don't share SRR CSV cache (~1ms overhead per config),
# but reduce total wall-clock by up to Nx for independent method/reference configs.
# Big O: reduces O(C × T_per_config) to O(ceil(C/P) × T_per_config)
# where C=configs, P=PARALLEL_CONFIGS, T_per_config=time per config.
PARALLEL_CONFIGS="${PARALLEL_CONFIGS:-auto}"

# Auto-detect available RAM (fallback: 24 GB)
# Pure bash: avoids awk fork. O(1) — reads ~25 lines then breaks.
# Matches run_concordance_analysis.sh pattern (Pass 30).
if [[ -f /proc/meminfo ]]; then
    while IFS=' ' read -r _key _val _; do
        if [[ "$_key" == "MemAvailable:" ]]; then
            AVAILABLE_RAM_GB=$(( _val / 1048576 ))
            break
        fi
    done < /proc/meminfo
    unset _key _val
elif command -v sysctl &>/dev/null; then
    _raw_bytes=$(sysctl -n hw.memsize 2>/dev/null)
    AVAILABLE_RAM_GB=$(( _raw_bytes / 1073741824 ))
    unset _raw_bytes
fi
# Fallback covers both unset and empty (e.g., MemAvailable line missing from /proc/meminfo)
AVAILABLE_RAM_GB="${AVAILABLE_RAM_GB:-24}"
# Ensure numeric (strip non-digits) before arithmetic
[[ "$AVAILABLE_RAM_GB" =~ ^[0-9]+$ ]] || AVAILABLE_RAM_GB=24
# Floor: ensure at least 4 GB to avoid starving R scripts
(( AVAILABLE_RAM_GB < 4 )) && AVAILABLE_RAM_GB=4

GPU_VRAM_GB=8

# ==============================================================================
# LOGGING AND OUTPUT
# ==============================================================================

CLEAR_LOGS="TRUE"
CLEAR_OUTPUT_FOLDER="TRUE"

# Figure resolution in DPI (300–600)
FIGURE_DPI="${FIGURE_DPI:-300}"

# HTML viewer: auto-generate an interactive results viewer in 3_POST_PROC/
# Set to "TRUE" to generate alignment_results_viewer.html after all configs run.
GENERATE_HTML_VIEWER="${GENERATE_HTML_VIEWER:-TRUE}"

# Single-config child mode: when invoked by parallel dispatch, process only one config.
# The parent sets __PP_SINGLE_CONFIG to the config path and re-invokes this script.
if [[ -n "${__PP_SINGLE_CONFIG:-}" ]]; then
    PIPELINE_CONFIGS=("$__PP_SINGLE_CONFIG")
    CLEAR_LOGS="FALSE"      # Don't clear shared logs from parent
    PARALLEL_CONFIGS=1       # Don't recurse into parallel dispatch
    # Unique RUN_ID per child to avoid log file contention
    _cfg_basename="${__PP_SINGLE_CONFIG##*/}"
    _cfg_basename="${_cfg_basename%.toml}"
    printf -v RUN_ID '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || RUN_ID=$(date +%Y%m%d_%H%M%S)
    RUN_ID="${RUN_ID}_${_cfg_basename}"
    unset _cfg_basename
    # Clear exported sentinel variables inherited from parent so modules re-source
    # their function definitions in this new bash process. (The parent exports these
    # sentinels for GNU Parallel subshells, but child processes started via
    # `bash "$0"` need to re-source the files since function definitions don't
    # survive across process boundaries.)
    unset PIPELINE_UTILS_SOURCED LOGGING_UTILS_SOURCED _TOML_PARSER_SOURCED LOGGING_INITIALIZED
else

PIPELINE_CONFIGS=(
    # ── Full — Eggplant_V4.1 ──
    "config/3_post_proc_configs/HPC_full_M1_Eggplant_V4.1.toml"    # M1 HISAT2 RefGuided   Eggplant_V4.1 genome
    "config/3_post_proc_configs/HPC_full_M2_Eggplant_V4.1.toml"    # M2 HISAT2 DeNovo      Eggplant_V4.1 transcript
    "config/3_post_proc_configs/HPC_full_M3_Eggplant_V4.1.toml"    # M3 STAR Align         Eggplant_V4.1 genome
    "config/3_post_proc_configs/HPC_full_M4_Eggplant_V4.1.toml"    # M4 Salmon SAF         Eggplant_V4.1 transcript
    "config/3_post_proc_configs/HPC_full_M5_Eggplant_V4.1.toml"    # M5 RSEM Bowtie2       Eggplant_V4.1 transcript

    # ── Full — GPE001970 ──
    "config/3_post_proc_configs/HPC_full_M1_GPE001970.toml"         # M1 HISAT2 RefGuided   GPE001970 genome
    "config/3_post_proc_configs/HPC_full_M2_GPE001970.toml"         # M2 HISAT2 DeNovo      GPE001970 transcript
    "config/3_post_proc_configs/HPC_full_M3_GPE001970.toml"         # M3 STAR Align         GPE001970 genome
    "config/3_post_proc_configs/HPC_full_M4_GPE001970.toml"         # M4 Salmon SAF         GPE001970 transcript
    "config/3_post_proc_configs/HPC_full_M5_GPE001970.toml"         # M5 RSEM Bowtie2       GPE001970 transcript
)

fi  # end of single-config child mode check

#===============================================================================
# PATHS AND UTILITIES
#===============================================================================

# Resolve script directory: parameter expansion avoids nested $(dirname) subshell
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"
BASE_DIR="$SCRIPT_DIR"
ANALYSIS_MODULES_DIR="$BASE_DIR/modules/c_post_processing/analysis_modules"
GENE_GROUPS_DIR="$BASE_DIR/inputs/3_post_proc_inputs/gene_groups_csv"
SRR_CSV_DIR="$BASE_DIR/inputs/3_post_proc_inputs/SRR_csv"
UTILITIES_DIR="$BASE_DIR/modules/c_post_processing/utilities"

source "$BASE_DIR/modules/logging/logging_utils.sh"
source "$UTILITIES_DIR/pipeline_utils.sh"
source "$BASE_DIR/config/shared/toml_parser.sh"

if [[ ! -d "$ANALYSIS_MODULES_DIR" ]]; then
    log_warn "Analysis modules directory not found: $ANALYSIS_MODULES_DIR"
    log_warn "R analysis scripts (Matrix_Creation, Basic_Heatmap, Heatmap_with_CV) will not run."
    log_warn "Restore from z_archive/modules/c_post_processing/analysis_modules/ if needed."
fi


#===============================================================================
# FUNCTIONS
#===============================================================================

# Map analysis name to output folder name.
# Uses case statement (not associative array) so the function works correctly
# when exported to GNU Parallel subshells (bash cannot export associative arrays).
get_output_folder_name() {
    case "$1" in
        "Matrix_Creation")                echo "0_Matrix_Creation" ;;
        "Basic_Heatmap")                  echo "I_Basic_Heatmap" ;;
        "Heatmap_with_CV")               echo "II_Heatmap_with_CV" ;;
        *)                                echo "" ;;
    esac
}

#===============================================================================
# INITIALIZATION
#===============================================================================

[[ ${#PIPELINE_CONFIGS[@]} -eq 0 ]] && { log_error "No configs enabled in PIPELINE_CONFIGS"; exit 1; }

# Skip conda hook (~0.3-0.5s) if already in the correct environment
if [[ "${CONDA_DEFAULT_ENV:-}" != "gea" ]]; then
    eval "$(conda shell.bash hook 2>/dev/null)" 2>/dev/null || true
    conda activate gea 2>/dev/null || log_warn "conda env 'gea' not found, using current env"
fi

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

log_step "Starting Post-Processing Pipeline (${#PIPELINE_CONFIGS[@]} config(s), parallel=$PARALLEL_CONFIGS)"

# Cache parallel availability once (avoids 3 PATH lookups per dataset iteration)
_HAS_PARALLEL=false
command -v parallel &>/dev/null && _HAS_PARALLEL=true

# Pre-compute parallel job count (THREADS and DESIRED_CPU_PER_JOB are set once at top)
# Two limits: CPU-based and RAM-based (each R process loads ggplot2/ComplexHeatmap ≈ 800MB)
JOBS=1
if [[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]]; then
    JOBS=$((THREADS / DESIRED_CPU_PER_JOB))
    (( JOBS < 1 )) && JOBS=1
    # Memory guard: cap concurrent R jobs so total < 75% of RAM
    # Each R figure-generation process uses ~800MB (ggplot2 + ComplexHeatmap + data)
    _R_MEM_MB=800
    _MAX_JOBS_BY_RAM=$(( AVAILABLE_RAM_GB * 1024 * 75 / 100 / _R_MEM_MB ))
    (( _MAX_JOBS_BY_RAM < 1 )) && _MAX_JOBS_BY_RAM=1
    if (( JOBS > _MAX_JOBS_BY_RAM )); then
        log_info "Capping parallel jobs from $JOBS to $_MAX_JOBS_BY_RAM (RAM limit: ${AVAILABLE_RAM_GB}GB, ~${_R_MEM_MB}MB/job)"
        JOBS=$_MAX_JOBS_BY_RAM
    fi
    unset _R_MEM_MB _MAX_JOBS_BY_RAM
fi

# Auto-detect parallel config count based on available RAM and CPU.
# Each config's R session uses ~2GB (ComplexHeatmap + ggplot2 + data).
if [[ "$PARALLEL_CONFIGS" == "auto" ]]; then
    _CFG_MEM_MB=2048
    PARALLEL_CONFIGS=$(( AVAILABLE_RAM_GB * 1024 * 60 / 100 / _CFG_MEM_MB ))
    (( PARALLEL_CONFIGS < 1 )) && PARALLEL_CONFIGS=1
    # Cap at number of configs
    (( PARALLEL_CONFIGS > ${#PIPELINE_CONFIGS[@]} )) && PARALLEL_CONFIGS=${#PIPELINE_CONFIGS[@]}
    # Cap at half of CPU threads (each config runs its own R jobs)
    _max_by_cpu=$(( THREADS / 2 ))
    (( _max_by_cpu < 1 )) && _max_by_cpu=1
    (( PARALLEL_CONFIGS > _max_by_cpu )) && PARALLEL_CONFIGS=$_max_by_cpu
    unset _CFG_MEM_MB _max_by_cpu
fi
# Ensure numeric
[[ "$PARALLEL_CONFIGS" =~ ^[0-9]+$ ]] || PARALLEL_CONFIGS=1

#===============================================================================
# PARALLEL CONFIG DISPATCH
#===============================================================================
# When PARALLEL_CONFIGS > 1 and multiple configs exist, dispatch each config
# as a child process for concurrent execution. Each child re-invokes this script
# with __PP_SINGLE_CONFIG set to process exactly one config.
# All functions and environment are inherited via the child bash process.
# Big O: O(ceil(C/P) × T_config) wall-clock instead of O(C × T_config).

_TOTAL_CONFIGS=${#PIPELINE_CONFIGS[@]}

if [[ "$PARALLEL_CONFIGS" -gt 1 && ${#PIPELINE_CONFIGS[@]} -gt 1 ]]; then
    log_step "Parallel Config Dispatch (${#PIPELINE_CONFIGS[@]} configs, max $PARALLEL_CONFIGS concurrent)"

    declare -a _cfg_pids=() _cfg_logs=()
    _cfg_active=0

    for _cfg_file in "${PIPELINE_CONFIGS[@]}"; do
        # Throttle: wait for a slot when at concurrency limit
        while (( _cfg_active >= PARALLEL_CONFIGS )); do
            if wait -n 2>/dev/null; then
                (( _cfg_active-- ))
            else
                # Fallback for bash < 4.3: wait for all and reset counter
                wait; _cfg_active=0
            fi
        done

        _cfg_tag="${_cfg_file##*/}"
        _cfg_tag="${_cfg_tag%.toml}"
        _cfg_log="${LOG_DIR}/parallel_config_${_cfg_tag}_${RUN_ID}.log"
        _cfg_logs+=("$_cfg_log")

        log_info "  Dispatching: ${_cfg_file##*/}"
        __PP_SINGLE_CONFIG="$_cfg_file" bash "$0" > "$_cfg_log" 2>&1 &
        _cfg_pids+=($!)
        (( _cfg_active++ ))
    done

    # Wait for all configs and collect exit statuses
    _cfg_failures=0
    for _ci in "${!_cfg_pids[@]}"; do
        if ! wait "${_cfg_pids[$_ci]}"; then
            log_warn "Config ${PIPELINE_CONFIGS[$_ci]##*/} failed (pid=${_cfg_pids[$_ci]})"
            (( _cfg_failures++ )) || true
        fi
    done

    # Aggregate per-config logs into main log file (single I/O operation)
    _existing_logs=()
    for _cl in "${_cfg_logs[@]}"; do
        [[ -f "$_cl" ]] && _existing_logs+=("$_cl")
    done
    if [[ ${#_existing_logs[@]} -gt 0 ]]; then
        cat "${_existing_logs[@]}" >> "$LOG_FILE"
        rm -f "${_existing_logs[@]}"
    fi

    if [[ $_cfg_failures -gt 0 ]]; then
        log_warn "$_cfg_failures of ${#PIPELINE_CONFIGS[@]} configs had failures"
    fi
    log_step "Parallel dispatch complete"

    # Skip the sequential config loop below (all configs already processed)
    PIPELINE_CONFIGS=()
fi

#===============================================================================
# CONFIG LOOP
#===============================================================================

# Global CSV cache persists across configs — avoids re-parsing the same SRR CSV
# when multiple configs reference the same datasets. O(datasets) instead of O(configs × datasets).
declare -A _GLOBAL_SRR_CACHE=()

# Export static variables once before the config loop — these never change between configs.
# Per-config variables (MASTER_REFERENCE, GENE_GROUPS_STR, etc.) are exported inside the loop.
export BASE_DIR THREADS ENABLE_GPU ANALYSIS_MODULES_DIR GENE_GROUPS_DIR UTILITIES_DIR SRR_CSV_DIR
export AVAILABLE_RAM_GB GPU_VRAM_GB FIGURE_DPI

# Static analysis->folder mapping — computed once, reused across all configs.
# Avoids re-declaring the associative array inside the CLEAR_OUTPUT_FOLDER block per config.
declare -A _FOLDER_NAME_MAP=(
    [Matrix_Creation]="0_Matrix_Creation"
    [Basic_Heatmap]="I_Basic_Heatmap"
    [Heatmap_with_CV]="II_Heatmap_with_CV"
)

for CONFIG_FILE in "${PIPELINE_CONFIGS[@]}"; do

    [[ "$CONFIG_FILE" != /* ]] && CONFIG_FILE="$SCRIPT_DIR/$CONFIG_FILE"
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_warn "Config file not found, skipping: $CONFIG_FILE"
        continue
    fi

    # Load config (sets METHODS, ANALYSES, GENE_GROUPS, SRR_DATASETS, etc.)
    # TOML keys are parsed as uppercase bash variables by load_toml
    load_toml "$CONFIG_FILE"
    log_step "Config: ${CONFIG_FILE##*/}"

    # Snapshot error/warning line count so we can report per-config delta.
    # Single wc -l fork: O(1) C-level counting vs O(L) bash read loop.
    _err_baseline=0
    if [[ -f "$ERROR_WARN_FILE" && -s "$ERROR_WARN_FILE" ]]; then
        _err_baseline=$(wc -l < "$ERROR_WARN_FILE")
    fi

    MASTER_REFERENCE="${MASTER_REFERENCES[0]}"
    if [[ ${#MASTER_REFERENCES[@]} -gt 1 ]]; then
        log_warn "MASTER_REFERENCES has ${#MASTER_REFERENCES[@]} entries but only the first ('$MASTER_REFERENCE') is used."
        log_warn "Use separate per-reference configs to process multiple references."
    fi
    [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]] && OVERWRITE_EXISTING="TRUE" || OVERWRITE_EXISTING="FALSE"
    export OVERWRITE_EXISTING

    # Build combined SRR list using global cache — avoids re-parsing CSVs across configs
    SRR_COMBINED_LIST=()
    declare -A _CACHED_SRR_LISTS=()
    for dataset in "${SRR_DATASETS[@]}"; do
        # Check global cross-config cache first
        if [[ -n "${_GLOBAL_SRR_CACHE[$dataset]+set}" ]]; then
            _cached="${_GLOBAL_SRR_CACHE[$dataset]}"
        else
            csv_file="$SRR_CSV_DIR/${dataset}.csv"
            if [[ -f "$csv_file" ]]; then
                _cached=$(parse_srr_csv "$csv_file")
                [[ -n "$_cached" ]] && _GLOBAL_SRR_CACHE["$dataset"]="$_cached"
            else
                log_warn "SRR CSV not found: $csv_file"
                _cached=""
            fi
        fi
        if [[ -n "$_cached" ]]; then
            _CACHED_SRR_LISTS["$dataset"]="$_cached"
            mapfile -t -O "${#SRR_COMBINED_LIST[@]}" SRR_COMBINED_LIST <<< "$_cached"
        else
            log_warn "No valid SRR entries for dataset: $dataset"
        fi
    done

    # Validate
    [[ -z "$MASTER_REFERENCE" ]]          && { log_error "No master reference in $CONFIG_FILE"; continue; }
    [[ ${#METHODS[@]} -eq 0 ]]            && { log_error "No methods in $CONFIG_FILE"; continue; }
    [[ ${#GENE_GROUPS[@]} -eq 0 ]]        && { log_error "No gene groups in $CONFIG_FILE"; continue; }
    [[ ${#SRR_COMBINED_LIST[@]} -eq 0 ]]  && { log_error "No SRR samples loaded from $CONFIG_FILE"; continue; }
    [[ ${#ANALYSES[@]} -eq 0 ]]           && { log_error "No analyses in $CONFIG_FILE"; continue; }

    # Clear output folders if requested (rm -rf + mkdir is faster than find -delete on deep trees)
    if [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]]; then
        log_info "Clearing output folders for $MASTER_REFERENCE..."
        # Uses _FOLDER_NAME_MAP declared once before the config loop (avoids per-config re-declaration)
        _clear_targets=()
        for method in "${METHODS[@]}"; do
            output_base="$BASE_DIR/3_POST_PROC/$method/Figure_Outputs"
            [[ -d "$output_base" ]] || continue
            for analysis in "${ANALYSES[@]}"; do
                folder_name="${_FOLDER_NAME_MAP[$analysis]}"
                target="$output_base/$folder_name/$MASTER_REFERENCE"
                [[ -n "$folder_name" && -d "$target" ]] && _clear_targets+=("$target")
            done
        done
        if [[ ${#_clear_targets[@]} -gt 0 ]]; then
            rm -rf "${_clear_targets[@]}"
            mkdir -p "${_clear_targets[@]}"
            log_info "  Cleared ${#_clear_targets[@]} output directories"
        fi
    fi

    # Export per-config variables for R scripts and subprocesses
    # (static vars like BASE_DIR, THREADS, etc. are exported once before the loop)
    export MASTER_REFERENCE
    export GENE_GROUPS_STR="${GENE_GROUPS[*]}" ANALYSES_STR="${ANALYSES[*]}" SRR_DATASETS_STR="${SRR_DATASETS[*]}"

    # Run summary
    log_info "Threads: $THREADS | GPU: $ENABLE_GPU (${GPU_VRAM_GB}GB VRAM) | Parallel: $ENABLE_GNU_PARALLEL (Jobs: $JOBS)"
    log_info "Reference: $MASTER_REFERENCE | Methods: ${METHODS[*]}"
    log_info "Gene Groups: ${GENE_GROUPS[*]} | Datasets: ${SRR_DATASETS[*]} (${#SRR_COMBINED_LIST[@]} samples)"
    log_info "Analyses: ${ANALYSES[*]}"

    # ── Phase 0: Setup method environments (once per config, not per dataset) ──
    # Method environments depend on METHODS[] and MASTER_REFERENCE which are
    # config-level, not dataset-level. Hoisting saves N_datasets × N_methods
    # redundant setup calls.
    _methods_setup_ok=true
    for method in "${METHODS[@]}"; do
        setup_method_env "$method" "$MASTER_REFERENCE" || {
            log_error "Failed to set up environment for $method — skipping config"
            _methods_setup_ok=false
            break
        }
    done
    $_methods_setup_ok || continue

    # Export utils once per config (avoids repeated export -f per dataset)
    if $_HAS_PARALLEL && [[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]]; then
        export_utils_for_parallel
        export SCRIPT_DIR LOG_FILE ERROR_WARN_FILE RUN_ID
    fi

    # ── Pre-compute analysis classification (config-level, not dataset-level) ──
    # Split ANALYSES into heavy (thread-bound) and figure (parallelizable) categories
    # once per config instead of re-classifying per dataset.
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

    # Pre-export figure analyses string and worker function once per config (not per dataset).
    # FIGURE_ANALYSES is config-level — re-exporting per dataset was redundant.
    if [[ ${#FIGURE_ANALYSES[@]} -gt 0 ]] && $_HAS_PARALLEL && [[ "$ENABLE_GNU_PARALLEL" == "TRUE" && $JOBS -gt 1 ]]; then
        export _FIGURE_ANALYSES_STR="${FIGURE_ANALYSES[*]}"
        _run_figures_for_method() {
            local -a _figs=()
            IFS=' ' read -ra _figs <<< "$_FIGURE_ANALYSES_STR"
            run_batched_analyses "$1" "$2" "${_figs[@]}"
        }
        export -f _run_figures_for_method
    fi

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

        # SRR_COMBINED_LIST_STR holds the current dataset's samples (not a cross-dataset merge)
        export CURRENT_DATASET="$dataset" SRR_COMBINED_LIST_STR="${CURRENT_SRR_LIST[*]}"
        log_step "Dataset: $dataset (${#CURRENT_SRR_LIST[@]} samples)"

        # ── Phase 1: Preprocessing (parallel across methods when possible) ──
        if [[ ${#METHODS[@]} -gt 1 && "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && $_HAS_PARALLEL; then
            log_info "Phase 1: Preprocessing (parallel across ${#METHODS[@]} methods)"
            parallel \
                -j "${#METHODS[@]}" \
                --halt soon,fail,30% \
                --joblog "$LOG_DIR/parallel_preproc_${dataset}.log" \
                run_method_preprocessing {} "$MASTER_REFERENCE" \
                < <(printf '%s\n' "${METHODS[@]}") \
                || log_warn "Phase 1: Some preprocessing tasks failed for dataset '$dataset' (see joblog)"
        else
            log_info "Phase 1: Preprocessing (sequential, all threads)"
            for method in "${METHODS[@]}"; do
                run_method_preprocessing "$method" "$MASTER_REFERENCE"
            done
        fi

        # ── Phase 2: Thread-heavy analyses (sequential — needs full threads) ──
        if [[ ${#HEAVY_ANALYSES[@]} -gt 0 ]]; then
            if [[ ${#METHODS[@]} -gt 1 && "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && $_HAS_PARALLEL; then
                # Each method's heavy analyses are independent of other methods.
                # Run methods in parallel, each getting THREADS/N_METHODS cores.
                # Memory guard: DESeq2/WGCNA/GSEA use ~2-4GB per R process;
                # cap concurrent methods so total < 75% of RAM
                _n_methods=${#METHODS[@]}
                _heavy_mem_per_method=3072  # ~3GB per heavy R analysis
                _max_by_ram=$(( AVAILABLE_RAM_GB * 1024 * 75 / 100 / _heavy_mem_per_method ))
                (( _max_by_ram < 1 )) && _max_by_ram=1
                (( _n_methods > _max_by_ram )) && _n_methods=$_max_by_ram
                _threads_per_method=$(( THREADS / _n_methods ))
                (( _threads_per_method < 1 )) && _threads_per_method=1
                log_info "Phase 2: Heavy analyses (parallel across ${_n_methods} methods, ${_threads_per_method} threads each): ${HEAVY_ANALYSES[*]}"
                # Worker function: run all heavy analyses for one method with reduced thread count
                _run_heavy_for_method() {
                    local _method="$1" _ref="$2"
                    shift 2
                    export THREADS="$_threads_per_method"
                    run_batched_analyses "$_method" "$_ref" "$@"
                }
                export -f _run_heavy_for_method
                export _threads_per_method
                parallel \
                    -j "$_n_methods" \
                    --halt soon,fail,30% \
                    --joblog "$LOG_DIR/parallel_heavy_${dataset}.log" \
                    _run_heavy_for_method {} "$MASTER_REFERENCE" "${HEAVY_ANALYSES[@]}" \
                    < <(printf '%s\n' "${METHODS[@]}") \
                    || log_warn "Phase 2: Some heavy analyses failed for dataset '$dataset' (see joblog)"
            else
                log_info "Phase 2: Heavy analyses (sequential, batched): ${HEAVY_ANALYSES[*]}"
                for method in "${METHODS[@]}"; do
                    log_step "Processing Method (heavy): $method"
                    run_batched_analyses "$method" "$MASTER_REFERENCE" "${HEAVY_ANALYSES[@]}"
                done
            fi
        fi

        # ── Phase 3: Figure generation (parallelisable — lightweight per job) ──
        # _FIGURE_ANALYSES_STR and _run_figures_for_method are pre-exported above
        # the dataset loop (config-level, not dataset-level).
        if [[ ${#FIGURE_ANALYSES[@]} -gt 0 ]]; then
            if [[ "$ENABLE_GNU_PARALLEL" == "TRUE" && $JOBS -gt 1 ]] && $_HAS_PARALLEL; then
                log_info "Phase 3: Figure generation (GNU Parallel, $JOBS jobs, batched per method): ${FIGURE_ANALYSES[*]}"
                parallel \
                    -j "$JOBS" \
                    --halt soon,fail,30% \
                    --joblog "$LOG_DIR/parallel_figures_${dataset}.log" \
                    _run_figures_for_method {} "$MASTER_REFERENCE" \
                    < <(printf '%s\n' "${METHODS[@]}") \
                    || log_warn "Phase 3: Some figure tasks failed for dataset '$dataset' (see joblog)"
            else
                log_info "Phase 3: Figure generation (sequential, batched): ${FIGURE_ANALYSES[*]}"
                for method in "${METHODS[@]}"; do
                    run_batched_analyses "$method" "$MASTER_REFERENCE" "${FIGURE_ANALYSES[@]}"
                done
            fi
        fi
    done

    # Config summary
    log_step "Config Complete: ${CONFIG_FILE##*/}"
    log_info "Datasets: ${SRR_DATASETS[*]} | Methods: ${#METHODS[@]} per dataset"
    # Single wc -l fork: O(1) C-level counting vs O(L) bash read loop.
    _err_total=0
    if [[ -f "$ERROR_WARN_FILE" && -s "$ERROR_WARN_FILE" ]]; then
        _err_total=$(wc -l < "$ERROR_WARN_FILE")
    fi
    _err_delta=$(( _err_total - _err_baseline ))
    if [[ $_err_delta -gt 0 ]]; then
        log_info "Errors/Warnings this config: $_err_delta (total: $_err_total, see $ERROR_WARN_FILE)"
    else
        log_info "No errors encountered for this config"
    fi

done

#===============================================================================
# HTML VIEWER GENERATION
#===============================================================================

if [[ "${GENERATE_HTML_VIEWER:-TRUE}" == "TRUE" && -z "${__PP_SINGLE_CONFIG:-}" ]]; then
    _viewer_script="$BASE_DIR/modules/c_post_processing/utilities/generate_html_viewer.py"
    _post_proc_dir="$BASE_DIR/3_POST_PROC"
    if command -v python3 &>/dev/null && [[ -f "$_viewer_script" ]]; then
        log_step "Generating HTML Results Viewer"
        _viewer_out=$(python3 "$_viewer_script" "$_post_proc_dir" 2>"$LOG_DIR/html_viewer_gen.log") \
            && log_info "HTML viewer: $_viewer_out" \
            || log_warn "HTML viewer generation failed — see $LOG_DIR/html_viewer_gen.log"
    else
        log_warn "HTML viewer skipped: python3 not found or viewer script missing ($_viewer_script)"
    fi
fi

#===============================================================================
# FINAL SUMMARY
#===============================================================================

log_step "All Configs Complete"
log_info "Configs run: ${_TOTAL_CONFIGS} | Log: $LOG_FILE | Time: $TIME_FILE"
# O(1) bash builtin — avoids $(date) subprocess fork
_final_ts=""; printf -v _final_ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _final_ts=$(date '+%Y-%m-%d %H:%M:%S')
log_step "Pipeline completed at $_final_ts"
