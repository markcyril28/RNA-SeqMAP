#!/bin/bash
#===============================================================================
# MASTER POST-PROCESSING SCRIPT - RNA-SEQ ANALYSIS PIPELINE
#===============================================================================

set -o pipefail   # -e/-u omitted intentionally (sourced functions use boolean returns)

# ==============================================================================
# PIPELINE CONFIGURATION
# ==============================================================================
PIPELINE_CONFIGS=(
    # ── Full — Eggplant_V4.1 ──
    #"config/3_post_proc_configs/HPC_full_M1_Eggplant_V4.1.toml"    # M1 HISAT2 RefGuided   Eggplant_V4.1 genome
    #"config/3_post_proc_configs/HPC_full_M2_Eggplant_V4.1.toml"    # M2 HISAT2 DeNovo      Eggplant_V4.1 transcript
    #"config/3_post_proc_configs/HPC_full_M3_Eggplant_V4.1.toml"    # M3 STAR Align         Eggplant_V4.1 genome
    #"config/3_post_proc_configs/HPC_full_M4_Eggplant_V4.1.toml"    # M4 Salmon SAF         Eggplant_V4.1 transcript
    #"config/3_post_proc_configs/HPC_full_M5_Eggplant_V4.1.toml"    # M5 RSEM Bowtie2       Eggplant_V4.1 transcript

    # ── Full — GPE001970 ──
    #"config/3_post_proc_configs/HPC_full_M1_GPE001970.toml"         # M1 HISAT2 RefGuided   GPE001970 genome
    #"config/3_post_proc_configs/HPC_full_M2_GPE001970.toml"         # M2 HISAT2 DeNovo      GPE001970 transcript
    "config/3_post_proc_configs/HPC_full_M3_GPE001970.toml"         # M3 STAR Align         GPE001970 genome
    #"config/3_post_proc_configs/HPC_full_M4_GPE001970.toml"         # M4 Salmon SAF         GPE001970 transcript
    #"config/3_post_proc_configs/HPC_full_M5_GPE001970.toml"         # M5 RSEM Bowtie2       GPE001970 transcript
)

# Script-level child PID array — must be declared before traps so cleanup can reap
# child processes spawned in parallel config dispatch (prevents orphans on SIGTERM/SIGINT)
declare -a _cfg_pids=()

# Clean up background logging processes on exit/signal (prevents zombie processes
# in Nextflow/Snakemake containers that would stall work-dir cleanup)
_postproc_cleanup() {
    local rc=$?
    # Kill any background child config processes to prevent orphans on HPC
    if [[ ${#_cfg_pids[@]} -gt 0 ]]; then
        kill "${_cfg_pids[@]}" 2>/dev/null || true
        wait "${_cfg_pids[@]}" 2>/dev/null || true
    fi
    type -t _logging_cleanup_bg &>/dev/null && _logging_cleanup_bg
    exit $rc
}
trap _postproc_cleanup EXIT
trap 'type -t _logging_cleanup_bg &>/dev/null && _logging_cleanup_bg; exit 143' TERM
trap 'type -t _logging_cleanup_bg &>/dev/null && _logging_cleanup_bg; exit 130' INT

# ==============================================================================
# SYSTEM RESOURCES
# ==============================================================================

# Respect pre-set THREADS from environment; check HPC scheduler vars before nproc
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-${PBS_NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 12)}}}"
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
    # macOS: hw.memsize is total RAM. Use vm_stat to estimate available (free+inactive).
    # Falls back to 75% of total if vm_stat parsing fails.
    _raw_bytes=$(sysctl -n hw.memsize 2>/dev/null)
    _total_gb=$(( _raw_bytes / 1073741824 ))
    # Page size: 4096 on Intel, 16384 on Apple Silicon — query dynamically
    _page_size=$(sysctl -n hw.pagesize 2>/dev/null)
    _page_size="${_page_size:-4096}"
    _vm_free_pages=$(vm_stat 2>/dev/null | awk '/Pages free|Pages inactive/ {gsub(/\./,"",$NF); s+=$NF} END {print s+0}')
    if [[ "${_vm_free_pages:-0}" -gt 0 ]]; then
        AVAILABLE_RAM_GB=$(( _vm_free_pages * _page_size / 1073741824 ))
    else
        AVAILABLE_RAM_GB=$(( _total_gb * 75 / 100 ))
    fi
    unset _raw_bytes _total_gb _vm_free_pages _page_size
fi
# Fallback covers both unset and empty (e.g., MemAvailable line missing from /proc/meminfo)
AVAILABLE_RAM_GB="${AVAILABLE_RAM_GB:-24}"
# Ensure numeric (reset to default if non-numeric) before arithmetic
[[ "$AVAILABLE_RAM_GB" =~ ^[0-9]+$ ]] || AVAILABLE_RAM_GB=24
# Floor: ensure at least 4 GB to avoid starving R scripts
(( AVAILABLE_RAM_GB < 4 )) && AVAILABLE_RAM_GB=4

GPU_VRAM_GB=8

# ==============================================================================
# LOGGING AND OUTPUT
# ==============================================================================

CLEAR_LOGS="${CLEAR_LOGS:-TRUE}"
CLEAR_OUTPUT_FOLDER="${CLEAR_OUTPUT_FOLDER:-TRUE}"
CLEAR_CACHE="${CLEAR_CACHE:-TRUE}"

# Figure resolution in DPI (300–600)
FIGURE_DPI="${FIGURE_DPI:-300}"

# HTML viewer: auto-generate an interactive results viewer in II_RESULTS/3_POST_PROC/
# Set to "TRUE" to generate alignment_results_viewer.html after all configs run.
GENERATE_HTML_VIEWER="${GENERATE_HTML_VIEWER:-FALSE}"

# Single-config child mode: when invoked by parallel dispatch, process only one config.
# The parent sets __PP_SINGLE_CONFIG to the config path and re-invokes this script.
if [[ -n "${__PP_SINGLE_CONFIG:-}" ]]; then
    PIPELINE_CONFIGS=("$__PP_SINGLE_CONFIG")
    CLEAR_LOGS="FALSE"      # Don't clear shared logs from parent
    CLEAR_CACHE="FALSE"      # Don't re-clear caches from child
    PARALLEL_CONFIGS=1       # Don't recurse into parallel dispatch
    # Unique RUN_ID per child to avoid log file contention
    _cfg_basename="${__PP_SINGLE_CONFIG##*/}"
    _cfg_basename="${_cfg_basename%.toml}"
    printf -v RUN_ID '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || RUN_ID=$(date +%Y%m%d_%H%M%S)
    RUN_ID="${RUN_ID}_${_cfg_basename}"
    unset _cfg_basename
    # Clear sentinel variables so modules re-source their function definitions
    # in this child process (function definitions don't cross process boundaries).
    unset PIPELINE_UTILS_SOURCED LOGGING_UTILS_SOURCED _TOML_PARSER_SOURCED LOGGING_INITIALIZED
else

PIPELINE_CONFIGS=(
    # ── Full — Eggplant_V4.1 ──
    #"config/3_post_proc_configs/HPC_full_M1_Eggplant_V4.1.toml"    # M1 HISAT2 RefGuided   Eggplant_V4.1 genome
    #"config/3_post_proc_configs/HPC_full_M2_Eggplant_V4.1.toml"    # M2 HISAT2 DeNovo      Eggplant_V4.1 transcript
    #"config/3_post_proc_configs/HPC_full_M3_Eggplant_V4.1.toml"    # M3 STAR Align         Eggplant_V4.1 genome
    #"config/3_post_proc_configs/HPC_full_M4_Eggplant_V4.1.toml"    # M4 Salmon SAF         Eggplant_V4.1 transcript
    #"config/3_post_proc_configs/HPC_full_M5_Eggplant_V4.1.toml"    # M5 RSEM Bowtie2       Eggplant_V4.1 transcript

    # ── Full — GPE001970 ──
    #"config/3_post_proc_configs/HPC_full_M1_GPE001970.toml"         # M1 HISAT2 RefGuided   GPE001970 genome
    #"config/3_post_proc_configs/HPC_full_M2_GPE001970.toml"         # M2 HISAT2 DeNovo      GPE001970 transcript
    "config/3_post_proc_configs/HPC_full_M3_GPE001970.toml"         # M3 STAR Align         GPE001970 genome
    #"config/3_post_proc_configs/HPC_full_M4_GPE001970.toml"         # M4 Salmon SAF         GPE001970 transcript
    #"config/3_post_proc_configs/HPC_full_M5_GPE001970.toml"         # M5 RSEM Bowtie2       GPE001970 transcript
)

fi  # end of single-config child mode check

#===============================================================================
# PATHS AND UTILITIES
#===============================================================================

# Resolve script directory: honour pre-set BASE_DIR from orchestrators (Nextflow/Snakemake)
if [[ -z "${BASE_DIR:-}" ]]; then
    SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
    [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
    SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] run_post_processing.sh: Failed to resolve script directory" >&2; exit 1; }
    BASE_DIR="$SCRIPT_DIR"
else
    # Validate BASE_DIR is absolute; resolve if relative (defensive against misconfigured orchestrators)
    if [[ "$BASE_DIR" != /* ]]; then
        BASE_DIR="$(cd "$BASE_DIR" 2>/dev/null && pwd)" || { echo "[ERROR] BASE_DIR is set but invalid: $BASE_DIR" >&2; exit 1; }
    fi
    SCRIPT_DIR="$BASE_DIR"
fi
_SELF_SCRIPT="$SCRIPT_DIR/${BASH_SOURCE[0]##*/}"
ANALYSIS_MODULES_DIR="${ANALYSIS_MODULES_DIR:-$BASE_DIR/modules_gea/c_post_processing/analysis_modules}"
GENE_GROUPS_DIR="${GENE_GROUPS_DIR:-$BASE_DIR/I_INPUTS/inputs/3_post_proc_inputs/gene_groups_csv}"
SRR_CSV_DIR="${SRR_CSV_DIR:-$BASE_DIR/I_INPUTS/inputs/3_post_proc_inputs/SRR_csv}"
UTILITIES_DIR="${UTILITIES_DIR:-$BASE_DIR/modules_gea/c_post_processing/utilities}"

source "$BASE_DIR/modules_gea/logging/logging_utils.sh" || {
    echo "[ERROR] Failed to source logging_utils.sh: $BASE_DIR/modules_gea/logging/logging_utils.sh" >&2
    exit 1
}
source "$UTILITIES_DIR/pipeline_utils.sh" || {
    echo "[ERROR] Failed to source pipeline_utils.sh: $UTILITIES_DIR/pipeline_utils.sh" >&2
    exit 1
}
source "$BASE_DIR/config/shared/toml_parser.sh" || {
    echo "[ERROR] Failed to source TOML parser: $BASE_DIR/config/shared/toml_parser.sh" >&2
    exit 1
}

if [[ ! -d "$ANALYSIS_MODULES_DIR" ]]; then
    log_warn "Analysis modules directory not found: $ANALYSIS_MODULES_DIR"
    log_warn "R analysis scripts (Matrix_Creation, Basic_Heatmap, Heatmap_with_CV) will not run."
    log_warn "Restore from z_archive/modules/c_post_processing/analysis_modules/ if needed."
fi


#===============================================================================
# FUNCTIONS
#===============================================================================

#===============================================================================
# INITIALIZATION
#===============================================================================

[[ ${#PIPELINE_CONFIGS[@]} -eq 0 ]] && { log_error "No configs enabled in PIPELINE_CONFIGS"; exit 1; }

# Absolutize relative config paths against BASE_DIR so the pipeline works when
# CWD differs from project root (e.g., Nextflow scratch workDir, Snakemake shadow).
for _i in "${!PIPELINE_CONFIGS[@]}"; do
    [[ "${PIPELINE_CONFIGS[$_i]}" != /* ]] && PIPELINE_CONFIGS[$_i]="${BASE_DIR}/${PIPELINE_CONFIGS[$_i]}"
done
unset _i

# Pre-extract primary gene group from the first config so log dirs are gene-group-scoped
# before setup_logging() is called. load_toml() is available (toml_parser.sh sourced above).
# Unset TOML vars afterward so the main config loop starts clean.
if [[ -z "${CURRENT_GENE_GROUP:-}" && ${#PIPELINE_CONFIGS[@]} -gt 0 ]]; then
    _pre_cfg="${PIPELINE_CONFIGS[0]}"
    if [[ -f "$_pre_cfg" ]]; then
        load_toml "$_pre_cfg"
        CURRENT_GENE_GROUP="${GENE_GROUPS[0]:-}"
        unset GENE_GROUPS METHODS ANALYSES SRR_DATASETS MASTER_REFERENCES
    fi
fi
CURRENT_GENE_GROUP="${CURRENT_GENE_GROUP:-_active}"
export CURRENT_GENE_GROUP

# Skip conda activation when orchestrator manages the environment (Nextflow/Snakemake)
# Skip conda hook (~0.3-0.5s) if already in the correct environment
if [[ -z "${WF_MANAGED_ENV:-}" && "${CONDA_DEFAULT_ENV:-}" != "gea" ]]; then
    eval "$(conda shell.bash hook 2>/dev/null)" 2>/dev/null || true
    conda activate gea 2>/dev/null || log_warn "conda env 'gea' not found, using current env"
fi

# Log dirs use absolute paths so subprocesses that change directories still resolve correctly
LOG_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/logs/log_files"
TIME_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/logs/time_logs"
SPACE_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/logs/space_logs"
SPACE_TIME_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/logs/space_time_logs"
ERROR_WARN_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/logs/software_catalogs"
GPU_LOG_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/logs/gpu_log"
export LOG_DIR TIME_DIR SPACE_DIR SPACE_TIME_DIR ERROR_WARN_DIR SOFTWARE_CATALOG_DIR GPU_LOG_DIR

# Ensure log directories exist before setup_logging (under Nextflow/Snakemake, output
# directories may not be pre-created by the orchestrator)
mkdir -p "$LOG_DIR" "$TIME_DIR" "$SPACE_DIR" "$SPACE_TIME_DIR" \
         "$ERROR_WARN_DIR" "$SOFTWARE_CATALOG_DIR" "$GPU_LOG_DIR" 2>/dev/null || true

# Mirror the full pipeline log to a top-level RNA-SeqMAP/logs/ directory so the
# canonical run log is reachable without descending into II_RESULTS/3_POST_PROC/.
# MIRROR_LOG_FILE is honoured by _logging_setup_redirect in logging_utils.sh.
MIRROR_LOG_DIR="$BASE_DIR/logs/log_files"
MIRROR_LOG_FILE="$MIRROR_LOG_DIR/pipeline_${RUN_ID}_full_log.log"
mkdir -p "$MIRROR_LOG_DIR" 2>/dev/null || true
if [[ "${CLEAR_LOGS^^}" == "TRUE" ]]; then
    # Clear root-level mirror logs too — setup_logging only clears under LOG_DIR
    find "$MIRROR_LOG_DIR" -maxdepth 1 -type f -name '*.log' -delete 2>/dev/null || true
fi
export MIRROR_LOG_FILE

setup_logging "$CLEAR_LOGS"
export LOG_FILE TIME_FILE SPACE_FILE SPACE_TIME_FILE ERROR_WARN_FILE SOFTWARE_FILE GPU_LOG_FILE

# Skip software catalog only if it already holds data rows. setup_logging() writes
# the CSV header via _init_csv_headers, so a plain [[ -s ]] check would mis-skip
# on first run and leave software_catalog.csv with header-only content.
_sw_rows=0
[[ -f "${SOFTWARE_FILE:-}" ]] && _sw_rows=$(wc -l < "${SOFTWARE_FILE}" 2>/dev/null || echo 0)
if (( _sw_rows < 2 )); then
	catalog_all_software
else
	log_info "Software catalog already populated, skipping: $SOFTWARE_FILE"
fi
unset _sw_rows

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
# Under orchestration (Nextflow/Snakemake), force serial dispatch unless explicitly
# overridden — the orchestrator manages parallelism at the workflow level.
if [[ -n "${WF_MANAGED_ENV:-}" && -z "${WF_ALLOW_INTERNAL_PARALLEL:-}" ]]; then
    PARALLEL_CONFIGS=1
fi

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
    # Detect wait -n support once (bash 4.3+); avoids conflating
    # "child exited with error" (non-zero rc) with "unsupported flag" (rc=2).
    _has_wait_n=false
    if (( BASH_VERSINFO[0] > 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] >= 3) )); then
        _has_wait_n=true
    fi

    for _cfg_file in "${PIPELINE_CONFIGS[@]}"; do
        # Throttle: wait for a slot when at concurrency limit
        while (( _cfg_active >= PARALLEL_CONFIGS )); do
            if $_has_wait_n; then
                # wait -n returns child's exit code; slot freed regardless of success/failure
                wait -n 2>/dev/null
                (( _cfg_active-- ))
            else
                # bash < 4.3: wait for all children and reset counter
                wait; _cfg_active=0
            fi
        done

        _cfg_tag="${_cfg_file##*/}"
        _cfg_tag="${_cfg_tag%.toml}"
        _cfg_log="${LOG_DIR}/parallel_config_${_cfg_tag}_${RUN_ID}.log"
        _cfg_logs+=("$_cfg_log")

        log_info "  Dispatching: ${_cfg_file##*/}"
        # Export orchestrator flag and all key paths so child process inherits them
        export WF_MANAGED_ENV BASE_DIR ANALYSIS_MODULES_DIR GENE_GROUPS_DIR SRR_CSV_DIR UTILITIES_DIR
        export THREADS ENABLE_GPU AVAILABLE_RAM_GB GPU_VRAM_GB FIGURE_DPI
        __PP_SINGLE_CONFIG="$_cfg_file" bash "$_SELF_SCRIPT" > "$_cfg_log" 2>&1 &
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

# Clear persistent R caches if requested (once before config loop — caches are shared)
if [[ "$CLEAR_CACHE" == "TRUE" ]]; then
    log_info "Clearing persistent pipeline caches..."
    _cache_count=0
    # Sample labels cache
    [[ -f "$SRR_CSV_DIR/.sample_labels_cache.rds" ]] && rm -f "$SRR_CSV_DIR/.sample_labels_cache.rds" && _cache_count=$((_cache_count + 1))
    # Gene name mapping caches (*.namemap.rds beside gene group CSVs)
    while IFS= read -r -d '' _f; do
        rm -f "$_f" && _cache_count=$((_cache_count + 1))
    done < <(find "$GENE_GROUPS_DIR" -name '*.namemap.rds' -print0 2>/dev/null)
    # GPU detection cache (R tempdir varies per session; search common temp roots)
    for _tmp_root in "${TMPDIR:-/tmp}" "${TEMP:-}" "${TMP:-}"; do
        [[ -z "$_tmp_root" || ! -d "$_tmp_root" ]] && continue
        while IFS= read -r -d '' _f; do
            rm -f "$_f" && _cache_count=$((_cache_count + 1))
        done < <(find "$_tmp_root" -maxdepth 2 -name '.gpu_detect_cache.rds' -print0 2>/dev/null)
    done
    log_info "  Cleared $_cache_count cache file(s)"
    unset _cache_count _f _tmp_root
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
    # Preprocessing analyses — no figure output folders (suppress spurious warnings)
    [Stringtie_Matrix]=""
    [Tximport_STAR]=""
    [Tximport_Salmon]=""
    [Tximport_RSEM]=""
)

for CONFIG_FILE in "${PIPELINE_CONFIGS[@]}"; do

    [[ "$CONFIG_FILE" != /* ]] && CONFIG_FILE="$BASE_DIR/$CONFIG_FILE"
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_warn "Config file not found, skipping: $CONFIG_FILE"
        continue
    fi

    # Load config (sets METHODS, ANALYSES, GENE_GROUPS, SRR_DATASETS, etc.)
    # TOML keys are parsed as uppercase bash variables by load_toml
    load_toml "$CONFIG_FILE"
    CURRENT_GENE_GROUP="${GENE_GROUPS[0]:-$CURRENT_GENE_GROUP}"
    export CURRENT_GENE_GROUP
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
    # Respect externally-set OVERWRITE_EXISTING; default to CLEAR_OUTPUT_FOLDER value
    OVERWRITE_EXISTING="${OVERWRITE_EXISTING:-$( [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]] && echo TRUE || echo FALSE )}"
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

    # ── Per-gene-group dispatch ───────────────────────────────────────────────
    # Each gene group becomes its own immediate child of 3_POST_PROC/ (its own
    # CURRENT_GENE_GROUP top-level folder). Processing ONE gene group per pass
    # keeps the deep heatmap folder collision-free as just {dataset} (the
    # "{gene_group}_in_" prefix is dropped in 2_processing_engine.R). Snapshot the
    # full list first: setup_method_env()->_rebuild_exported_arrays() rewrites
    # GENE_GROUPS from GENE_GROUPS_STR (single) inside the loop body, but a bash
    # for-loop already captured this expansion at entry — the copy just makes the
    # intent explicit and robust.
    _ALL_GENE_GROUPS=("${GENE_GROUPS[@]}")
    for CURRENT_GENE_GROUP in "${_ALL_GENE_GROUPS[@]}"; do
    export CURRENT_GENE_GROUP
    export GENE_GROUPS_STR="$CURRENT_GENE_GROUP"   # single group → top-level folder
    log_step "Gene group: $CURRENT_GENE_GROUP (top-level folder)"

    # Clear output folders if requested (rm -rf + mkdir is faster than find -delete on deep trees)
    if [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]]; then
        log_info "Clearing output folders for $MASTER_REFERENCE..."
        # Uses _FOLDER_NAME_MAP declared once before the config loop (avoids per-config re-declaration)
        _clear_targets=()
        for method in "${METHODS[@]}"; do
            output_base="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP/$method/Figure_Outputs"
            [[ -d "$output_base" ]] || continue
            for analysis in "${ANALYSES[@]}"; do
                if [[ -z "${_FOLDER_NAME_MAP[$analysis]+x}" ]]; then
                    log_warn "Unrecognized analysis '$analysis' — not in folder name map, skipping clear"
                    continue
                fi
                folder_name="${_FOLDER_NAME_MAP[$analysis]}"
                # Preprocessing analyses have empty folder name — no figures to clear
                [[ -z "$folder_name" ]] && continue
                target="$output_base/$folder_name/$MASTER_REFERENCE"
                [[ -d "$target" ]] && _clear_targets+=("$target")
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
    # GENE_GROUPS_STR is exported per gene group at the top of the dispatch loop above.
    export ANALYSES_STR="${ANALYSES[*]}" SRR_DATASETS_STR="${SRR_DATASETS[*]}"

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

    # Export utils once (function defs don't change between configs)
    if $_HAS_PARALLEL && [[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && [[ "${_UTILS_EXPORTED:-}" != "true" ]]; then
        export_utils_for_parallel
        _UTILS_EXPORTED=true
    fi
    if $_HAS_PARALLEL && [[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]]; then
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

    # Pre-compute Phase 1 & 2 parallelism limits once per config (not per dataset).
    # METHODS, AVAILABLE_RAM_GB, and THREADS are config-level constants.
    if [[ ${#METHODS[@]} -gt 1 && "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && $_HAS_PARALLEL; then
        # Phase 1 RAM cap: preprocessing R scripts use ~1-2GB each
        _preproc_mem_per_method=1536  # ~1.5GB per preprocessing R process
        _p1_max_by_ram=$(( AVAILABLE_RAM_GB * 1024 * 75 / 100 / _preproc_mem_per_method ))
        (( _p1_max_by_ram < 1 )) && _p1_max_by_ram=1
        _p1_n_methods=${#METHODS[@]}
        (( _p1_n_methods > _p1_max_by_ram )) && _p1_n_methods=$_p1_max_by_ram
        export _p1_n_methods

        # Phase 2 thread/RAM calculation
        if [[ ${#HEAVY_ANALYSES[@]} -gt 0 ]]; then
            _n_methods=${#METHODS[@]}
            _heavy_mem_per_method=3072  # ~3GB per heavy R analysis
            _max_by_ram=$(( AVAILABLE_RAM_GB * 1024 * 75 / 100 / _heavy_mem_per_method ))
            (( _max_by_ram < 1 )) && _max_by_ram=1
            (( _n_methods > _max_by_ram )) && _n_methods=$_max_by_ram
            _threads_per_method=$(( THREADS / _n_methods ))
            (( _threads_per_method < 1 )) && _threads_per_method=1
            # Worker function: run all heavy analyses for one method with reduced thread count
            _run_heavy_for_method() {
                local _method="$1" _ref="$2"
                shift 2
                export THREADS="$_threads_per_method"
                run_batched_analyses "$_method" "$_ref" "$@"
            }
            export -f _run_heavy_for_method
            export _threads_per_method _n_methods
        fi
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
            log_info "Phase 1: Preprocessing (parallel across ${_p1_n_methods} methods, RAM-capped)"
            parallel \
                -j "$_p1_n_methods" \
                --halt soon,fail,30% \
                --joblog "$LOG_DIR/parallel_preproc_${dataset}_${BASHPID:-$$}.log" \
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
        # _n_methods, _threads_per_method, and _run_heavy_for_method are hoisted
        # to config-level above the dataset loop to avoid per-dataset re-declaration.
        if [[ ${#HEAVY_ANALYSES[@]} -gt 0 ]]; then
            if [[ ${#METHODS[@]} -gt 1 && "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && $_HAS_PARALLEL; then
                log_info "Phase 2: Heavy analyses (parallel across ${_n_methods} methods, ${_threads_per_method} threads each): ${HEAVY_ANALYSES[*]}"
                parallel \
                    -j "$_n_methods" \
                    --halt soon,fail,30% \
                    --joblog "$LOG_DIR/parallel_heavy_${dataset}_${BASHPID:-$$}.log" \
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
                    --joblog "$LOG_DIR/parallel_figures_${dataset}_${BASHPID:-$$}.log" \
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

    done  # ── end per-gene-group dispatch ──

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
    _viewer_script="$BASE_DIR/modules_gea/c_post_processing/utilities/generate_html_viewer.py"
    _post_proc_dir="$BASE_DIR/II_RESULTS/3_POST_PROC/$CURRENT_GENE_GROUP"
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
