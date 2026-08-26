#!/bin/bash
#===============================================================================
# CROSS-METHOD CONCORDANCE ANALYSIS
#===============================================================================
# Compares gene expression quantification across all 5 alignment/quantification
# methods (M1-M5) for a given reference genome. Produces correlation matrices,
# discordant gene lists, and a unified report.
#
# Steps:
#   1. Load & harmonize TPM matrices from all methods
#   2. Compute pairwise Spearman correlations
#   3. Generate unified Markdown report
#
# Usage:
#   bash run_concordance_analysis.sh [config_file|config_cross|config_cross_dir]
#
#   If no config_file is provided, uses internal defaults for GPE001970.
#
#   Config crosses can be stored in: config/4_concordance_combination/
#   and selected via:
#     - positional arg: cross filename or cross basename
#     - env var: CONCORDANCE_CONFIG_CROSSES="cross1,cross2,..."
#
# Prerequisites:
#   - Alignment results for M1-M5 in II_RESULTS/2_ALIGNMENT_RESULTs/
#   - Post-processing matrices in II_RESULTS/3_POST_PROC/ (for M4/M5 fallbacks)
#   - R packages: ComplexHeatmap, circlize, grid
#   - conda env "gea" with all dependencies
#===============================================================================

set -euo pipefail

# Script-level child PID array — must be declared before traps so cleanup can reap
# child processes spawned in multi-config mode (prevents orphans on SIGTERM/SIGINT)
declare -a _cfg_pids=()

# Trap handler: log errors on unexpected exit (aids debugging in orchestrated contexts)
_concordance_cleanup() {
    local rc=$?
    # Kill any background child config processes to prevent orphans on HPC
    if [[ ${#_cfg_pids[@]} -gt 0 ]]; then
        kill "${_cfg_pids[@]}" 2>/dev/null || true
        wait "${_cfg_pids[@]}" 2>/dev/null || true
    fi
    if [[ $rc -ne 0 ]]; then
        echo "[ERROR] run_concordance_analysis.sh exited with code $rc" >&2
        # Use log_error if available (may not be sourced yet at early failure)
        type -t log_error &>/dev/null && log_error "Concordance analysis failed (exit $rc)"
    fi
    # Clean up background tee/sed processes from logging redirections (prevents zombie
    # processes in Nextflow/Snakemake containers that would stall work-dir cleanup)
    type -t _logging_cleanup_bg &>/dev/null && _logging_cleanup_bg
    exit $rc  # propagate original exit code so orchestrators (Nextflow/Snakemake) see failures
}
trap _concordance_cleanup EXIT
trap 'type -t _logging_cleanup_bg &>/dev/null && _logging_cleanup_bg; exit 143' TERM
trap 'type -t _logging_cleanup_bg &>/dev/null && _logging_cleanup_bg; exit 130' INT

# Resolve script directory: honour pre-set BASE_DIR from orchestrators (Nextflow/Snakemake)
if [[ -z "${BASE_DIR:-}" ]]; then
    SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
    [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
    SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] run_concordance_analysis.sh: Failed to resolve script directory" >&2; exit 1; }
    BASE_DIR="${SCRIPT_DIR}"
else
    # Validate BASE_DIR is absolute; resolve if relative (defensive against misconfigured orchestrators)
    if [[ "$BASE_DIR" != /* ]]; then
        BASE_DIR="$(cd "$BASE_DIR" 2>/dev/null && pwd)" || { echo "[ERROR] BASE_DIR is set but invalid: $BASE_DIR" >&2; exit 1; }
    fi
    SCRIPT_DIR="$BASE_DIR"
fi
_SELF_SCRIPT="$SCRIPT_DIR/${BASH_SOURCE[0]##*/}"
REPORT_BASE="${REPORT_BASE:-${BASE_DIR}/II_RESULTS/4_CONCORDANCE_ANALYSIS}"

#===============================================================================
# CONDA ENVIRONMENT
#===============================================================================

# Activate conda env only if not already active — child processes (bash "$0")
# inherit CONDA_DEFAULT_ENV but not the conda PATH entries from the parent shell,
# so we must always run the hook. Skip only when PATH is already configured.
# O(1) string match avoids ~0.3-0.5s conda hook overhead per child process.
# Skip conda activation when orchestrator manages the environment (Nextflow/Snakemake)
if [[ -z "${WF_MANAGED_ENV:-}" ]]; then
    if [[ "${CONDA_DEFAULT_ENV:-}" != "gea" ]] || ! command -v Rscript &>/dev/null; then
        # In non-interactive shells (e.g. WSL2 child processes), conda is not in PATH
        # because ~/.bashrc is not sourced. Bootstrap it from known install locations.
        if ! command -v conda &>/dev/null; then
            for _conda_prefix in \
                "${HOME}/miniconda3" \
                "${HOME}/anaconda3" \
                "/opt/conda" \
                "/opt/miniconda3" \
                "/opt/anaconda3"; do
                if [[ -f "${_conda_prefix}/etc/profile.d/conda.sh" ]]; then
                    # shellcheck source=/dev/null
                    source "${_conda_prefix}/etc/profile.d/conda.sh"
                    break
                fi
            done
        fi
        eval "$(conda shell.bash hook 2>/dev/null)" 2>/dev/null || true
        conda activate gea 2>/dev/null || true
    fi
fi

# Source logging utilities for consistent output
source "${SCRIPT_DIR}/modules_gea/logging/logging_utils.sh" 2>/dev/null || {
    # Minimal fallback if logging module unavailable
    log_info()  { echo "[INFO]  $*"; }
    log_warn()  { echo "[WARN]  $*"; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_step()  { echo ""; echo "==> $*"; }
}

# Source TOML parser
source "${SCRIPT_DIR}/config/shared/toml_parser.sh" || {
    log_error "Failed to source TOML parser at ${SCRIPT_DIR}/config/shared/toml_parser.sh"
    exit 1
}

#===============================================================================
# CONFIGURATION (override via env, config file, or associative array below)
#===============================================================================

CONFIG_INPUT="${1:-}"
# Absolutize relative config path against BASE_DIR so the pipeline works when
# CWD differs from project root (e.g., Nextflow scratch workDir, Snakemake shadow).
[[ -n "$CONFIG_INPUT" && "$CONFIG_INPUT" != /* ]] && CONFIG_INPUT="${BASE_DIR}/${CONFIG_INPUT}"
CONFIG_CROSS_DIR="${CONCORDANCE_CONFIG_CROSS_DIR:-${BASE_DIR}/config/4_concordance_combination}"
if [[ ! -d "$CONFIG_CROSS_DIR" ]]; then
    log_error "Concordance config directory not found: $CONFIG_CROSS_DIR"
    log_error "Ensure BASE_DIR is correct or set CONCORDANCE_CONFIG_CROSS_DIR explicitly."
    exit 1
fi

# Analysis steps to run (comment out entries to skip)
# ─────────────────────────────────────────────────────
# Guard: only set defaults if ANALYSES is not already defined as an array
# O(1) builtin attribute check — avoids $(declare -p) subshell fork
if ! [[ -v ANALYSES && "${ANALYSES@a}" == *a* ]]; then
    ANALYSES=(
        "Load_Matrices"                 # Step 1: Load & harmonize matrices (required by all others)
        "Quantification_Concordance"    # Step 2: Compare quantification across methods
        "Generate_Report"               # Step 3: Produce concordance report
    )
fi

# Figure resolution in DPI (300–600)
FIGURE_DPI="${FIGURE_DPI:-300}"

# Clear previous outputs before running (respect parent env in child dispatch)
CLEAR_LOGS="${CLEAR_LOGS:-TRUE}"
CLEAR_OUTPUT_FOLDER="${CLEAR_OUTPUT_FOLDER:-TRUE}"
CLEAR_CACHE="${CLEAR_CACHE:-FALSE}"

# Optional curated config list (comment in/out as needed).
# Entries can be:
#   - absolute/relative file paths
#   - cross basenames from config/4_concordance_combination (with or without .sh)
# Load order matters: later entries override earlier ones.
CONCORDANCE_CONFIGS=(
    #"defaults.toml"
    "cross_genomes_vs_genomes.toml"
    "cross_equivalent_gene_between_genomes.toml"
    "cross_methods_vs_methods.toml"
    "cross_within_gene_group_vs_within_gene_group.toml"
    #"cross_full_factorial_example.toml"
)

REINVOKE_ARGS=()
if [[ -n "$CONFIG_INPUT" ]]; then
    REINVOKE_ARGS+=("$CONFIG_INPUT")
fi

source_config_file() {
    local cfg_path="$1"
    if [[ -f "$cfg_path" ]]; then
        log_info "Loading config: $cfg_path"
        if [[ "$cfg_path" == *.toml ]]; then
            load_toml "$cfg_path"
        else
            # Legacy .sh fallback
            # shellcheck disable=SC1090
            source "$cfg_path"
        fi
        return 0
    fi
    return 1
}

load_config_entry() {
    local entry="$1"

    # Trim whitespace-only entries.
    entry="${entry#${entry%%[![:space:]]*}}"
    entry="${entry%${entry##*[![:space:]]}}"
    [[ -z "$entry" ]] && return 0

    # Loop over candidate paths — try .toml first, then .sh fallback
    local path
    for path in "$entry" "${CONFIG_CROSS_DIR}/${entry}" "${CONFIG_CROSS_DIR}/${entry}.toml" "${CONFIG_CROSS_DIR}/${entry}.sh"; do
        [[ -f "$path" ]] && { source_config_file "$path" || { log_error "Failed to parse config: $path"; return 1; }; return 0; }
    done
    log_warn "Config entry not found: ${entry}"
}

# Load manually curated config list.
# In single-config child mode, load only the specified config; otherwise load all.
if [[ -n "${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" ]]; then
    load_config_entry "${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG}"
elif [[ -v CONCORDANCE_CONFIGS && "${CONCORDANCE_CONFIGS@a}" == *a* && ${#CONCORDANCE_CONFIGS[@]} -gt 0 ]]; then
    for _cfg_entry in "${CONCORDANCE_CONFIGS[@]}"; do
        load_config_entry "$_cfg_entry"
    done
    unset _cfg_entry
fi

# Load optional config file/cross/dir (first positional argument)
if [[ -n "$CONFIG_INPUT" ]]; then
    if [[ -f "$CONFIG_INPUT" ]]; then
        source_config_file "$CONFIG_INPUT" || { log_error "Failed to parse config: $CONFIG_INPUT"; exit 1; }
    elif [[ -d "$CONFIG_INPUT" ]]; then
        # Bash glob + mapfile avoids find+sort subprocess pair. O(F) where F = config files.
        _cfg_files=()
        for _cfg in "$CONFIG_INPUT"/*.toml "$CONFIG_INPUT"/*.sh; do
            [[ -f "$_cfg" ]] && _cfg_files+=("$_cfg")
        done
        # Sort for deterministic order (globs are locale-sorted but toml/sh interleave)
        mapfile -t _cfg_files < <(printf '%s\n' "${_cfg_files[@]}" | sort)
        for _cfg in "${_cfg_files[@]}"; do
            source_config_file "$_cfg" || { log_error "Failed to parse config: $_cfg"; exit 1; }
        done
        unset _cfg_files
    elif [[ -f "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}" ]]; then
        source_config_file "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}" || { log_error "Failed to parse config: ${CONFIG_CROSS_DIR}/${CONFIG_INPUT}"; exit 1; }
    elif [[ -f "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.toml" ]]; then
        source_config_file "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.toml" || { log_error "Failed to parse config: ${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.toml"; exit 1; }
    elif [[ -f "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.sh" ]]; then
        source_config_file "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.sh" || { log_error "Failed to parse config: ${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.sh"; exit 1; }
    else
        log_warn "Config input not found: ${CONFIG_INPUT} (continuing with defaults)"
    fi
fi

# Optionally source additional config crosses from config/4_concordance_combination
# Example: CONCORDANCE_CONFIG_CROSSES="defaults,cross_methods_vs_methods"
if [[ -n "${CONCORDANCE_CONFIG_CROSSES:-}" ]]; then
    IFS=',' read -r -a _cfg_crosses <<< "$CONCORDANCE_CONFIG_CROSSES"
    for _cfg_cross in "${_cfg_crosses[@]}"; do
        load_config_entry "$_cfg_cross"
    done
    unset _cfg_crosses _cfg_cross
fi

# Internal child-run overrides (applied after config sourcing)
[[ -n "${__CONCORDANCE_OVERRIDE_MASTER_REFERENCE:-}" ]] && MASTER_REFERENCE="${__CONCORDANCE_OVERRIDE_MASTER_REFERENCE}"
[[ -n "${__CONCORDANCE_OVERRIDE_METHODS:-}" ]] && METHODS="${__CONCORDANCE_OVERRIDE_METHODS}"
[[ -n "${__CONCORDANCE_OVERRIDE_GENE_GROUPS:-}" ]] && GENE_GROUPS="${__CONCORDANCE_OVERRIDE_GENE_GROUPS}"
[[ -n "${__CONCORDANCE_OVERRIDE_REPORT_BASE:-}" ]] && REPORT_BASE="${__CONCORDANCE_OVERRIDE_REPORT_BASE}"
[[ -n "${__CONCORDANCE_OVERRIDE_RUN_ALL_MASTER_REFERENCES:-}" ]] && RUN_ALL_MASTER_REFERENCES="${__CONCORDANCE_OVERRIDE_RUN_ALL_MASTER_REFERENCES}"
[[ -n "${__CONCORDANCE_OVERRIDE_RUN_ALL_METHOD_COMBINATIONS:-}" ]] && RUN_ALL_METHOD_COMBINATIONS="${__CONCORDANCE_OVERRIDE_RUN_ALL_METHOD_COMBINATIONS}"
[[ -n "${__CONCORDANCE_OVERRIDE_RUN_ALL_GENE_GROUP_COMBINATIONS:-}" ]] && RUN_ALL_GENE_GROUP_COMBINATIONS="${__CONCORDANCE_OVERRIDE_RUN_ALL_GENE_GROUP_COMBINATIONS}"

sanitize_tag() {
    # Pure bash via nameref: replace separators with _, strip non-alnum.
    # Nameref avoids $(sanitize_tag ...) subshell fork — caller passes output variable name.
    # Usage: sanitize_tag "input_string" result_var
    local -n _st_out=$2
    _st_out="${1//[,|\/[:space:]]/_}"
    _st_out="${_st_out//[^[:alnum:]_.-]/}"
    [[ -z "$_st_out" ]] && _st_out="combo"
}

# Join array elements with a custom separator via nameref — avoids $(IFS=...) subshell fork.
# Usage: _join_with result_var "separator" "${array[@]}"
_join_with() {
    local -n _jw_out=$1; local IFS="$2"; shift 2; _jw_out="$*"
}

#===============================================================================
# PATHS (must be set before dispatch blocks so parent logging works)
#===============================================================================

REPORT_BASE="${REPORT_BASE:-${BASE_DIR}/II_RESULTS/4_CONCORDANCE_ANALYSIS}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPORT_BASE}}"
ALIGNMENT_BASE="${ALIGNMENT_BASE:-${BASE_DIR}/II_RESULTS/2_ALIGNMENT_RESULTs}"
POST_PROC_BASE="${POST_PROC_BASE:-${BASE_DIR}/II_RESULTS/3_POST_PROC/${CURRENT_GENE_GROUP:-_active}}"
ANALYSIS_MODULES_DIR="${BASE_DIR}/modules_gea/c_post_processing/analysis_modules"
UTILITIES_DIR="${BASE_DIR}/modules_gea/c_post_processing/utilities"
CONCORDANCE_SCRIPT_DIR="${BASE_DIR}/modules_gea/c_post_processing/cross_method_concordance"

# Clear previous outputs if requested
if [[ "${CLEAR_OUTPUT_FOLDER:-FALSE}" == "TRUE" && -d "${OUTPUT_DIR}" ]]; then
    log_info "Clearing previous output folder: ${OUTPUT_DIR}"
    rm -rf "${OUTPUT_DIR}"
fi
# Logs always live under the top-level REPORT_BASE (II_RESULTS/4_CONCORDANCE_ANALYSIS/logs/),
# even when child processes override REPORT_BASE to per-config subdirectories.
# NOTE: When OUTPUT_DIR == REPORT_BASE (default), the rm -rf above already deleted logs.
# This block only has effect when OUTPUT_DIR is a subdirectory of REPORT_BASE.
CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE:-${REPORT_BASE}/logs}"
if [[ "${CLEAR_LOGS:-FALSE}" == "TRUE" && -d "${CONCORDANCE_LOG_BASE}" ]]; then
    log_info "Clearing previous logs: ${CONCORDANCE_LOG_BASE}"
    rm -rf "${CONCORDANCE_LOG_BASE}"
fi

# Clear persistent R caches if requested
if [[ "${CLEAR_CACHE:-FALSE}" == "TRUE" ]]; then
    log_info "Clearing persistent pipeline caches..."
    _cache_count=0
    _gg_dir="${BASE_DIR}/I_INPUTS/inputs/3_post_proc_inputs/gene_groups_csv"
    # Gene name mapping caches (*.namemap.rds beside gene group CSVs)
    while IFS= read -r -d '' _f; do
        rm -f "$_f" && _cache_count=$((_cache_count + 1))
    done < <(find "$_gg_dir" -name '*.namemap.rds' -print0 2>/dev/null)
    # GPU detection cache (R tempdir varies per session; search common temp roots)
    for _tmp_root in "${TMPDIR:-/tmp}" "${TEMP:-}" "${TMP:-}"; do
        [[ -z "$_tmp_root" || ! -d "$_tmp_root" ]] && continue
        while IFS= read -r -d '' _f; do
            rm -f "$_f" && _cache_count=$((_cache_count + 1))
        done < <(find "$_tmp_root" -maxdepth 2 -name '.gpu_detect_cache.rds' -print0 2>/dev/null)
    done
    log_info "  Cleared $_cache_count cache file(s)"
    unset _cache_count _f _gg_dir _tmp_root
fi

mkdir -p "${OUTPUT_DIR}" || {
    log_error "Failed to create output directory: ${OUTPUT_DIR}"
    exit 1
}

# Set up structured logging (mirrors run_post_processing.sh)
LOG_DIR="${CONCORDANCE_LOG_BASE}/log_files"
TIME_DIR="${CONCORDANCE_LOG_BASE}/time_logs"
SPACE_DIR="${CONCORDANCE_LOG_BASE}/space_logs"
SPACE_TIME_DIR="${CONCORDANCE_LOG_BASE}/space_time_logs"
ERROR_WARN_DIR="${CONCORDANCE_LOG_BASE}/error_warn_logs"
SOFTWARE_CATALOG_DIR="${CONCORDANCE_LOG_BASE}/software_catalogs"
GPU_LOG_DIR="${CONCORDANCE_LOG_BASE}/gpu_log"
export LOG_DIR TIME_DIR SPACE_DIR SPACE_TIME_DIR ERROR_WARN_DIR SOFTWARE_CATALOG_DIR GPU_LOG_DIR

# Ensure log directories exist before setup_logging (which may skip mkdir
# if LOGGING_INITIALIZED was inherited from a parent process).
mkdir -p "$LOG_DIR" "$TIME_DIR" "$SPACE_DIR" "$SPACE_TIME_DIR" \
         "$ERROR_WARN_DIR" "$SOFTWARE_CATALOG_DIR" "$GPU_LOG_DIR" 2>/dev/null || true

# Reset so setup_logging re-derives file paths for this concordance run
unset LOGGING_INITIALIZED

if declare -f setup_logging &>/dev/null; then
    setup_logging "$CLEAR_LOGS"
    export LOG_FILE TIME_FILE SPACE_FILE SPACE_TIME_FILE ERROR_WARN_FILE SOFTWARE_FILE GPU_LOG_FILE
fi

# Maximum concurrent background concordance jobs (prevents CPU/memory exhaustion on
# constrained systems when factorial modes launch many combinations).
# O(min(N, MAX_CONCORDANCE_JOBS)) wall-clock vs O(N) unthrottled.
MAX_CONCORDANCE_JOBS="${MAX_CONCORDANCE_JOBS:-4}"

# Detect wait -n support once (bash 4.3+); avoids conflating
# "child exited with error" (non-zero rc) with "unsupported flag" (rc=2).
_has_wait_n=false
if (( BASH_VERSINFO[0] > 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] >= 3) )); then
    _has_wait_n=true
fi

# Throttle helper: waits until the number of tracked PIDs drops below MAX_CONCORDANCE_JOBS.
# Usage: _throttle_pids <array_name>
# Complexity: O(J) per call where J = active jobs; total O(N) amortized across N launches.
_throttle_pids() {
    local -n _pids_ref=$1
    while [[ ${#_pids_ref[@]} -ge $MAX_CONCORDANCE_JOBS ]]; do
        # Use `wait -n` (Bash 4.3+) to block until any child exits — avoids busy-wait polling.
        # Falls back to poll+sleep for older Bash versions.
        if $_has_wait_n; then
            # wait -n returns child's exit code; slot freed regardless of success/failure
            wait -n "${_pids_ref[@]}" 2>/dev/null || true
        else
            # bash < 4.3: poll+sleep to avoid busy-wait
            sleep 0.5
        fi
        # Compact the PID array: keep only still-running PIDs
        local _still_running=()
        for _p in "${_pids_ref[@]}"; do
            if kill -0 "$_p" 2>/dev/null; then
                _still_running+=("$_p")
            else
                wait "$_p" 2>/dev/null || true
            fi
        done
        _pids_ref=("${_still_running[@]}")
    done
}

# Optional multi-config mode:
# When multiple CONCORDANCE_CONFIGS are active and we're not already in single-config
# child mode, dispatch a separate child process per config with an isolated output folder.
if [[ -z "${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" ]] && \
   [[ -v CONCORDANCE_CONFIGS && "${CONCORDANCE_CONFIGS@a}" == *a* ]] && \
   [[ ${#CONCORDANCE_CONFIGS[@]} -gt 1 ]]; then
    _parent_report_base="${REPORT_BASE}"
    _overall_rc=0

    log_step "MULTI-CONFIG CONCORDANCE MODE"
    log_info "Found ${#CONCORDANCE_CONFIGS[@]} config crosses"

    _cfg_pids=()
    for _cfg in "${CONCORDANCE_CONFIGS[@]}"; do
        _throttle_pids _cfg_pids
        _cfg_basename="${_cfg##*/}"
        _cfg_tag="${_cfg_basename%.*}_concordance"
        log_step "Launching concordance for config: ${_cfg_basename}"
        # Give each child a unique RUN_ID to prevent log file collisions
        # when multiple children start within the same second.
        printf -v _child_run_id '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || _child_run_id=$(date +%Y%m%d_%H%M%S)
        _child_run_id="${_child_run_id}_${_cfg_tag}"
        __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/${_cfg_tag}" \
        __CONCORDANCE_OVERRIDE_SINGLE_CONFIG="${_cfg}" \
        CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE}" \
        RUN_ID="${_child_run_id}" LOGGING_INITIALIZED="" \
        CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE CLEAR_CACHE=FALSE \
        bash "$_SELF_SCRIPT" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
        _cfg_pids+=($!)
    done
    for _pid in "${_cfg_pids[@]}"; do
        wait "$_pid" || _overall_rc=1
    done

    if [[ $_overall_rc -ne 0 ]]; then
        log_error "One or more config crosses failed in multi-config mode"
        exit 1
    fi
    log_info "All config crosses completed successfully"
    exit 0
fi

# Optional multi-reference mode:
# - RUN_ALL_MASTER_REFERENCES=TRUE: iterate over MASTER_REFERENCES and run once per reference
# - Default FALSE: run only one reference (MASTER_REFERENCE or first MASTER_REFERENCES entry)
RUN_ALL_MASTER_REFERENCES="${RUN_ALL_MASTER_REFERENCES:-FALSE}"
if [[ "${RUN_ALL_MASTER_REFERENCES^^}" == "TRUE" ]]; then
    if [[ -v MASTER_REFERENCES && "${MASTER_REFERENCES@a}" == *a* && ${#MASTER_REFERENCES[@]} -gt 0 ]]; then
        _parent_report_base="${REPORT_BASE}"
        _overall_rc=0

        log_step "MULTI-REFERENCE CONCORDANCE MODE"
        log_info "Found ${#MASTER_REFERENCES[@]} references in MASTER_REFERENCES"

        # Run references in parallel — each writes to isolated output dir
        # Reduces wall-clock from O(R × time) to O(max(1, R/MAX_CONCORDANCE_JOBS) × time)
        _ref_pids=()
        for _ref in "${MASTER_REFERENCES[@]}"; do
            _throttle_pids _ref_pids
            log_step "Launching concordance for reference: ${_ref}"
            printf -v _child_run_id '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || _child_run_id=$(date +%Y%m%d_%H%M%S)
            _child_run_id="${_child_run_id}_${_ref}"
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/${_ref}" \
            __CONCORDANCE_OVERRIDE_MASTER_REFERENCE="${_ref}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_MASTER_REFERENCES="FALSE" \
            __CONCORDANCE_OVERRIDE_SINGLE_CONFIG="${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" \
            CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE}" \
            RUN_ID="${_child_run_id}" LOGGING_INITIALIZED="" \
            CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE CLEAR_CACHE=FALSE \
            bash "$_SELF_SCRIPT" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
            _ref_pids+=($!)
        done
        for _pid in "${_ref_pids[@]}"; do
            wait "$_pid" || _overall_rc=1
        done

        if [[ $_overall_rc -ne 0 ]]; then
            log_error "One or more references failed in multi-reference mode"
            exit 1
        fi
        log_info "All references completed successfully"
        exit 0
    else
        log_warn "RUN_ALL_MASTER_REFERENCES=TRUE but MASTER_REFERENCES array is empty/unset; continuing single-reference run"
    fi
fi

# Optional method-combination mode:
# - RUN_ALL_METHOD_COMBINATIONS=TRUE: iterate over METHOD_COMBINATIONS entries
# - METHOD_COMBINATIONS entry format: comma or pipe separated method IDs
#   e.g., "M1_HISAT2_RefGuided,M3_STAR_Align"
RUN_ALL_METHOD_COMBINATIONS="${RUN_ALL_METHOD_COMBINATIONS:-FALSE}"
if [[ "${RUN_ALL_METHOD_COMBINATIONS^^}" == "TRUE" ]]; then
    if [[ -v METHOD_COMBINATIONS && "${METHOD_COMBINATIONS@a}" == *a* && ${#METHOD_COMBINATIONS[@]} -gt 0 ]]; then
        _parent_report_base="${REPORT_BASE}"
        _overall_rc=0

        log_step "METHOD-COMBINATION CONCORDANCE MODE"
        log_info "Found ${#METHOD_COMBINATIONS[@]} method combinations"

        # Run method combinations in parallel — each writes to isolated output dir
        # Reduces wall-clock from O(C × time) to O(max(1, C/MAX_CONCORDANCE_JOBS) × time)
        # Throttled to MAX_CONCORDANCE_JOBS concurrent processes to limit CPU/memory pressure.
        _combo_pids=()
        for _combo in "${METHOD_COMBINATIONS[@]}"; do
            _throttle_pids _combo_pids
            # Single-pass: replace both comma and pipe separators with space
            _methods="${_combo//[,|]/ }"
            sanitize_tag "$_combo" _combo_tag
            log_step "Launching concordance for methods: ${_methods}"
            printf -v _child_run_id '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || _child_run_id=$(date +%Y%m%d_%H%M%S)
            _child_run_id="${_child_run_id}_methods_${_combo_tag}"
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/methods_${_combo_tag}" \
            __CONCORDANCE_OVERRIDE_METHODS="${_methods}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_METHOD_COMBINATIONS="FALSE" \
            __CONCORDANCE_OVERRIDE_SINGLE_CONFIG="${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" \
            CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE}" \
            RUN_ID="${_child_run_id}" LOGGING_INITIALIZED="" \
            CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE CLEAR_CACHE=FALSE \
            bash "$_SELF_SCRIPT" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
            _combo_pids+=($!)
        done
        for _pid in "${_combo_pids[@]}"; do
            wait "$_pid" || _overall_rc=1
        done

        if [[ $_overall_rc -ne 0 ]]; then
            log_error "One or more method combinations failed"
            exit 1
        fi
        log_info "All method combinations completed successfully"
        exit 0
    else
        log_warn "RUN_ALL_METHOD_COMBINATIONS=TRUE but METHOD_COMBINATIONS array is empty/unset; continuing standard run"
    fi
fi

# Optional gene-group-combination mode:
# - RUN_ALL_GENE_GROUP_COMBINATIONS=TRUE: iterate over GENE_GROUP_COMBINATIONS entries
# - GENE_GROUP_COMBINATIONS entry format: comma or pipe separated gene-group basenames
#   e.g., "SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF"
RUN_ALL_GENE_GROUP_COMBINATIONS="${RUN_ALL_GENE_GROUP_COMBINATIONS:-FALSE}"
if [[ "${RUN_ALL_GENE_GROUP_COMBINATIONS^^}" == "TRUE" ]]; then
    if [[ -v GENE_GROUP_COMBINATIONS && "${GENE_GROUP_COMBINATIONS@a}" == *a* && ${#GENE_GROUP_COMBINATIONS[@]} -gt 0 ]]; then
        _parent_report_base="${REPORT_BASE}"
        _overall_rc=0

        log_step "GENE-GROUP-COMBINATION CONCORDANCE MODE"
        log_info "Found ${#GENE_GROUP_COMBINATIONS[@]} gene-group combinations"

        # Run gene-group combinations in parallel — each writes to isolated output dir
        # Reduces wall-clock from O(C × time) to O(max(1, C/MAX_CONCORDANCE_JOBS) × time)
        # Throttled to MAX_CONCORDANCE_JOBS concurrent processes to limit CPU/memory pressure.
        _gg_combo_pids=()
        for _combo in "${GENE_GROUP_COMBINATIONS[@]}"; do
            _throttle_pids _gg_combo_pids
            _groups="${_combo//|/,}"
            sanitize_tag "$_combo" _combo_tag
            log_step "Launching concordance for gene groups: ${_groups}"
            printf -v _child_run_id '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || _child_run_id=$(date +%Y%m%d_%H%M%S)
            _child_run_id="${_child_run_id}_genes_${_combo_tag}"
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/genes_${_combo_tag}" \
            __CONCORDANCE_OVERRIDE_GENE_GROUPS="${_groups}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_GENE_GROUP_COMBINATIONS="FALSE" \
            __CONCORDANCE_OVERRIDE_SINGLE_CONFIG="${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" \
            CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE}" \
            RUN_ID="${_child_run_id}" LOGGING_INITIALIZED="" \
            CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE CLEAR_CACHE=FALSE \
            bash "$_SELF_SCRIPT" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
            _gg_combo_pids+=($!)
        done
        for _pid in "${_gg_combo_pids[@]}"; do
            wait "$_pid" || _overall_rc=1
        done

        if [[ $_overall_rc -ne 0 ]]; then
            log_error "One or more gene-group combinations failed"
            exit 1
        fi
        log_info "All gene-group combinations completed successfully"
        exit 0
    else
        log_warn "RUN_ALL_GENE_GROUP_COMBINATIONS=TRUE but GENE_GROUP_COMBINATIONS array is empty/unset; continuing standard run"
    fi
fi

# Reference genome to compare across methods
# If MASTER_REFERENCE was not explicitly provided, and a sourced config set
# MASTER_REFERENCES as an array, take the first element
# (matches run_post_processing.sh behaviour).
if [[ -z "${MASTER_REFERENCE:-}" ]] && [[ -v MASTER_REFERENCES && "${MASTER_REFERENCES@a}" == *a* ]]; then
    MASTER_REFERENCE="${MASTER_REFERENCES[0]}"
fi
# Default must match 0_concordance_config.R and 0_shared_config.R ("Eggplant_V4.1")
MASTER_REFERENCE="${MASTER_REFERENCE:-Eggplant_V4.1}"

# Methods to compare (space-separated string).
# If a sourced config set METHODS as an array, flatten it to a string.
# Must unset array before reassignment — bash arrays cannot be exported to child processes.
if [[ -v METHODS && "${METHODS@a}" == *a* ]]; then
    _methods_str="${METHODS[*]}"
    unset METHODS
    METHODS="$_methods_str"
    unset _methods_str
fi
METHODS="${METHODS:-
    M1_HISAT2_RefGuided 
    M2_HISAT2_DeNovo 
    M3_STAR_Align 
    M4_Salmon_Saf 
    M5_RSEM_Bowtie2
}"

# Enforce method/reference compatibility (cross_method mode only):
# - *_genome references:      M1, M3
# - *transcript* references:  M2, M4, M5
# Set ENFORCE_REFERENCE_METHOD_COMPATIBILITY=FALSE to disable filtering.
ENFORCE_REFERENCE_METHOD_COMPATIBILITY="${ENFORCE_REFERENCE_METHOD_COMPATIBILITY:-TRUE}"
if [[ "${ENFORCE_REFERENCE_METHOD_COMPATIBILITY^^}" == "TRUE" && "${CONCORDANCE_MODE:-cross_method}" == "cross_method" ]]; then
    _methods_filtered=""
    _methods_dropped=""
    _has_rule=0

    if [[ "${MASTER_REFERENCE}" == *_genome ]]; then
        _has_rule=1
        _allowed_methods=("M1_HISAT2_RefGuided" "M3_STAR_Align")
    elif [[ "${MASTER_REFERENCE}" == *transcript* ]]; then
        _has_rule=1
        _allowed_methods=("M2_HISAT2_DeNovo" "M4_Salmon_Saf" "M5_RSEM_Bowtie2")
    fi

    if [[ $_has_rule -eq 1 ]]; then
        # O(M) single-pass with pattern match instead of O(M × A) nested loop
        _allowed_pat=""
        printf -v _allowed_pat '|%s' "${_allowed_methods[@]}"
        _allowed_pat="@(${_allowed_pat:1})"  # extglob pattern: @(M1_...|M3_...)
        shopt -s extglob
        for _m in ${METHODS}; do
            # shellcheck disable=SC2053
            if [[ "$_m" == $_allowed_pat ]]; then
                _methods_filtered="${_methods_filtered:+${_methods_filtered} }${_m}"
            else
                _methods_dropped="${_methods_dropped:+${_methods_dropped} }${_m}"
            fi
        done
        shopt -u extglob

        if [[ -z "$_methods_filtered" ]]; then
            log_error "No compatible methods left for MASTER_REFERENCE='${MASTER_REFERENCE}' after filtering"
            exit 1
        fi

        if [[ -n "$_methods_dropped" ]]; then
            log_warn "Dropped incompatible methods for ${MASTER_REFERENCE}: ${_methods_dropped}"
            log_info "Using compatible methods: ${_methods_filtered}"
        fi

        METHODS="$_methods_filtered"
    fi

    unset _methods_filtered _methods_dropped _has_rule _allowed_methods _allowed_pat _m
fi

# Method-specific reference directory names
# M1/M3 align to genome; M2/M4/M5 align to transcriptome
# Derive defaults from MASTER_REFERENCE instead of hardcoding GPE001970,
# so that sourcing a non-GPE001970 config produces correct paths.
#
# Supports three MASTER_REFERENCE forms:
#   "GPE001970_genome"       → genome=GPE001970_genome,      transcript=GPE001970_transcripts
#   "GPE001970_transcripts"  → genome=GPE001970_genome,      transcript=GPE001970_transcripts
#   "GPE001970"              → genome=GPE001970_genome,      transcript=GPE001970_transcripts
# The bare-accession form makes cross-method concordance reference-agnostic:
# each method resolves to its appropriate reference type automatically.
if [[ "${MASTER_REFERENCE}" == *_genome ]]; then
    _genome_ref="${MASTER_REFERENCE}"
    _transcript_ref="${MASTER_REFERENCE%_genome}_transcripts"
elif [[ "${MASTER_REFERENCE}" == *transcript* ]]; then
    _transcript_ref="${MASTER_REFERENCE}"
    _genome_ref="${MASTER_REFERENCE%%_transcript*}_genome"
else
    # Bare accession (e.g., "GPE001970"): derive both suffixed forms
    _genome_ref="${MASTER_REFERENCE}_genome"
    _transcript_ref="${MASTER_REFERENCE}_transcripts"
fi
# Auto-detect transcript reference by checking multiple method directories.
# Handles non-standard naming like Eggplant_V4.1_transcripts.function.
_base_pattern="${_genome_ref%_genome}"
[[ "$_base_pattern" == "$_genome_ref" ]] && _base_pattern="$_genome_ref"
_align_base="${ALIGNMENT_BASE:-${BASE_DIR}/II_RESULTS/2_ALIGNMENT_RESULTs}"
_transcript_ref_dirs=(
    "$_align_base/M2_HISAT2_DeNovo/stringtie_WD"
    "$_align_base/M4_Salmon_Saf/Salmon_Quant"
    "$_align_base/M5_RSEM_Bowtie2/RSEM_Quant_WD"
)
# Bash glob replaces find|head subprocess — O(1) vs O(N) directory scan
for _probe_dir in "${_transcript_ref_dirs[@]}"; do
    if [[ -d "$_probe_dir" && ! -d "$_probe_dir/$_transcript_ref" ]]; then
        for _candidate in "$_probe_dir"/${_base_pattern}*transcript*/; do
            [[ -d "$_candidate" ]] || continue
            _found_ref="${_candidate%/}"
            _found_ref="${_found_ref##*/}"
            _parent_dir="${_probe_dir%/*}"
            log_info "Auto-derived transcript ref '${_transcript_ref}' not found in ${_parent_dir##*/}; using '${_found_ref}'"
            _transcript_ref="$_found_ref"
            break 2
        done
    fi
done
unset _base_pattern _found_ref _align_base _transcript_ref_dirs _probe_dir
declare -A METHOD_REF_DIRS
METHOD_REF_DIRS[M1_HISAT2_RefGuided]="${M1_REF_DIR:-$_genome_ref}"
METHOD_REF_DIRS[M2_HISAT2_DeNovo]="${M2_REF_DIR:-$_transcript_ref}"
METHOD_REF_DIRS[M3_STAR_Align]="${M3_REF_DIR:-$_genome_ref}"
METHOD_REF_DIRS[M4_Salmon_Saf]="${M4_REF_DIR:-$_transcript_ref}"
METHOD_REF_DIRS[M5_RSEM_Bowtie2]="${M5_REF_DIR:-$_transcript_ref}"
unset _genome_ref _transcript_ref

# Gene groups for concordance analysis (comma-separated basenames without .csv)
# If a sourced config set GENE_GROUPS as a bash array, join with commas
# (the R concordance config expects comma-separated, not space-separated).
# O(G) array join via IFS — avoids O(G²) string concatenation in a loop
if [[ -v GENE_GROUPS && "${GENE_GROUPS@a}" == *a* ]]; then
    _join_with _gg_joined ',' "${GENE_GROUPS[@]}"
    unset GENE_GROUPS
    GENE_GROUPS="$_gg_joined"
    unset _gg_joined
fi
GENE_GROUPS="${GENE_GROUPS:-SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF}"

# Helper: check if an analysis is enabled — O(1) associative array lookup
declare -A _ANALYSES_SET=()
for _a in "${ANALYSES[@]}"; do _ANALYSES_SET["$_a"]=1; done
unset _a
analysis_enabled() { [[ -n "${_ANALYSES_SET[$1]:-}" ]]; }

# Gene groups directory (reference-specific — strip _genome/_transcripts suffix to match dir name)
_GG_REF_TAG="${MASTER_REFERENCE%%_genome*}"
_GG_REF_TAG="${_GG_REF_TAG%%_transcripts*}"
_GG_REF_DIR="${BASE_DIR}/I_INPUTS/inputs/3_post_proc_inputs/gene_groups_csv/experimental/${_GG_REF_TAG}"
if [[ ! -d "$_GG_REF_DIR" ]]; then
    log_warn "Reference-specific gene groups dir not found: $_GG_REF_DIR"
    _GG_REF_DIR="${BASE_DIR}/I_INPUTS/inputs/3_post_proc_inputs/gene_groups_csv"
    log_warn "Falling back to generic gene groups dir: $_GG_REF_DIR"
fi
GENE_GROUPS_DIR="${GENE_GROUPS_DIR:-$_GG_REF_DIR}"
unset _GG_REF_TAG _GG_REF_DIR

# SRR CSV directory for sample labels
SRR_CSV_DIR="${SRR_CSV_DIR:-${BASE_DIR}/I_INPUTS/inputs/3_post_proc_inputs/SRR_csv}"

# System resources (auto-detect with sane fallbacks)
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-${PBS_NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 12)}}}"
ENABLE_GPU="${ENABLE_GPU:-FALSE}"
if [[ -z "${AVAILABLE_RAM_GB:-}" ]]; then
    if [[ -f /proc/meminfo ]]; then
        # Pure bash: avoids awk fork.  O(1) — reads ~25 lines then breaks.
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
        # Matches run_post_processing.sh pattern.
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
fi
# Fallback covers both unset and empty (e.g., MemAvailable line missing from /proc/meminfo)
[[ -z "${AVAILABLE_RAM_GB:-}" ]] && AVAILABLE_RAM_GB=24
# Ensure numeric (reset to default if non-numeric) before arithmetic
[[ "$AVAILABLE_RAM_GB" =~ ^[0-9]+$ ]] || AVAILABLE_RAM_GB=24
# Floor: ensure at least 4 GB to avoid starving R scripts
(( AVAILABLE_RAM_GB < 4 )) && AVAILABLE_RAM_GB=4
GPU_VRAM_GB="${GPU_VRAM_GB:-8}"

#===============================================================================
# EXPORT ENVIRONMENT FOR R SCRIPTS
#===============================================================================

# Build METHOD_REF_DIRS_STR from associative array for R consumption
# Format: "M1_HISAT2_RefGuided=GPE001970_genome;M2_HISAT2_DeNovo=GPE001970_transcripts;..."
# O(M) array build + single join — avoids O(M²) string concatenation
declare -a _mrd_parts=()
for method in ${METHODS}; do
    ref_dir="${METHOD_REF_DIRS[$method]:-}"
    [[ -n "$ref_dir" ]] && _mrd_parts+=("${method}=${ref_dir}")
done
_join_with METHOD_REF_DIRS_STR ';' "${_mrd_parts[@]}"
unset _mrd_parts

# Concordance mode: cross_method (default), cross_genome, cross_gene_group
CONCORDANCE_MODE="${CONCORDANCE_MODE:-cross_method}"
FIXED_METHOD="${FIXED_METHOD:-}"

# Build CONCORDANCE_GENOMES_STR from array for R (semicolon-separated)
# Must unset array before reassignment — bash arrays cannot be exported to child processes.
# O(G) array join via IFS — avoids O(G²) string concatenation
CONCORDANCE_GENOMES_STR=""
if [[ -v CONCORDANCE_GENOMES && "${CONCORDANCE_GENOMES@a}" == *a* ]]; then
    _join_with CONCORDANCE_GENOMES_STR ';' "${CONCORDANCE_GENOMES[@]}"
    unset CONCORDANCE_GENOMES
fi

export WF_MANAGED_ENV BASE_DIR MASTER_REFERENCE METHODS METHOD_REF_DIRS_STR
export GENE_GROUPS GENE_GROUPS_DIR SRR_CSV_DIR
export THREADS ENABLE_GPU AVAILABLE_RAM_GB GPU_VRAM_GB
export CONCORDANCE_SCRIPT_DIR ANALYSIS_MODULES_DIR UTILITIES_DIR
export OUTPUT_DIR ALIGNMENT_BASE POST_PROC_BASE REPORT_BASE
export FIGURE_DPI
# Sync CONCORDANCE_GENOMES (exported string for R) and CONCORDANCE_GENOMES_STR:
# - Array from TOML: already joined above → assign string form
# - Plain string from env var: if-block was skipped → preserve it and sync _STR
if [[ -n "$CONCORDANCE_GENOMES_STR" ]]; then
    CONCORDANCE_GENOMES="${CONCORDANCE_GENOMES_STR}"
elif [[ -n "${CONCORDANCE_GENOMES:-}" ]]; then
    CONCORDANCE_GENOMES_STR="${CONCORDANCE_GENOMES}"
fi
# Per-genome gene group CSV mapping for cross-genome orthology (positional row correspondence)
# GENOME_GENE_GROUPS_MAP is set by TOML as a flat string; export as-is for R parsing.
export GENOME_GENE_GROUPS_MAP="${GENOME_GENE_GROUPS_MAP:-}"
export CONCORDANCE_MODE FIXED_METHOD CONCORDANCE_GENOMES_STR CONCORDANCE_GENOMES

# Set CURRENT_METHOD for 0_shared_config.R compatibility (avoids warning)
if [[ -z "${CURRENT_METHOD:-}" ]]; then
    if [[ -n "${FIXED_METHOD:-}" ]]; then
        export CURRENT_METHOD="${FIXED_METHOD}"
    elif [[ -n "${METHODS:-}" ]]; then
        CURRENT_METHOD="${METHODS%% *}"
        export CURRENT_METHOD
    fi
fi

#===============================================================================
# RUN ANALYSIS PIPELINE
#===============================================================================

log_step "CONCORDANCE ANALYSIS (mode: ${CONCORDANCE_MODE})"
log_info "Reference:    ${MASTER_REFERENCE}"
[[ -n "$FIXED_METHOD" ]] && log_info "Fixed method: ${FIXED_METHOD}"
log_info "Methods:      ${METHODS}"
log_info "Gene groups:  ${GENE_GROUPS}"
log_info "Analyses:     ${ANALYSES[*]}"
log_info "Output:       ${OUTPUT_DIR}"
log_info "Threads:      ${THREADS}"

# Select loader and concordance scripts based on mode
case "${CONCORDANCE_MODE}" in
    cross_genome)
        STEP1_SCRIPT="1_load_matrices_cross_genome.R"
        STEP2_SCRIPT="2_quantification_concordance.R"
        ;;
    cross_equivalent_gene)
        STEP1_SCRIPT="1_load_matrices_cross_genome.R"
        STEP2_SCRIPT="2_equivalent_gene_concordance.R"
        ;;
    cross_gene_group)
        STEP1_SCRIPT="1_load_matrices_cross_gene_group.R"
        STEP2_SCRIPT="2_gene_group_concordance.R"
        ;;
    cross_method)
        STEP1_SCRIPT="1_load_matrices.R"
        STEP2_SCRIPT="2_quantification_concordance.R"
        ;;
    *)
        log_warn "Unknown CONCORDANCE_MODE '${CONCORDANCE_MODE}'; defaulting to cross_method"
        STEP1_SCRIPT="1_load_matrices.R"
        STEP2_SCRIPT="2_quantification_concordance.R"
        ;;
esac

run_step() {
    local step_num="$1"
    local step_name="$2"
    local script="$3"

    log_step "[STEP ${step_num}/3] ${step_name}"

    if ! Rscript "${CONCORDANCE_SCRIPT_DIR}/${script}"; then
        log_error "Step ${step_num} (${step_name}) failed!"
        exit 1
    fi
}

# ── Batch dispatch: run all enabled steps in a single R session ──
# Saves ~2-3s per additional step by eliminating redundant R interpreter init,
# package loading (ComplexHeatmap ~1.5s), and config sourcing overhead.
# Falls back to individual run_step() when batch dispatcher is missing or
# only 1 step is enabled.
# Big O: reduces O(S × T_startup) to O(T_startup) where S = enabled steps.
_batch_dispatcher="${CONCORDANCE_SCRIPT_DIR}/concordance_batch_dispatcher.R"
_s1_enabled=false; _s2_enabled=false; _s3_enabled=false
analysis_enabled "Load_Matrices"              && _s1_enabled=true
analysis_enabled "Quantification_Concordance" && _s2_enabled=true
analysis_enabled "Generate_Report"            && _s3_enabled=true

# Count enabled steps
_n_enabled=0
$_s1_enabled && (( _n_enabled++ )) || true
$_s2_enabled && (( _n_enabled++ )) || true
$_s3_enabled && (( _n_enabled++ )) || true

if [[ $_n_enabled -ge 2 && -f "$_batch_dispatcher" ]]; then
    # Build batch dispatcher arguments
    _batch_args=()
    $_s1_enabled && _batch_args+=("--step1=${STEP1_SCRIPT}")
    $_s2_enabled && _batch_args+=("--step2=${STEP2_SCRIPT}")
    $_s3_enabled && _batch_args+=("--step3=4_generate_report.R")

    log_step "Concordance batch dispatch (${_n_enabled} steps in single R session)"
    log_info "  Steps: ${_batch_args[*]}"
    if ! Rscript "$_batch_dispatcher" "${_batch_args[@]}"; then
        log_error "Concordance batch dispatch failed!"
        exit 1
    fi
else
    # Sequential fallback: individual Rscript calls per step
    # Step 1 must complete first (produces HARMONIZED_RDS consumed by steps 2-3)
    if $_s1_enabled; then
        run_step 1 "Load & Harmonize Matrices"    "${STEP1_SCRIPT}"
    else
        log_info "Skipping Step 1 (Load_Matrices)"
    fi

    # Check if Step 1 wrote a skip sentinel (e.g., zero common genes in cross_genome mode)
    if [[ -f "${OUTPUT_DIR}/.skip_sentinel" ]]; then
        _skip_reason=$(<"${OUTPUT_DIR}/.skip_sentinel")
        log_info "Analysis skipped by Step 1: ${_skip_reason}"
        log_step "CONCORDANCE ANALYSIS SKIPPED (mode: ${CONCORDANCE_MODE})"
        exit 0
    fi

    # Step 2: Quantification Concordance
    if $_s2_enabled; then
        run_step 2 "Quantification Concordance"   "${STEP2_SCRIPT}"
    else
        log_info "Skipping Step 2 (Quantification_Concordance)"
    fi

    # Step 3: Generate Report
    if $_s3_enabled; then
        run_step 3 "Generate Report"              "4_generate_report.R"
    else
        log_info "Skipping Step 3 (Generate_Report)"
    fi
fi

# Check skip sentinel (may be set by batch dispatcher's Step 1)
if [[ -f "${OUTPUT_DIR}/.skip_sentinel" ]]; then
    _skip_reason=$(<"${OUTPUT_DIR}/.skip_sentinel")
    log_info "Analysis skipped by Step 1: ${_skip_reason}"
    log_step "CONCORDANCE ANALYSIS SKIPPED (mode: ${CONCORDANCE_MODE})"
    exit 0
fi

log_step "CONCORDANCE ANALYSIS COMPLETE (mode: ${CONCORDANCE_MODE})"
# Report filename is always concordance_report.md (generated by 4_generate_report.R)
log_info "Report:  ${REPORT_BASE}/concordance_report.md"
[[ -d "${OUTPUT_DIR}/figures" ]] && log_info "Figures: ${OUTPUT_DIR}/figures/"
[[ -d "${OUTPUT_DIR}/tables"  ]] && log_info "Tables:  ${OUTPUT_DIR}/tables/"

# Copy logs into the output folder for self-contained results.
# In multi-config child mode, only copy this child's own log files (by RUN_ID prefix)
# to avoid a race condition where sibling children concurrently cp -r the shared log dir,
# producing incomplete/stale log snapshots.
if [[ -d "${CONCORDANCE_LOG_BASE}" ]]; then
    mkdir -p "${OUTPUT_DIR}/logs"
    if [[ -n "${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" && -n "${RUN_ID:-}" ]]; then
        # Child mode: copy only files matching this child's RUN_ID to avoid race.
        # Log files live in subdirectories (log_files/, time_logs/, etc.), so iterate
        # over each subdirectory and copy matching files preserving structure.
        # Batch cp per subdirectory: collect matching files, mkdir once, cp once.
        # Saves N-1 mkdir syscalls + reduces N cp forks to 1 per subdir.
        for _subdir in "${CONCORDANCE_LOG_BASE}"/*/; do
            [[ -d "$_subdir" ]] || continue
            _subname="${_subdir%/}"; _subname="${_subname##*/}"
            _matched=()
            for _lf in "$_subdir"*"${RUN_ID}"*; do
                [[ -f "$_lf" ]] && _matched+=("$_lf")
            done
            if [[ ${#_matched[@]} -gt 0 ]]; then
                mkdir -p "${OUTPUT_DIR}/logs/${_subname}"
                cp "${_matched[@]}" "${OUTPUT_DIR}/logs/${_subname}/" 2>/dev/null || true
            fi
        done
        unset _lf _subdir _subname _matched
    else
        # Single-config or parent mode: safe to copy everything
        cp -r "${CONCORDANCE_LOG_BASE}/." "${OUTPUT_DIR}/logs/" 2>/dev/null || true
    fi
    log_info "Logs:    ${OUTPUT_DIR}/logs/"
fi
