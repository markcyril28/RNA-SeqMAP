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
#   2. Compute pairwise Spearman/Pearson correlations; identify discordant genes
#   3. Generate unified Markdown report
#
# Usage:
#   bash run_concordance.sh [config_file|config_cross|config_cross_dir]
#
#   If no config_file is provided, uses internal defaults for GPE001970.
#
#   Config crosses can be stored in: config/4_concordance_combination/
#   and selected via:
#     - positional arg: cross filename or cross basename
#     - env var: CONCORDANCE_CONFIG_CROSSES="cross1,cross2,..."
#
# Prerequisites:
#   - Alignment results for M1-M5 in 2_ALIGNMENT_RESULTs/
#   - Post-processing matrices in 3_POST_PROC/ (for M4/M5 fallbacks)
#   - R packages: ComplexHeatmap, circlize, grid
#   - conda env "gea" with all dependencies
#===============================================================================

set -euo pipefail

# Resolve script directory: parameter expansion avoids nested $(dirname) subshell
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"
BASE_DIR="${SCRIPT_DIR}"
REPORT_BASE="${REPORT_BASE:-${BASE_DIR}/4_CONCORDANCE_ANALYSIS}"

#===============================================================================
# CONDA ENVIRONMENT
#===============================================================================

# Activate conda env only if not already active — child processes (bash "$0")
# inherit CONDA_DEFAULT_ENV but not the conda PATH entries from the parent shell,
# so we must always run the hook. Skip only when PATH is already configured.
# O(1) string match avoids ~0.3-0.5s conda hook overhead per child process.
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

# Source logging utilities for consistent output
source "${SCRIPT_DIR}/modules/logging/logging_utils.sh" 2>/dev/null || {
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
CONFIG_CROSS_DIR="${CONCORDANCE_CONFIG_CROSS_DIR:-${BASE_DIR}/config/4_concordance_combination}"

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
        [[ -f "$path" ]] && { source_config_file "$path"; return; }
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
        source_config_file "$CONFIG_INPUT"
    elif [[ -d "$CONFIG_INPUT" ]]; then
        # Bash glob + mapfile avoids find+sort subprocess pair. O(F) where F = config files.
        _cfg_files=()
        for _cfg in "$CONFIG_INPUT"/*.toml "$CONFIG_INPUT"/*.sh; do
            [[ -f "$_cfg" ]] && _cfg_files+=("$_cfg")
        done
        # Sort for deterministic order (globs are locale-sorted but toml/sh interleave)
        mapfile -t _cfg_files < <(printf '%s\n' "${_cfg_files[@]}" | sort)
        for _cfg in "${_cfg_files[@]}"; do
            source_config_file "$_cfg"
        done
        unset _cfg_files
    elif [[ -f "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}" ]]; then
        source_config_file "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}"
    elif [[ -f "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.toml" ]]; then
        source_config_file "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.toml"
    elif [[ -f "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.sh" ]]; then
        source_config_file "${CONFIG_CROSS_DIR}/${CONFIG_INPUT}.sh"
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
    # Pure bash: replace separators with _, strip non-alnum (avoids 3 subshell spawns)
    local clean="${1//[,|\/[:space:]]/_}"
    clean="${clean//[^[:alnum:]_.-]/}"
    [[ -n "$clean" ]] && printf "%s" "$clean" || printf "combo"
}

#===============================================================================
# PATHS (must be set before dispatch blocks so parent logging works)
#===============================================================================

REPORT_BASE="${REPORT_BASE:-${BASE_DIR}/4_CONCORDANCE_ANALYSIS}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPORT_BASE}}"
ALIGNMENT_BASE="${ALIGNMENT_BASE:-${BASE_DIR}/2_ALIGNMENT_RESULTs}"
POST_PROC_BASE="${POST_PROC_BASE:-${BASE_DIR}/3_POST_PROC}"
ANALYSIS_MODULES_DIR="${BASE_DIR}/modules/c_post_processing/analysis_modules"
UTILITIES_DIR="${BASE_DIR}/modules/c_post_processing/utilities"
CONCORDANCE_SCRIPT_DIR="${BASE_DIR}/modules/c_post_processing/cross_method_concordance"

# Clear previous outputs if requested
if [[ "${CLEAR_OUTPUT_FOLDER:-FALSE}" == "TRUE" && -d "${OUTPUT_DIR}" ]]; then
    log_info "Clearing previous output folder: ${OUTPUT_DIR}"
    rm -rf "${OUTPUT_DIR}"
fi
# Logs always live under the top-level REPORT_BASE (4_CONCORDANCE_ANALYSIS/logs/),
# even when child processes override REPORT_BASE to per-config subdirectories.
CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE:-${REPORT_BASE}/logs}"
if [[ "${CLEAR_LOGS:-FALSE}" == "TRUE" && -d "${CONCORDANCE_LOG_BASE}" ]]; then
    log_info "Clearing previous logs: ${CONCORDANCE_LOG_BASE}"
    rm -rf "${CONCORDANCE_LOG_BASE}"
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

# Throttle helper: waits until the number of tracked PIDs drops below MAX_CONCORDANCE_JOBS.
# Usage: _throttle_pids <array_name>
# Complexity: O(J) per call where J = active jobs; total O(N) amortized across N launches.
_throttle_pids() {
    local -n _pids_ref=$1
    while [[ ${#_pids_ref[@]} -ge $MAX_CONCORDANCE_JOBS ]]; do
        # Use `wait -n` (Bash 4.3+) to block until any child exits — avoids busy-wait polling.
        # Falls back to poll+sleep for older Bash versions.
        if wait -n "${_pids_ref[@]}" 2>/dev/null; then true; fi
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
        CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE \
        bash "$0" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
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
        _parent_report_base="${REPORT_BASE:-${BASE_DIR}/4_CONCORDANCE_ANALYSIS}"
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
            CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE \
            bash "$0" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
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
        _parent_report_base="${REPORT_BASE:-${BASE_DIR}/4_CONCORDANCE_ANALYSIS}"
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
            _combo_tag="$(sanitize_tag "$_combo")"
            log_step "Launching concordance for methods: ${_methods}"
            printf -v _child_run_id '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || _child_run_id=$(date +%Y%m%d_%H%M%S)
            _child_run_id="${_child_run_id}_methods_${_combo_tag}"
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/methods_${_combo_tag}" \
            __CONCORDANCE_OVERRIDE_METHODS="${_methods}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_METHOD_COMBINATIONS="FALSE" \
            __CONCORDANCE_OVERRIDE_SINGLE_CONFIG="${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" \
            CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE}" \
            RUN_ID="${_child_run_id}" LOGGING_INITIALIZED="" \
            CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE \
            bash "$0" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
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
        _parent_report_base="${REPORT_BASE:-${BASE_DIR}/4_CONCORDANCE_ANALYSIS}"
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
            _combo_tag="$(sanitize_tag "$_combo")"
            log_step "Launching concordance for gene groups: ${_groups}"
            printf -v _child_run_id '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || _child_run_id=$(date +%Y%m%d_%H%M%S)
            _child_run_id="${_child_run_id}_genes_${_combo_tag}"
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/genes_${_combo_tag}" \
            __CONCORDANCE_OVERRIDE_GENE_GROUPS="${_groups}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_GENE_GROUP_COMBINATIONS="FALSE" \
            __CONCORDANCE_OVERRIDE_SINGLE_CONFIG="${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" \
            CONCORDANCE_LOG_BASE="${CONCORDANCE_LOG_BASE}" \
            RUN_ID="${_child_run_id}" LOGGING_INITIALIZED="" \
            CLEAR_LOGS=FALSE CLEAR_OUTPUT_FOLDER=FALSE \
            bash "$0" ${REINVOKE_ARGS[@]+"${REINVOKE_ARGS[@]}"} &
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
_align_base="${ALIGNMENT_BASE:-${BASE_DIR}/2_ALIGNMENT_RESULTs}"
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
    _gg_joined="$(IFS=','; printf '%s' "${GENE_GROUPS[*]}")"
    unset GENE_GROUPS
    GENE_GROUPS="$_gg_joined"
    unset _gg_joined
fi
GENE_GROUPS="${GENE_GROUPS:-SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF}"

# Helper: check if an analysis is enabled
analysis_enabled() {
    local target="$1"
    local a
    for a in "${ANALYSES[@]}"; do
        [[ "$a" == "$target" ]] && return 0
    done
    return 1
}

# Gene groups directory (reference-specific — strip _genome/_transcripts suffix to match dir name)
_GG_REF_TAG="${MASTER_REFERENCE%%_genome*}"
_GG_REF_TAG="${_GG_REF_TAG%%_transcripts*}"
_GG_REF_DIR="${BASE_DIR}/inputs/3_post_proc_inputs/gene_groups_csv/experimental/${_GG_REF_TAG}"
if [[ ! -d "$_GG_REF_DIR" ]]; then
    log_warn "Reference-specific gene groups dir not found: $_GG_REF_DIR"
    _GG_REF_DIR="${BASE_DIR}/inputs/3_post_proc_inputs/gene_groups_csv"
    log_warn "Falling back to generic gene groups dir: $_GG_REF_DIR"
fi
GENE_GROUPS_DIR="${GENE_GROUPS_DIR:-$_GG_REF_DIR}"
unset _GG_REF_TAG _GG_REF_DIR

# SRR CSV directory for sample labels
SRR_CSV_DIR="${SRR_CSV_DIR:-${BASE_DIR}/inputs/3_post_proc_inputs/SRR_csv}"

# System resources (auto-detect with sane fallbacks)
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 12)}"
ENABLE_GPU="${ENABLE_GPU:-FALSE}"
if [[ -z "${AVAILABLE_RAM_GB:-}" ]]; then
    if [[ -f /proc/meminfo ]]; then
        AVAILABLE_RAM_GB=$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)
    elif command -v sysctl &>/dev/null; then
        AVAILABLE_RAM_GB=$(sysctl -n hw.memsize 2>/dev/null | awk '{printf "%d", $1/1073741824}')
    fi
fi
# Fallback covers both unset and empty (e.g., MemAvailable line missing from /proc/meminfo)
[[ -z "${AVAILABLE_RAM_GB:-}" ]] && AVAILABLE_RAM_GB=24
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
METHOD_REF_DIRS_STR="$(IFS=';'; printf '%s' "${_mrd_parts[*]}")"
unset _mrd_parts

# Concordance mode: cross_method (default), cross_genome, cross_gene_group
CONCORDANCE_MODE="${CONCORDANCE_MODE:-cross_method}"
FIXED_METHOD="${FIXED_METHOD:-}"

# Build CONCORDANCE_GENOMES_STR from array for R (semicolon-separated)
# Must unset array before reassignment — bash arrays cannot be exported to child processes.
# O(G) array join via IFS — avoids O(G²) string concatenation
CONCORDANCE_GENOMES_STR=""
if [[ -v CONCORDANCE_GENOMES && "${CONCORDANCE_GENOMES@a}" == *a* ]]; then
    CONCORDANCE_GENOMES_STR="$(IFS=';'; printf '%s' "${CONCORDANCE_GENOMES[*]}")"
    unset CONCORDANCE_GENOMES
fi

export BASE_DIR MASTER_REFERENCE METHODS METHOD_REF_DIRS_STR
export GENE_GROUPS GENE_GROUPS_DIR SRR_CSV_DIR
export THREADS ENABLE_GPU AVAILABLE_RAM_GB GPU_VRAM_GB
export CONCORDANCE_SCRIPT_DIR ANALYSIS_MODULES_DIR UTILITIES_DIR
export OUTPUT_DIR ALIGNMENT_BASE POST_PROC_BASE REPORT_BASE
export FIGURE_DPI
CONCORDANCE_GENOMES="${CONCORDANCE_GENOMES_STR}"
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

    log_step "[STEP ${step_num}/4] ${step_name}"

    if ! Rscript "${CONCORDANCE_SCRIPT_DIR}/${script}"; then
        log_error "Step ${step_num} (${step_name}) failed!"
        exit 1
    fi
}

# Step 1 must complete first (produces HARMONIZED_RDS consumed by steps 2-4)
if analysis_enabled "Load_Matrices"; then
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
if analysis_enabled "Quantification_Concordance"; then
    run_step 2 "Quantification Concordance"   "${STEP2_SCRIPT}"
else
    log_info "Skipping Step 2 (Quantification_Concordance)"
fi

# Step 3: Generate Report
if analysis_enabled "Generate_Report"; then
    run_step 3 "Generate Report"              "4_generate_report.R"
else
    log_info "Skipping Step 3 (Generate_Report)"
fi

log_step "CONCORDANCE ANALYSIS COMPLETE (mode: ${CONCORDANCE_MODE})"
# Report filename matches R prefix: cross_method, cross_genome, or cross_gene_group
_report_prefix="${CONCORDANCE_MODE:-cross_method}"
log_info "Report:  ${REPORT_BASE}/${_report_prefix}_concordance_report.md"
[[ -d "${OUTPUT_DIR}/figures" ]] && log_info "Figures: ${OUTPUT_DIR}/figures/"
[[ -d "${OUTPUT_DIR}/tables"  ]] && log_info "Tables:  ${OUTPUT_DIR}/tables/"

# Copy logs into the output folder for self-contained results.
# In multi-config child mode, only copy this child's own log files (by RUN_ID prefix)
# to avoid a race condition where sibling children concurrently cp -r the shared log dir,
# producing incomplete/stale log snapshots.
if [[ -d "${CONCORDANCE_LOG_BASE}" ]]; then
    mkdir -p "${OUTPUT_DIR}/logs"
    if [[ -n "${__CONCORDANCE_OVERRIDE_SINGLE_CONFIG:-}" && -n "${RUN_ID:-}" ]]; then
        # Child mode: copy only files matching this child's RUN_ID to avoid race
        for _lf in "${CONCORDANCE_LOG_BASE}"/*"${RUN_ID}"*; do
            [[ -f "$_lf" ]] && cp "$_lf" "${OUTPUT_DIR}/logs/" 2>/dev/null || true
        done
        unset _lf
    else
        # Single-config or parent mode: safe to copy everything
        cp -r "${CONCORDANCE_LOG_BASE}/." "${OUTPUT_DIR}/logs/" 2>/dev/null || true
    fi
    log_info "Logs:    ${OUTPUT_DIR}/logs/"
fi
