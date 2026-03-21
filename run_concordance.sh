#!/bin/bash
#===============================================================================
# CROSS-METHOD CONCORDANCE ANALYSIS
#===============================================================================
# Compares gene expression quantification across all 5 alignment/quantification
# methods (M1-M5) for a given reference genome. Produces correlation matrices,
# discordant gene lists, ranking stability analysis, and a unified report.
#
# Steps:
#   1. Load & harmonize TPM matrices from all methods
#   2. Compute pairwise Spearman/Pearson correlations; identify discordant genes
#   3. Assess ranking stability for gene groups of interest
#   4. Generate unified Markdown report
#
# Usage:
#   bash run_concordance.sh [config_file|config_class|config_class_dir]
#
#   If no config_file is provided, uses internal defaults for GPE001970.
#
#   Config classes can be stored in: config/4_concordance_combination/
#   and selected via:
#     - positional arg: class filename or class basename
#     - env var: CONCORDANCE_CONFIG_CLASSES="class1,class2,..."
#
# Prerequisites:
#   - Alignment results for M1-M5 in 2_ALIGNMENT_RESULTs/
#   - Post-processing matrices in 3_POST_PROC/ (for M4/M5 fallbacks)
#   - R packages: ComplexHeatmap, circlize, grid
#   - conda env "gea" with all dependencies
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}"
REPORT_BASE="${REPORT_BASE:-${BASE_DIR}/4_CONCORDANCE_ANALYSIS}"

#===============================================================================
# CONDA ENVIRONMENT
#===============================================================================

# Skip conda hook if already in the gea environment (~0.3-0.5s saved per invocation)
if [[ "${CONDA_DEFAULT_ENV:-}" != "gea" ]]; then
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
CONFIG_CLASS_DIR="${CONCORDANCE_CONFIG_CLASS_DIR:-${BASE_DIR}/config/4_concordance_combination}"

# Analysis steps to run (comment out entries to skip)
# ─────────────────────────────────────────────────────
if [[ ! "$(declare -p ANALYSES 2>/dev/null)" == "declare -a"* ]]; then
    ANALYSES=(
        "Load_Matrices"                 # Step 1: Load & harmonize matrices (required by all others)
        "Quantification_Concordance"    # Step 2: Compare quantification across methods
        #"Ranking_Stability"             # Step 3: Assess gene ranking stability
        "Generate_Report"               # Step 4: Produce concordance report
    )
fi

# Figure resolution in DPI (300–600)
FIGURE_DPI="${FIGURE_DPI:-300}"

# Clear previous outputs before running
CLEAR_LOGS="TRUE"
CLEAR_OUTPUT_FOLDER="TRUE"

# Optional curated config list (comment in/out as needed).
# Entries can be:
#   - absolute/relative file paths
#   - class basenames from config/4_concordance_combination (with or without .sh)
# Load order matters: later entries override earlier ones.
CONCORDANCE_CONFIGS=(
    #"defaults.toml"
    # "class_genomes_vs_genomes.toml"
    # "class_methods_vs_methods.toml"
    # "class_genes_vs_genes.toml"
    "class_full_factorial_example.toml"
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
    for path in "$entry" "${CONFIG_CLASS_DIR}/${entry}" "${CONFIG_CLASS_DIR}/${entry}.toml" "${CONFIG_CLASS_DIR}/${entry}.sh"; do
        [[ -f "$path" ]] && { source_config_file "$path"; return; }
    done
    log_warn "Config entry not found: ${entry}"
}

# Load manually curated config list first (for easy comment-in/out workflow).
if [[ "$(declare -p CONCORDANCE_CONFIGS 2>/dev/null)" == "declare -a"* && ${#CONCORDANCE_CONFIGS[@]} -gt 0 ]]; then
    for _cfg_entry in "${CONCORDANCE_CONFIGS[@]}"; do
        load_config_entry "$_cfg_entry"
    done
    unset _cfg_entry
fi

# Load optional config file/class/dir (first positional argument)
if [[ -n "$CONFIG_INPUT" ]]; then
    if [[ -f "$CONFIG_INPUT" ]]; then
        source_config_file "$CONFIG_INPUT"
    elif [[ -d "$CONFIG_INPUT" ]]; then
        while IFS= read -r _cfg; do
            source_config_file "$_cfg"
        done < <(find "$CONFIG_INPUT" -maxdepth 1 -type f \( -name "*.toml" -o -name "*.sh" \) | sort)
    elif [[ -f "${CONFIG_CLASS_DIR}/${CONFIG_INPUT}" ]]; then
        source_config_file "${CONFIG_CLASS_DIR}/${CONFIG_INPUT}"
    elif [[ -f "${CONFIG_CLASS_DIR}/${CONFIG_INPUT}.toml" ]]; then
        source_config_file "${CONFIG_CLASS_DIR}/${CONFIG_INPUT}.toml"
    elif [[ -f "${CONFIG_CLASS_DIR}/${CONFIG_INPUT}.sh" ]]; then
        source_config_file "${CONFIG_CLASS_DIR}/${CONFIG_INPUT}.sh"
    else
        log_warn "Config input not found: ${CONFIG_INPUT} (continuing with defaults)"
    fi
fi

# Optionally source additional config classes from config/4_concordance_combination
# Example: CONCORDANCE_CONFIG_CLASSES="defaults,class_methods_vs_methods"
if [[ -n "${CONCORDANCE_CONFIG_CLASSES:-}" ]]; then
    IFS=',' read -r -a _cfg_classes <<< "$CONCORDANCE_CONFIG_CLASSES"
    for _cfg_class in "${_cfg_classes[@]}"; do
        load_config_entry "$_cfg_class"
    done
    unset _cfg_classes _cfg_class
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
        # Wait for any one child to finish, then compact the PID array
        local _still_running=()
        for _p in "${_pids_ref[@]}"; do
            if kill -0 "$_p" 2>/dev/null; then
                _still_running+=("$_p")
            else
                wait "$_p" 2>/dev/null || true
            fi
        done
        _pids_ref=("${_still_running[@]}")
        # If still at capacity, sleep briefly to avoid busy-waiting
        [[ ${#_pids_ref[@]} -ge $MAX_CONCORDANCE_JOBS ]] && sleep 0.5
    done
}

# Optional multi-reference mode:
# - RUN_ALL_MASTER_REFERENCES=TRUE: iterate over MASTER_REFERENCES and run once per reference
# - Default FALSE: run only one reference (MASTER_REFERENCE or first MASTER_REFERENCES entry)
RUN_ALL_MASTER_REFERENCES="${RUN_ALL_MASTER_REFERENCES:-FALSE}"
if [[ "${RUN_ALL_MASTER_REFERENCES^^}" == "TRUE" ]]; then
    if [[ "$(declare -p MASTER_REFERENCES 2>/dev/null)" == "declare -a"* && ${#MASTER_REFERENCES[@]} -gt 0 ]]; then
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
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/${_ref}" \
            __CONCORDANCE_OVERRIDE_MASTER_REFERENCE="${_ref}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_MASTER_REFERENCES="FALSE" \
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
    if [[ "$(declare -p METHOD_COMBINATIONS 2>/dev/null)" == "declare -a"* && ${#METHOD_COMBINATIONS[@]} -gt 0 ]]; then
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
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/methods_${_combo_tag}" \
            __CONCORDANCE_OVERRIDE_METHODS="${_methods}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_METHOD_COMBINATIONS="FALSE" \
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
    if [[ "$(declare -p GENE_GROUP_COMBINATIONS 2>/dev/null)" == "declare -a"* && ${#GENE_GROUP_COMBINATIONS[@]} -gt 0 ]]; then
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
            __CONCORDANCE_OVERRIDE_REPORT_BASE="${_parent_report_base}/genes_${_combo_tag}" \
            __CONCORDANCE_OVERRIDE_GENE_GROUPS="${_groups}" \
            __CONCORDANCE_OVERRIDE_RUN_ALL_GENE_GROUP_COMBINATIONS="FALSE" \
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
if [[ -z "${MASTER_REFERENCE:-}" && "$(declare -p MASTER_REFERENCES 2>/dev/null)" == "declare -a"* ]]; then
    MASTER_REFERENCE="${MASTER_REFERENCES[0]}"
fi
# Default must match 0_concordance_config.R and 0_shared_config.R ("Eggplant_V4.1")
MASTER_REFERENCE="${MASTER_REFERENCE:-Eggplant_V4.1}"

# Methods to compare (space-separated string).
# If a sourced config set METHODS as an array, flatten it to a string.
if [[ "$(declare -p METHODS 2>/dev/null)" == "declare -a"* ]]; then
    METHODS="${METHODS[*]}"
fi
METHODS="${METHODS:-
    M1_HISAT2_RefGuided 
    M2_HISAT2_DeNovo 
    M3_STAR_Align 
    M4_Salmon_Saf 
    M5_RSEM_Bowtie2
}"

# Enforce method/reference compatibility:
# - *_genome references:      M1, M3
# - *transcript* references:  M2, M4, M5
# Set ENFORCE_REFERENCE_METHOD_COMPATIBILITY=FALSE to disable filtering.
ENFORCE_REFERENCE_METHOD_COMPATIBILITY="${ENFORCE_REFERENCE_METHOD_COMPATIBILITY:-TRUE}"
if [[ "${ENFORCE_REFERENCE_METHOD_COMPATIBILITY^^}" == "TRUE" ]]; then
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
_genome_ref="${MASTER_REFERENCE}"
# Build transcript reference by swapping _genome → _transcripts (or appending _transcripts)
if [[ "$_genome_ref" == *_genome ]]; then
    _transcript_ref="${_genome_ref%_genome}_transcripts"
else
    _transcript_ref="${_genome_ref}_transcripts"
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
            log_info "Auto-derived transcript ref '${_transcript_ref}' not found in $(basename "$(dirname "$_probe_dir")"); using '${_found_ref}'"
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

# Gene groups for ranking stability (comma-separated basenames without .csv)
# If a sourced config set GENE_GROUPS as a bash array, join with commas
# (the R concordance config expects comma-separated, not space-separated).
if [[ "$(declare -p GENE_GROUPS 2>/dev/null)" == "declare -a"* ]]; then
    _gg_joined=""
    for _gg in "${GENE_GROUPS[@]}"; do
        _gg_joined="${_gg_joined:+${_gg_joined},}${_gg}"
    done
    GENE_GROUPS="$_gg_joined"
    unset _gg_joined _gg
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
# PATHS
#===============================================================================

REPORT_BASE="${REPORT_BASE:-${BASE_DIR}/4_CONCORDANCE_ANALYSIS}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPORT_BASE}/cross_method_concordance}"
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
if [[ "${CLEAR_LOGS:-FALSE}" == "TRUE" && -d "${REPORT_BASE}/logs" ]]; then
    log_info "Clearing previous logs: ${REPORT_BASE}/logs"
    rm -rf "${REPORT_BASE}/logs"
fi

mkdir -p "${OUTPUT_DIR}/figures" "${OUTPUT_DIR}/tables" || {
    log_error "Failed to create output directories under ${OUTPUT_DIR}"
    exit 1
}

# Set up structured logging (mirrors run_post_processing.sh)
LOG_DIR="${REPORT_BASE}/logs/log_files"
TIME_DIR="${REPORT_BASE}/logs/time_logs"
SPACE_DIR="${REPORT_BASE}/logs/space_logs"
SPACE_TIME_DIR="${REPORT_BASE}/logs/space_time_logs"
ERROR_WARN_DIR="${REPORT_BASE}/logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="${REPORT_BASE}/logs/software_catalogs"
GPU_LOG_DIR="${REPORT_BASE}/logs/gpu_log"
export LOG_DIR TIME_DIR SPACE_DIR SPACE_TIME_DIR ERROR_WARN_DIR SOFTWARE_CATALOG_DIR GPU_LOG_DIR

if declare -f setup_logging &>/dev/null; then
    setup_logging "$CLEAR_LOGS"
    export LOG_FILE TIME_FILE SPACE_FILE SPACE_TIME_FILE ERROR_WARN_FILE SOFTWARE_FILE GPU_LOG_FILE
fi

#===============================================================================
# EXPORT ENVIRONMENT FOR R SCRIPTS
#===============================================================================

# Build METHOD_REF_DIRS_STR from associative array for R consumption
# Format: "M1_HISAT2_RefGuided=GPE001970_genome;M2_HISAT2_DeNovo=GPE001970_transcripts;..."
METHOD_REF_DIRS_STR=""
for method in ${METHODS}; do
    ref_dir="${METHOD_REF_DIRS[$method]:-}"
    if [[ -n "$ref_dir" ]]; then
        METHOD_REF_DIRS_STR="${METHOD_REF_DIRS_STR:+${METHOD_REF_DIRS_STR};}${method}=${ref_dir}"
    fi
done

export BASE_DIR MASTER_REFERENCE METHODS METHOD_REF_DIRS_STR
export GENE_GROUPS GENE_GROUPS_DIR SRR_CSV_DIR
export THREADS ENABLE_GPU AVAILABLE_RAM_GB GPU_VRAM_GB
export CONCORDANCE_SCRIPT_DIR ANALYSIS_MODULES_DIR UTILITIES_DIR
export OUTPUT_DIR ALIGNMENT_BASE POST_PROC_BASE REPORT_BASE
export FIGURE_DPI

#===============================================================================
# RUN ANALYSIS PIPELINE
#===============================================================================

log_step "CROSS-METHOD CONCORDANCE ANALYSIS"
log_info "Reference:    ${MASTER_REFERENCE}"
log_info "Methods:      ${METHODS}"
log_info "Gene groups:  ${GENE_GROUPS}"
log_info "Analyses:     ${ANALYSES[*]}"
log_info "Output:       ${OUTPUT_DIR}"
log_info "Threads:      ${THREADS}"

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
    run_step 1 "Load & Harmonize Matrices"    "1_load_matrices.R"
else
    log_info "Skipping Step 1 (Load_Matrices)"
fi

# Steps 2 and 3 are independent (both read HARMONIZED_RDS, write separate outputs).
# Run enabled steps in parallel to halve wall-clock time for this phase.
_run_step2=false; _run_step3=false
analysis_enabled "Quantification_Concordance" && _run_step2=true
analysis_enabled "Ranking_Stability"          && _run_step3=true

if [[ "$_run_step2" == true || "$_run_step3" == true ]]; then
    _label=""
    $_run_step2 && _label="Quantification Concordance"
    $_run_step3 && _label="${_label:+${_label} & }Ranking Stability"
    _parallel_note=""
    [[ "$_run_step2" == true && "$_run_step3" == true ]] && _parallel_note=" (parallel)"
    log_step "[STEPS 2+3] ${_label}${_parallel_note}"

    _step2_log="${OUTPUT_DIR}/step2.log"
    _step3_log="${OUTPUT_DIR}/step3.log"

    # Ensure temp logs are cleaned up on early exit (SIGINT/SIGTERM)
    _concordance_cleanup() { rm -f "$_step2_log" "$_step3_log"; }
    trap '_concordance_cleanup' EXIT

    _pid2=""; _pid3=""
    if [[ "$_run_step2" == true ]]; then
        Rscript "${CONCORDANCE_SCRIPT_DIR}/2_quantification_concordance.R" > "$_step2_log" 2>&1 &
        _pid2=$!
    fi
    if [[ "$_run_step3" == true ]]; then
        Rscript "${CONCORDANCE_SCRIPT_DIR}/3_ranking_stability.R" > "$_step3_log" 2>&1 &
        _pid3=$!
    fi

    _step2_rc=0; _step3_rc=0
    [[ -n "$_pid2" ]] && { wait "$_pid2" || _step2_rc=$?; }
    [[ -n "$_pid3" ]] && { wait "$_pid3" || _step3_rc=$?; }

    # Stream logs to stdout for visibility
    cat "$_step2_log" "$_step3_log" 2>/dev/null
    rm -f "$_step2_log" "$_step3_log"

    _any_failed=0
    if [[ $_step2_rc -ne 0 ]]; then
        log_error "Step 2 (Quantification Concordance) failed (exit=$_step2_rc)!"
        _any_failed=1
    fi
    if [[ $_step3_rc -ne 0 ]]; then
        log_error "Step 3 (Ranking Stability) failed (exit=$_step3_rc)!"
        _any_failed=1
    fi
    [[ $_any_failed -ne 0 ]] && exit 1
    log_info "Steps 2+3 completed successfully"
else
    log_info "Skipping Steps 2+3 (Quantification_Concordance, Ranking_Stability)"
fi

# Step 4 reads outputs from both steps 2 and 3
if analysis_enabled "Generate_Report"; then
    run_step 4 "Generate Report"              "4_generate_report.R"
else
    log_info "Skipping Step 4 (Generate_Report)"
fi

log_step "CONCORDANCE ANALYSIS COMPLETE"
log_info "Report:  ${REPORT_BASE}/cross_method_concordance_report.md"
log_info "Figures: ${OUTPUT_DIR}/figures/"
log_info "Tables:  ${OUTPUT_DIR}/tables/"

# Copy logs into the output folder for self-contained results
if [[ -d "${REPORT_BASE}/logs" ]]; then
    mkdir -p "${OUTPUT_DIR}/logs"
    cp -r "${REPORT_BASE}/logs/." "${OUTPUT_DIR}/logs/" 2>/dev/null || true
    log_info "Logs:    ${OUTPUT_DIR}/logs/"
fi
