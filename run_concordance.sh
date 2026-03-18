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
#   bash run_cross_method_concordance.sh [config_file]
#
#   If no config_file is provided, uses internal defaults for GPE001970.
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

#===============================================================================
# CONDA ENVIRONMENT
#===============================================================================

eval "$(conda shell.bash hook 2>/dev/null)" 2>/dev/null || true
conda activate gea 2>/dev/null || true

# Source logging utilities for consistent output
source "${SCRIPT_DIR}/modules/logging/logging_utils.sh" 2>/dev/null || {
    # Minimal fallback if logging module unavailable
    log_info()  { echo "[INFO]  $*"; }
    log_warn()  { echo "[WARN]  $*"; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_step()  { echo ""; echo "==> $*"; }
}

#===============================================================================
# CONFIGURATION (override via env, config file, or associative array below)
#===============================================================================

# Load optional config file (first positional argument)
if [[ -n "${1:-}" && -f "$1" ]]; then
    log_info "Loading config from: $1"
    source "$1"
fi

# Reference genome to compare across methods
# If a sourced config set MASTER_REFERENCES as an array, take the first element
# (matches run_post_processing.sh behaviour).
if [[ "$(declare -p MASTER_REFERENCES 2>/dev/null)" == "declare -a"* ]]; then
    MASTER_REFERENCE="${MASTER_REFERENCES[0]}"
fi
MASTER_REFERENCE="${MASTER_REFERENCE:-GPE001970_genome}"

# Methods to compare (space-separated string).
# If a sourced config set METHODS as an array, flatten it to a string.
if [[ "$(declare -p METHODS 2>/dev/null)" == "declare -a"* ]]; then
    METHODS="${METHODS[*]}"
fi
METHODS="${METHODS:-M1_HISAT2_RefGuided M2_HISAT2_DeNovo M3_STAR_Align M4_Salmon_Saf M5_RSEM_Bowtie2}"

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
# Verify the auto-derived transcript reference exists for at least one method's
# alignment directory.  If not, search for the actual transcript reference dir
# (handles non-standard naming like Eggplant_V4.1_transcripts.function).
_m2_stringtie_base="${ALIGNMENT_BASE:-${BASE_DIR}/2_ALIGNMENT_RESULTs}/M2_HISAT2_DeNovo/stringtie_WD"
if [[ -d "$_m2_stringtie_base" && ! -d "$_m2_stringtie_base/$_transcript_ref" ]]; then
    _base_pattern="${_genome_ref%_genome}"
    [[ "$_base_pattern" == "$_genome_ref" ]] && _base_pattern="$_genome_ref"
    _found_ref=$(find "$_m2_stringtie_base" -maxdepth 1 -type d -name "${_base_pattern}*transcript*" -printf '%f\n' 2>/dev/null | head -1)
    if [[ -n "$_found_ref" ]]; then
        log_info "Auto-derived transcript ref '${_transcript_ref}' not found; using '${_found_ref}'"
        _transcript_ref="$_found_ref"
    fi
    unset _base_pattern _found_ref
fi
unset _m2_stringtie_base
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

# Gene groups directory (reference-specific — strip _genome/_transcripts suffix to match dir name)
_GG_REF_TAG="${MASTER_REFERENCE%%_genome}"
_GG_REF_TAG="${_GG_REF_TAG%%_transcripts}"
_GG_REF_DIR="${BASE_DIR}/inputs/gene_groups_csv/experimental/${_GG_REF_TAG}"
if [[ ! -d "$_GG_REF_DIR" ]]; then
    log_warn "Reference-specific gene groups dir not found: $_GG_REF_DIR"
    _GG_REF_DIR="${BASE_DIR}/inputs/gene_groups_csv"
    log_warn "Falling back to generic gene groups dir: $_GG_REF_DIR"
fi
GENE_GROUPS_DIR="${GENE_GROUPS_DIR:-$_GG_REF_DIR}"
unset _GG_REF_TAG _GG_REF_DIR

# SRR CSV directory for sample labels
SRR_CSV_DIR="${SRR_CSV_DIR:-${BASE_DIR}/inputs/SRR_csv}"

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
AVAILABLE_RAM_GB="${AVAILABLE_RAM_GB:-24}"
GPU_VRAM_GB="${GPU_VRAM_GB:-8}"

#===============================================================================
# PATHS
#===============================================================================

OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/3_POST_PROC/cross_method_concordance}"
ALIGNMENT_BASE="${ALIGNMENT_BASE:-${BASE_DIR}/2_ALIGNMENT_RESULTs}"
POST_PROC_BASE="${POST_PROC_BASE:-${BASE_DIR}/3_POST_PROC}"
ANALYSIS_MODULES_DIR="${BASE_DIR}/modules/c_post_processing/analysis_modules"
UTILITIES_DIR="${BASE_DIR}/modules/c_post_processing/utilities"
CONCORDANCE_SCRIPT_DIR="${BASE_DIR}/modules/c_post_processing/cross_method_concordance"

mkdir -p "${OUTPUT_DIR}/figures" "${OUTPUT_DIR}/tables" || {
    log_error "Failed to create output directories under ${OUTPUT_DIR}"
    exit 1
}

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
export OUTPUT_DIR ALIGNMENT_BASE POST_PROC_BASE

#===============================================================================
# RUN ANALYSIS PIPELINE
#===============================================================================

log_step "CROSS-METHOD CONCORDANCE ANALYSIS"
log_info "Reference:    ${MASTER_REFERENCE}"
log_info "Methods:      ${METHODS}"
log_info "Gene groups:  ${GENE_GROUPS}"
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
run_step 1 "Load & Harmonize Matrices"    "1_load_matrices.R"

# Steps 2 and 3 are independent (both read HARMONIZED_RDS, write separate outputs).
# Run them in parallel to halve wall-clock time for this phase.
log_step "[STEPS 2+3] Quantification Concordance & Ranking Stability (parallel)"
_step2_log="${OUTPUT_DIR}/step2.log"
_step3_log="${OUTPUT_DIR}/step3.log"

Rscript "${CONCORDANCE_SCRIPT_DIR}/2_quantification_concordance.R" > "$_step2_log" 2>&1 &
_pid2=$!
Rscript "${CONCORDANCE_SCRIPT_DIR}/3_ranking_stability.R" > "$_step3_log" 2>&1 &
_pid3=$!

_step2_rc=0; _step3_rc=0
wait "$_pid2" || _step2_rc=$?
wait "$_pid3" || _step3_rc=$?

# Stream logs to stdout for visibility
cat "$_step2_log" "$_step3_log" 2>/dev/null
rm -f "$_step2_log" "$_step3_log"

if [[ $_step2_rc -ne 0 ]]; then
    log_error "Step 2 (Quantification Concordance) failed (exit=$_step2_rc)!"
    exit 1
fi
if [[ $_step3_rc -ne 0 ]]; then
    log_error "Step 3 (Ranking Stability) failed (exit=$_step3_rc)!"
    exit 1
fi
log_info "Steps 2 and 3 completed successfully"

# Step 4 reads outputs from both steps 2 and 3
run_step 4 "Generate Report"              "4_generate_report.R"

log_step "CONCORDANCE ANALYSIS COMPLETE"
log_info "Report:  ${POST_PROC_BASE}/cross_method_concordance_report.md"
log_info "Figures: ${OUTPUT_DIR}/figures/"
log_info "Tables:  ${OUTPUT_DIR}/tables/"
