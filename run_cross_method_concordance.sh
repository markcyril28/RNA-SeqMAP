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

eval "$(conda shell.bash hook)"
conda activate gea 2>/dev/null || echo "Warning: conda env 'gea' not found, using current env"

#===============================================================================
# CONFIGURATION (override via env, config file, or associative array below)
#===============================================================================

# Load optional config file (first positional argument)
if [[ -n "${1:-}" && -f "$1" ]]; then
    echo "[CONFIG] Loading config from: $1"
    source "$1"
fi

# Reference genome to compare across methods
MASTER_REFERENCE="${MASTER_REFERENCE:-GPE001970_genome}"

# Methods to compare (space-separated)
METHODS="${METHODS:-M1_HISAT2_RefGuided M2_HISAT2_DeNovo M3_STAR_Align M4_Salmon_Saf M5_RSEM_Bowtie2}"

# Method-specific reference directory names
# M1/M3 align to genome; M2/M4/M5 align to transcriptome
declare -A METHOD_REF_DIRS
METHOD_REF_DIRS[M1_HISAT2_RefGuided]="${M1_REF_DIR:-GPE001970_genome}"
METHOD_REF_DIRS[M2_HISAT2_DeNovo]="${M2_REF_DIR:-GPE001970_transcripts}"
METHOD_REF_DIRS[M3_STAR_Align]="${M3_REF_DIR:-GPE001970_genome}"
METHOD_REF_DIRS[M4_Salmon_Saf]="${M4_REF_DIR:-GPE001970_transcripts}"
METHOD_REF_DIRS[M5_RSEM_Bowtie2]="${M5_REF_DIR:-GPE001970_transcripts}"

# Gene groups for ranking stability (comma-separated basenames without .csv)
GENE_GROUPS="${GENE_GROUPS:-SmelDMPs_v5,SmelGRF-GIF_with_Control}"

# Gene groups directory (reference-specific)
GENE_GROUPS_DIR="${GENE_GROUPS_DIR:-${BASE_DIR}/inputs/gene_groups_csv/experimental/GPE001970}"

# SRR CSV directory for sample labels
SRR_CSV_DIR="${SRR_CSV_DIR:-${BASE_DIR}/inputs/SRR_csv}"

# System resources
THREADS="${THREADS:-64}"
ENABLE_GPU="${ENABLE_GPU:-FALSE}"
AVAILABLE_RAM_GB="${AVAILABLE_RAM_GB:-128}"
GPU_VRAM_GB="${GPU_VRAM_GB:-8}"

#===============================================================================
# PATHS
#===============================================================================

OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/3_POST_PROC/cross_method_concordance}"
ALIGNMENT_BASE="${ALIGNMENT_BASE:-${BASE_DIR}/2_ALIGNMENT_RESULTs}"
POST_PROC_BASE="${POST_PROC_BASE:-${BASE_DIR}/3_POST_PROC}"
ANALYSIS_MODULES_DIR="${BASE_DIR}/modules/c_post_processing/analysis_modules"
CONCORDANCE_SCRIPT_DIR="${BASE_DIR}/modules/c_post_processing/cross_method_concordance"

mkdir -p "${OUTPUT_DIR}/figures" "${OUTPUT_DIR}/tables"

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
export CONCORDANCE_SCRIPT_DIR ANALYSIS_MODULES_DIR
export OUTPUT_DIR ALIGNMENT_BASE POST_PROC_BASE

#===============================================================================
# RUN ANALYSIS PIPELINE
#===============================================================================

echo "============================================================"
echo "  CROSS-METHOD CONCORDANCE ANALYSIS"
echo "============================================================"
echo "  Reference:    ${MASTER_REFERENCE}"
echo "  Methods:      ${METHODS}"
echo "  Gene groups:  ${GENE_GROUPS}"
echo "  Output:       ${OUTPUT_DIR}"
echo "  Threads:      ${THREADS}"
echo "============================================================"
echo ""

run_step() {
    local step_num="$1"
    local step_name="$2"
    local script="$3"

    echo "------------------------------------------------------------"
    echo "  [STEP ${step_num}/4] ${step_name}"
    echo "------------------------------------------------------------"

    if ! Rscript "${CONCORDANCE_SCRIPT_DIR}/${script}"; then
        echo "[ERROR] Step ${step_num} (${step_name}) failed!"
        exit 1
    fi
    echo ""
}

run_step 1 "Load & Harmonize Matrices"    "1_load_matrices.R"
run_step 2 "Quantification Concordance"   "2_quantification_concordance.R"
run_step 3 "Ranking Stability"            "3_ranking_stability.R"
run_step 4 "Generate Report"              "4_generate_report.R"

echo ""
echo "============================================================"
echo "  CONCORDANCE ANALYSIS COMPLETE"
echo "============================================================"
echo "  Report:  ${POST_PROC_BASE}/cross_method_concordance_report.md"
echo "  Figures: ${OUTPUT_DIR}/figures/"
echo "  Tables:  ${OUTPUT_DIR}/tables/"
echo "============================================================"
