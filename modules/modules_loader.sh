#!/bin/bash
# ==============================================================================
# MODULES LOADER - SINGLE ENTRY POINT
# ==============================================================================
# Loads all HeatSeq pipeline modules in correct dependency order
# Usage: source "modules/modules_loader.sh"
#
# Structure:
#   logging/           - Logging utilities
#   a_preprocessing/   - Download, trimming, QC functions
#   b_main_methods/    - GEA analysis methods (M1-M5)
#   c_post_processing/   - Post-processing analysis, concordance, utilities
# ==============================================================================

# Guard against double-sourcing
[[ "${MODULES_LOADER_SOURCED:-}" == "true" ]] && return 0
export MODULES_LOADER_SOURCED="true"

# Export MODULES_DIR so child modules can derive SCRIPT_DIR without cd+dirname+pwd subshells
# Resolve without nested dirname subshell — one cd instead of two forks
MODULES_DIR="${BASH_SOURCE[0]%/*}"
[[ "$MODULES_DIR" == "${BASH_SOURCE[0]}" ]] && MODULES_DIR="."
MODULES_DIR="$(cd "$MODULES_DIR" && pwd)"
export MODULES_DIR

# ==============================================================================
# LOAD MODULES IN DEPENDENCY ORDER
# ==============================================================================

# 1. Logging utilities (no dependencies)
source "$MODULES_DIR/logging/logging_utils.sh"
# NOTE: gpu_utils.sh is NOT sourced here — R does its own GPU detection in
# 0_shared_config.R, and no active pipeline code calls any gpu_utils function.
# Source it explicitly if needed: source "$MODULES_DIR/logging/gpu_utils.sh"

# 2. Preprocessing modules
source "$MODULES_DIR/a_preprocessing/global_config_preproc.sh"
source "$MODULES_DIR/a_preprocessing/shared_utils_preproc.sh"
source "$MODULES_DIR/a_preprocessing/download.sh"
source "$MODULES_DIR/a_preprocessing/trimming.sh"
source "$MODULES_DIR/a_preprocessing/quality_checks.sh"

# 3. Main methods modules
source "$MODULES_DIR/b_main_methods/global_config_method.sh"
source "$MODULES_DIR/b_main_methods/shared_utils_method.sh"
source "$MODULES_DIR/b_main_methods/methods_loader.sh"


# ==============================================================================
# AVAILABLE FUNCTIONS (for reference)
# ==============================================================================
# Preprocessing: download_srrs, download_srrs_parallel, trim_srrs,
#                trim_srrs_trimmomatic, download_and_trim_srrs,
#                run_quality_control, run_quality_control_all
#
# Methods: hisat2_ref_guided_pipeline (M1), hisat2_de_novo_pipeline (M2),
#          star_alignment_pipeline (M3), salmon_saf_pipeline (M4),
#          bowtie2_rsem_pipeline (M5)
# ==============================================================================
