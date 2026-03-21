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
#   0_input_information/ - Sample metadata files
# ==============================================================================

# Guard against double-sourcing
[[ "${MODULES_LOADER_SOURCED:-}" == "true" ]] && return 0
export MODULES_LOADER_SOURCED="true"

# Export MODULES_DIR so child modules can derive SCRIPT_DIR without cd+dirname+pwd subshells
MODULES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MODULES_DIR

# ==============================================================================
# LOAD MODULES IN DEPENDENCY ORDER
# ==============================================================================

# 1. Logging and GPU utilities (no dependencies)
source "$MODULES_DIR/logging/logging_utils.sh"
source "$MODULES_DIR/logging/gpu_utils.sh"

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
# INITIALIZE DIRECTORIES
# ==============================================================================

init_directories() {
	mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT"
	init_method_directories
	# Per-stage log directories (1_SRRs, 2_ALIGNMENT_RESULTs, 3_POST_PROC)
	local stage_dir
	for stage_dir in "1_SRRs" "2_ALIGNMENT_RESULTs" "3_POST_PROC"; do
		mkdir -p "${stage_dir}/logs/log_files" "${stage_dir}/logs/time_logs" "${stage_dir}/logs/space_logs" \
			"${stage_dir}/logs/space_time_logs" "${stage_dir}/logs/error_warn_logs" "${stage_dir}/logs/software_catalogs" \
			"${stage_dir}/logs/gpu_log"
	done
}

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
