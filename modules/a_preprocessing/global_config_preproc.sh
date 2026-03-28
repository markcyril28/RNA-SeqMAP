#!/bin/bash
# ==============================================================================
# PREPROCESSING GLOBAL CONFIGURATION
# ==============================================================================
# Centralized configuration for all preprocessing operations
# Sourced by: download.sh, trimming.sh, quality_checks.sh
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${PREPROC_CONFIG_SOURCED:-}" == "true" ]] && return 0
PREPROC_CONFIG_SOURCED="true"

# ==============================================================================
# IMPORTANT PARAMETERS (tweak here)
# ==============================================================================

# Total CPU threads available to the pipeline (auto-detect if not set)
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-${PBS_NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 12)}}}"
# Number of parallel jobs (GNU Parallel); THREADS_PER_JOB is auto-calculated
# Set to "auto" to calculate from THREADS / OPTIMAL_THREADS_PER_JOB
JOBS="${JOBS:-2}"
# Optimal threads per job for Stage 1 (TrimGalore = 4 threads optimal)
OPTIMAL_THREADS_PER_JOB="${OPTIMAL_THREADS_PER_JOB:-4}"
# Resolve JOBS="auto" → numeric (uses _resolve_auto_jobs from runtime_defaults.sh
# if already sourced, otherwise inline resolution for standalone use)
if [[ "${JOBS}" == "auto" || "${JOBS}" == "AUTO" ]]; then
	(( OPTIMAL_THREADS_PER_JOB < 1 )) && OPTIMAL_THREADS_PER_JOB=4
	JOBS=$(( THREADS / OPTIMAL_THREADS_PER_JOB ))
	(( JOBS < 1 )) && JOBS=1
fi
# Guard: ensure JOBS >= 1 before division (prevents bash arithmetic error if externally set to 0)
(( JOBS < 1 )) && JOBS=1
THREADS_PER_JOB="${THREADS_PER_JOB:-$((THREADS / JOBS))}"
[[ $THREADS_PER_JOB -lt 1 ]] && THREADS_PER_JOB=1

# Enable GNU Parallel for parallel sample processing (TRUE/FALSE)
USE_GNU_PARALLEL="${USE_GNU_PARALLEL:-FALSE}"

# Delete raw SRR files after successful download + trimming (TRUE/FALSE)
DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING="${DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING:-FALSE}"

# ==============================================================================
# DIRECTORY STRUCTURE
# ==============================================================================
# Orchestrated execution (Nextflow/Snakemake) must set PROJECT_ROOT explicitly;
# the "." fallback only works when CWD is the repo root (standalone bash mode).
if [[ -n "${WF_MANAGED_ENV:-}" && -z "${PROJECT_ROOT:-}" ]]; then
	echo "ERROR: WF_MANAGED_ENV is set but PROJECT_ROOT is not. Orchestrators must export PROJECT_ROOT." >&2
	return 1
fi
RAW_DIR_ROOT="${RAW_DIR_ROOT:-${PROJECT_ROOT:-.}/1_SRRs/A_RAW_SRR}"
TRIM_DIR_ROOT="${TRIM_DIR_ROOT:-${PROJECT_ROOT:-.}/1_SRRs/B_TRIMMED_SRR}"
FASTQC_ROOT="${FASTQC_ROOT:-${PROJECT_ROOT:-.}/1_SRRs/C_FastQC}"

# ==============================================================================
# TRIMMING PARAMETER PROFILES
# ==============================================================================
# Each profile contains: HEADCROP_BASES, TAILCROP_BASES, MINLEN, SLIDINGWINDOW
# Format: "HEADCROP:TAILCROP:MINLEN:SLIDINGWINDOW_SIZE:SLIDINGWINDOW_QUAL"

# Default trimming profile (used for all datasets — currently identical parameters)
# Format: HEADCROP:TAILCROP:MINLEN:SLIDINGWINDOW_SIZE:SLIDINGWINDOW_QUAL
# To add dataset-specific profiles, define TRIM_PROFILE_<DATASET> and update init_srr_trim_profiles()
TRIM_PROFILE_DEFAULT="12:0:36:4:20"

# Declare associative array to map SRR IDs to their trim profiles
declare -gA SRR_TRIM_PROFILE_MAP 2>/dev/null || declare -A SRR_TRIM_PROFILE_MAP

# ==============================================================================
# SRR TO TRIM PROFILE MAPPING
# ==============================================================================

# O(S) — pre-populate all known SRR IDs with TRIM_PROFILE_DEFAULT.
# When dataset-specific profiles are needed, override individual entries after this loop.
init_srr_trim_profiles() {
	local _all_known_srrs=(
		# PRJNA328564 (Main Dataset)
		SRR3884685 SRR3884677 SRR3884675 SRR3884690 SRR3884689 SRR3884684
		SRR3884686 SRR3884687 SRR3884597 SRR3884679 SRR3884608 SRR3884620
		SRR3884631 SRR3884642 SRR3884653 SRR3884664 SRR3884680 SRR3884681 SRR3884678
		# SAMN28540077 (Chinese Dataset)
		SRR20722232 SRR20722226 SRR20722234 SRR20722228 SRR4243802
		SRR20722233 SRR20722230 SRR20722227 SRR20722229
		# SAMN28540068 (Chinese Dataset)
		SRR20722387 SRR20722297 SRR20722385 SRR20722296 SRR20722386
		SRR20722383 SRR20722384 SRR31755282
		# PRJNA865018 (SmelDMP GEA Set 1)
		SRR21010466 SRR21010456 SRR21010454 SRR21010462 SRR21010460
		SRR21010458 SRR21010452 SRR21010450 SRR21010464
		# PRJNA941250 (SmelDMP GEA Set 2)
		SRR23909869 SRR23909870 SRR23909871 SRR23909866 SRR23909867
		SRR23909868 SRR23909863 SRR23909864 SRR23909865
		# Other
		SRR34564302 SRR34848077 SRR3479277
	)
	local srr
	for srr in "${_all_known_srrs[@]}"; do
		SRR_TRIM_PROFILE_MAP["$srr"]="$TRIM_PROFILE_DEFAULT"
	done
}

# Function to get trim parameters for a specific SRR
# Usage: get_trim_params SRR_ID
# Returns: Sets global variables HEADCROP_BASES, TAILCROP_BASES, MINLEN, SW_SIZE, SW_QUAL
get_trim_params() {
	local srr="$1"
	local profile="${SRR_TRIM_PROFILE_MAP[$srr]:-$TRIM_PROFILE_DEFAULT}"
	
	IFS=':' read -r HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL <<< "$profile"
	export HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL
}

# Initialize the mapping on source
init_srr_trim_profiles
