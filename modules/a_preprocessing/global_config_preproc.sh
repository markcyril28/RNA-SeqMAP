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
export PREPROC_CONFIG_SOURCED="true"

# ==============================================================================
# IMPORTANT PARAMETERS - MODIFY THESE
# ==============================================================================

# Thread and Job Configuration
THREADS="${THREADS:-12}"
JOBS="${JOBS:-2}"
THREADS_PER_JOB="${THREADS_PER_JOB:-$((THREADS / JOBS))}"
[[ $THREADS_PER_JOB -lt 1 ]] && THREADS_PER_JOB=1

# GNU Parallel Configuration
USE_GNU_PARALLEL="${USE_GNU_PARALLEL:-FALSE}"

# Cleanup Configuration
DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING="${DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING:-FALSE}"

# ==============================================================================
# DIRECTORY STRUCTURE
# ==============================================================================
RAW_DIR_ROOT="${RAW_DIR_ROOT:-1_RAW_SRR}"
TRIM_DIR_ROOT="${TRIM_DIR_ROOT:-2_TRIMMED_SRR}"
FASTQC_ROOT="${FASTQC_ROOT:-3_FastQC}"

# ==============================================================================
# TRIMMING PARAMETER PROFILES
# ==============================================================================
# Each profile contains: HEADCROP_BASES, TAILCROP_BASES, MINLEN, SLIDINGWINDOW
# Format: "HEADCROP:TAILCROP:MINLEN:SLIDINGWINDOW_SIZE:SLIDINGWINDOW_QUAL"

# Default trimming profile (fallback)
TRIM_PROFILE_DEFAULT="12:0:36:4:20"

# Profile for PRJNA328564 (Main Dataset - older sequencing)
TRIM_PROFILE_PRJNA328564="12:0:36:4:20"

# Profile for SAMN28540077 (Chinese Dataset - high quality)
TRIM_PROFILE_SAMN28540077="12:0:36:4:20"

# Profile for SAMN28540068 (Chinese Dataset - high quality)
TRIM_PROFILE_SAMN28540068="12:0:36:4:20"

# Profile for PRJNA865018 (SmelDMP GEA Set 1)
TRIM_PROFILE_PRJNA865018="12:0:36:4:20"

# Profile for PRJNA941250 (SmelDMP GEA Set 2)
TRIM_PROFILE_PRJNA941250="12:0:36:4:20"

# Profile for OTHER_SRR_LIST
TRIM_PROFILE_OTHER="12:0:36:4:20"

# Declare associative array to map SRR IDs to their trim profiles
declare -gA SRR_TRIM_PROFILE_MAP

# ==============================================================================
# SRR TO TRIM PROFILE MAPPING
# ==============================================================================

init_srr_trim_profiles() {
	# Map PRJNA328564 SRRs
	local prjna328564_srrs=(
		SRR3884685 SRR3884677 SRR3884675 SRR3884690 SRR3884689 SRR3884684
		SRR3884686 SRR3884687 SRR3884597 SRR3884679 SRR3884608 SRR3884620
		SRR3884631 SRR3884642 SRR3884653 SRR3884664 SRR3884680 SRR3884681 SRR3884678
	)
	for srr in "${prjna328564_srrs[@]}"; do
		SRR_TRIM_PROFILE_MAP["$srr"]="$TRIM_PROFILE_PRJNA328564"
	done

	# Map SAMN28540077 SRRs
	local samn28540077_srrs=(
		SRR20722232 SRR20722226 SRR20722234 SRR20722228 SRR4243802
		SRR20722233 SRR20722230 SRR20722227 SRR20722229
	)
	for srr in "${samn28540077_srrs[@]}"; do
		SRR_TRIM_PROFILE_MAP["$srr"]="$TRIM_PROFILE_SAMN28540077"
	done

	# Map SAMN28540068 SRRs
	local samn28540068_srrs=(
		SRR20722387 SRR20722297 SRR20722385 SRR20722296 SRR20722386
		SRR20722383 SRR20722384 SRR31755282
	)
	for srr in "${samn28540068_srrs[@]}"; do
		SRR_TRIM_PROFILE_MAP["$srr"]="$TRIM_PROFILE_SAMN28540068"
	done

	# Map PRJNA865018 SRRs
	local prjna865018_srrs=(
		SRR21010466 SRR21010456 SRR21010454 SRR21010462 SRR21010460
		SRR21010458 SRR21010452 SRR21010450 SRR21010464
	)
	for srr in "${prjna865018_srrs[@]}"; do
		SRR_TRIM_PROFILE_MAP["$srr"]="$TRIM_PROFILE_PRJNA865018"
	done

	# Map PRJNA941250 SRRs
	local prjna941250_srrs=(
		SRR23909869 SRR23909870 SRR23909871 SRR23909866 SRR23909867
		SRR23909868 SRR23909863 SRR23909864 SRR23909865
	)
	for srr in "${prjna941250_srrs[@]}"; do
		SRR_TRIM_PROFILE_MAP["$srr"]="$TRIM_PROFILE_PRJNA941250"
	done

	# Map OTHER SRRs
	local other_srrs=(SRR34564302 SRR34848077 SRR3479277)
	for srr in "${other_srrs[@]}"; do
		SRR_TRIM_PROFILE_MAP["$srr"]="$TRIM_PROFILE_OTHER"
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

# Function to display trim profile for an SRR
show_trim_profile() {
	local srr="$1"
	local profile="${SRR_TRIM_PROFILE_MAP[$srr]:-$TRIM_PROFILE_DEFAULT}"
	echo "Trim profile for $srr: $profile"
}

# Initialize the mapping on source
init_srr_trim_profiles

# ==============================================================================
# INITIALIZE PREPROCESSING DIRECTORIES
# ==============================================================================
init_preproc_directories() {
	mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT"
}
