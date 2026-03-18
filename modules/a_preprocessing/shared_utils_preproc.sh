#!/bin/bash
# ==============================================================================
# PREPROCESSING SHARED UTILITIES
# ==============================================================================
# Common helper functions used across preprocessing scripts
# Sourced by: download.sh, trimming.sh, quality_checks.sh
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${PREPROC_SHARED_SOURCED:-}" == "true" ]] && return 0
export PREPROC_SHARED_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/global_config_preproc.sh"
source "$SCRIPT_DIR/../logging/logging_utils.sh"

# ==============================================================================
# SHARED COMPRESSION DETECTION (preprocessing)
# ==============================================================================
# Detect pigz once at module load so download.sh/trimming.sh can use it without
# spawning `command -v pigz` per sample.  shared_utils_method.sh has its own
# detection — both export the same variable names so whichever loads first wins.
if [[ -z "${_SHARED_GZIP_C:-}" ]]; then
	if command -v pigz &>/dev/null; then
		_SHARED_HAS_PIGZ="true"
		_SHARED_GZIP_DC="pigz -dc"
		_SHARED_GZIP_C="pigz"
	else
		_SHARED_HAS_PIGZ="false"
		_SHARED_GZIP_DC="gzip -dc"
		_SHARED_GZIP_C="gzip"
	fi
	export _SHARED_HAS_PIGZ _SHARED_GZIP_DC _SHARED_GZIP_C
fi

# ==============================================================================
# FASTQ FILE DETECTION FUNCTIONS
# ==============================================================================

# Find trimmed FASTQ files for a given SRR ID
# Sets: trimmed1, trimmed2 (empty if single-end)
find_trimmed_fastq() {
	local SRR="$1"
	local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
	trimmed1="" trimmed2=""
	
	# Paired-end patterns (compressed first — TrimGalore default output is .fq.gz)
	if [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq.gz" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq.gz" ]]; then
		trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq.gz"
		trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq.gz"
	elif [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq" ]]; then
		trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
		trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
	else
		# Glob fallback for non-standard paired-end names (no subshell via compgen)
		local files1=("$TrimGalore_DIR"/${SRR}*val_1.*)
		if [[ -f "${files1[0]:-}" ]]; then
			local files2=("$TrimGalore_DIR"/${SRR}*val_2.*)
			trimmed1="${files1[0]}"
			[[ -f "${files2[0]:-}" ]] && trimmed2="${files2[0]}"
		# Single-end patterns (compressed first)
		elif [[ -f "$TrimGalore_DIR/${SRR}_trimmed.fq.gz" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_trimmed.fq.gz"
		elif [[ -f "$TrimGalore_DIR/${SRR}_trimmed.fq" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_trimmed.fq"
		else
			local files=("$TrimGalore_DIR"/${SRR}*trimmed.fq*)
			[[ -f "${files[0]:-}" ]] && trimmed1="${files[0]}"
		fi
	fi
}

# Find raw FASTQ files for a given SRR ID
# Sets: raw1, raw2 (empty if single-end or not found)
find_raw_fastq() {
	local SRR="$1"
	local raw_dir="$RAW_DIR_ROOT/$SRR"
	raw1="" raw2=""
	
	if [[ -f "$raw_dir/${SRR}_1.fastq" && -f "$raw_dir/${SRR}_2.fastq" ]]; then
		raw1="$raw_dir/${SRR}_1.fastq"
		raw2="$raw_dir/${SRR}_2.fastq"
	elif [[ -f "$raw_dir/${SRR}_1.fastq.gz" && -f "$raw_dir/${SRR}_2.fastq.gz" ]]; then
		raw1="$raw_dir/${SRR}_1.fastq.gz"
		raw2="$raw_dir/${SRR}_2.fastq.gz"
	# Single-end patterns
	elif [[ -f "$raw_dir/${SRR}.fastq" ]]; then
		raw1="$raw_dir/${SRR}.fastq"
	elif [[ -f "$raw_dir/${SRR}.fastq.gz" ]]; then
		raw1="$raw_dir/${SRR}.fastq.gz"
	fi
}

# ==============================================================================
# TRIMMING VERIFICATION
# ==============================================================================

# Verify trimming success and optionally cleanup raw files
verify_trimming_and_cleanup() {
	local SRR="$1"
	local trimmed1="$2"
	local trimmed2="$3"
	local raw1="${4:-}"
	local raw2="${5:-}"
	
	local success=false
	if { [[ -f "$trimmed1" && -s "$trimmed1" ]] || [[ -f "${trimmed1}.gz" && -s "${trimmed1}.gz" ]]; } && \
	   { [[ -z "$trimmed2" ]] || [[ -f "$trimmed2" && -s "$trimmed2" ]] || [[ -f "${trimmed2}.gz" && -s "${trimmed2}.gz" ]]; }; then
		success=true
		log_info "Trimming completed for $SRR"
		
		if [[ "$DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING" == "TRUE" && -n "$raw1" ]]; then
			log_info "Cleaning up raw files for $SRR..."
			rm -f "$raw1" "$raw2"
			rm -rf "$RAW_DIR_ROOT/$SRR/$SRR"
			rmdir "$RAW_DIR_ROOT/$SRR" 2>/dev/null || true
		fi
	else
		log_warn "Trimming may have failed for $SRR - keeping raw files"
	fi
	
	[[ "$success" == "true" ]]
}

# ==============================================================================
# READ LENGTH DETECTION
# ==============================================================================

# Detect read length from FASTQ file
detect_read_length() {
	local fastq="$1"
	local default_length="${2:-150}"
	
	[[ ! -f "$fastq" ]] && { echo "$default_length"; return 1; }
	
	# Use cached pigz detection (set at module load) instead of per-call command -v
	local decompress_cmd="cat"
	case "$fastq" in
		*.gz) decompress_cmd="${_SHARED_GZIP_DC:-zcat}" ;;
		*.bz2) decompress_cmd="bzcat" ;;
	esac
	
	local avg_length=$($decompress_cmd "$fastq" 2>/dev/null | \
		awk 'NR%4==2 {sum+=length($0); count++} count==1000 {print int(sum/count); exit} END {if (count>0 && count<1000) print int(sum/count)}')
	
	if [[ -z "$avg_length" || $avg_length -lt 50 || $avg_length -gt 300 ]]; then
		echo "$default_length"
		return 1
	fi
	echo "$avg_length"
}

# ==============================================================================
# PARALLEL PROCESSING UTILITIES
# ==============================================================================

# Check if GNU Parallel should be used
# Returns 0 (true) if parallel should be used, 1 (false) otherwise
should_use_parallel() {
	if [[ "${USE_GNU_PARALLEL:-FALSE}" != "TRUE" ]]; then
		return 1
	fi
	if ! command -v parallel &>/dev/null; then
		log_warn "GNU Parallel requested but not installed. Falling back to sequential processing."
		return 1
	fi
	if [[ "${JOBS:-1}" -le 1 ]]; then
		return 1
	fi
	return 0
}

# ==============================================================================
# COMPRESSION UTILITIES
# ==============================================================================

gzip_trimmed_fastq_files() {
	log_info "Compressing trimmed FASTQ files in $TRIM_DIR_ROOT..."
	local _compress_cmd="gzip" _parallel_jobs
	# Use cached pigz detection (set at module load) instead of per-call command -v
	if [[ "${_SHARED_HAS_PIGZ:-false}" == "true" ]]; then
		_compress_cmd="pigz -p ${THREADS_PER_JOB:-4}"
		_parallel_jobs="${JOBS:-2}"
		log_info "Using pigz for multi-threaded compression (${_parallel_jobs} jobs x ${THREADS_PER_JOB:-4} threads)"
	else
		# gzip is single-threaded: use all available threads as parallel jobs
		_parallel_jobs="${THREADS:-$(nproc 2>/dev/null || echo 4)}"
		log_info "Using gzip with ${_parallel_jobs} parallel jobs"
	fi
	local _compress_rc=0
	# No -I {}: lets xargs batch multiple files per invocation (fewer process spawns)
	find "$TRIM_DIR_ROOT" -type f -name "*.fq" -print0 | \
		xargs -0 -P "$_parallel_jobs" $_compress_cmd 2>/dev/null || _compress_rc=$?
	if [[ $_compress_rc -ne 0 ]]; then
		log_warn "Compression finished with errors (exit code: $_compress_rc) — some .fq files may not have been compressed."
	else
		log_info "Compression completed successfully."
	fi
}

# ==============================================================================
# CLEANUP UTILITIES
# ==============================================================================

# Delete trimmed FASTQ files for a list of SRR IDs
# Usage: delete_trimmed_fastq_by_srr_list SRR1 SRR2 SRR3 ...
delete_trimmed_fastq_by_srr_list() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_warn "No SRR IDs provided for deletion"; return 1; }
	
	log_info "Deleting trimmed FASTQ files for ${#SRR_LIST[@]} SRR(s)..."
	local deleted_count=0
	
	for SRR in "${SRR_LIST[@]}"; do
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		if [[ -d "$trim_dir" ]]; then
			log_info "Deleting trimmed files for $SRR..."
			rm -rf "$trim_dir"
			((deleted_count++)) || true
		else
			log_warn "Trimmed directory not found for $SRR: $trim_dir"
		fi
	done
	
	log_info "Deleted trimmed files for $deleted_count SRR(s)."
}

# Delete raw FASTQ files for a list of SRR IDs
# Usage: delete_raw_srr_by_srr_list SRR1 SRR2 SRR3 ...
delete_raw_srr_by_srr_list() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_warn "No SRR IDs provided for deletion"; return 1; }
	
	log_info "Deleting raw SRR files for ${#SRR_LIST[@]} SRR(s)..."
	local deleted_count=0
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		if [[ -d "$raw_dir" ]]; then
			log_info "Deleting raw files for $SRR..."
			rm -rf "$raw_dir"
			((deleted_count++)) || true
		else
			log_warn "Raw directory not found for $SRR: $raw_dir"
		fi
	done
	
	log_info "Deleted raw files for $deleted_count SRR(s)."
}
