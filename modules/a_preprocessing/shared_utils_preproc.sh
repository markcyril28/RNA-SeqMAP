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
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/a_preprocessing}"
if [[ -z "$SCRIPT_DIR" ]]; then SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."; fi
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

# Cache conda profile path at module load — avoids dirname subshell per parallel worker.
if [[ -z "${_CONDA_PROFILE_SCRIPT+x}" && -n "${CONDA_EXE:-}" ]]; then
	_CONDA_PROFILE_SCRIPT="${CONDA_EXE%/*}/../etc/profile.d/conda.sh"
	export _CONDA_PROFILE_SCRIPT
fi

# Cache CPU count at module load — avoids nproc subprocess spawn per function call.
# Used by gzip_trimmed_fastq_files() and other utilities as thread count fallback.
if [[ -z "${_CACHED_NPROC:-}" ]]; then
	_CACHED_NPROC=$(nproc 2>/dev/null || echo 4)
	export _CACHED_NPROC
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
		# Single glob expansion for both val_1 and val_2 — O(1) readdir vs O(2).
		local _all_vals=("$TrimGalore_DIR"/${SRR}*val_[12].*)
		if [[ -f "${_all_vals[0]:-}" ]]; then
			local _v
			for _v in "${_all_vals[@]}"; do
				[[ "$_v" == *val_1* && -f "$_v" ]] && trimmed1="$_v"
				[[ "$_v" == *val_2* && -f "$_v" ]] && trimmed2="$_v"
				# Early exit once both files found — avoids O(n) unnecessary glob iterations
				[[ -n "$trimmed1" && -n "$trimmed2" ]] && break
			done
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
			if [[ -z "${RAW_DIR_ROOT:-}" || -z "$SRR" ]]; then
				log_error "RAW_DIR_ROOT or SRR is empty — refusing to delete"
				return 1
			fi
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

# Global read length cache — avoids redundant FASTQ decompression when multiple
# methods (M1, M2, M3, M5) detect read length from the same trimmed files.
# Key: file path, Value: detected read length. Persists across method calls.
declare -gA _READ_LENGTH_CACHE 2>/dev/null || declare -A _READ_LENGTH_CACHE

# Detect read length from FASTQ file
detect_read_length() {
	local fastq="$1"
	local default_length="${2:-150}"

	[[ ! -f "$fastq" ]] && { echo "$default_length"; return 1; }

	# O(1) cache check — skip decompression if already detected for this file
	if [[ -n "${_READ_LENGTH_CACHE[$fastq]:-}" ]]; then
		echo "${_READ_LENGTH_CACHE[$fastq]}"
		return 0
	fi

	# Use cached pigz detection (set at module load) instead of per-call command -v
	local decompress_cmd="cat"
	case "$fastq" in
		*.gz) decompress_cmd="${_SHARED_GZIP_DC:-zcat}" ;;
		*.bz2) decompress_cmd="bzcat" ;;
	esac

	# Sample 100 reads (400 lines) instead of 1000 — statistically equivalent for
	# read length detection but 10x fewer lines decompressed. O(400) vs O(4000).
	# head merged into awk (NR>400{exit}) — eliminates 1 subprocess per cache miss.
	local avg_length=$($decompress_cmd "$fastq" 2>/dev/null | \
		awk 'NR>400{exit} NR%4==2 {sum+=length($0); count++} END {if (count>0) print int(sum/count)}')

	if [[ -z "$avg_length" || $avg_length -lt 50 || $avg_length -gt 300 ]]; then
		echo "$default_length"
		return 1
	fi

	# Cache result for subsequent calls with the same file
	_READ_LENGTH_CACHE["$fastq"]="$avg_length"
	echo "$avg_length"
}

# ==============================================================================
# PARALLEL PROCESSING UTILITIES
# ==============================================================================

# Check if GNU Parallel should be used
# Returns 0 (true) if parallel should be used, 1 (false) otherwise
# O(1) check using cached binary availability from shared_utils_method.sh or local cache.
# Avoids per-call `command -v parallel` subprocess spawn.
should_use_parallel() {
	[[ "${USE_GNU_PARALLEL:-FALSE}" != "TRUE" ]] && return 1
	# Use cached detection if available, else cache now (first call)
	if [[ -z "${_SHARED_HAS_PARALLEL:-}" ]]; then
		_SHARED_HAS_PARALLEL=false
		command -v parallel &>/dev/null && _SHARED_HAS_PARALLEL=true
	fi
	if [[ "$_SHARED_HAS_PARALLEL" != "true" ]]; then
		log_warn "GNU Parallel requested but not installed. Falling back to sequential processing."
		return 1
	fi
	[[ "${JOBS:-1}" -le 1 ]] && return 1
	return 0
}

# ==============================================================================
# COMPRESSION UTILITIES
# ==============================================================================

gzip_trimmed_fastq_files() {
	# Single find pass: early exit if no .fq files exist (xargs -r / --no-run-if-empty).
	# Eliminates previous double-find pattern (one for check, one for compression).
	# O(tree) single traversal vs O(2 × tree).
	log_info "Compressing trimmed FASTQ files in $TRIM_DIR_ROOT..."
	local _compress_cmd="gzip" _parallel_jobs
	# Use cached pigz detection (set at module load) instead of per-call command -v
	if [[ "${_SHARED_HAS_PIGZ:-false}" == "true" ]]; then
		_compress_cmd="pigz -p ${THREADS_PER_JOB:-4}"
		_parallel_jobs="${JOBS:-2}"
		log_info "Using pigz for multi-threaded compression (${_parallel_jobs} jobs x ${THREADS_PER_JOB:-4} threads)"
	else
		# gzip is single-threaded: use all available threads as parallel jobs
		_parallel_jobs="${THREADS:-${_CACHED_NPROC:-4}}"
		log_info "Using gzip with ${_parallel_jobs} parallel jobs"
	fi
	local _compress_rc=0
	# -r (--no-run-if-empty): xargs exits 0 without spawning compress if find yields nothing.
	# No -I {}: lets xargs batch multiple files per invocation (fewer process spawns)
	find "$TRIM_DIR_ROOT" -type f -name "*.fq" -print0 | \
		xargs -0 -r -P "$_parallel_jobs" $_compress_cmd 2>/dev/null || _compress_rc=$?
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
	[[ -z "${RAW_DIR_ROOT:-}" ]] && { log_error "RAW_DIR_ROOT is empty — refusing to delete"; return 1; }

	log_info "Deleting raw SRR files for ${#SRR_LIST[@]} SRR(s)..."
	local deleted_count=0

	for SRR in "${SRR_LIST[@]}"; do
		[[ -z "$SRR" ]] && { log_warn "Empty SRR ID — skipping"; continue; }
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
