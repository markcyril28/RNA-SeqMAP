#!/bin/bash
# ==============================================================================
# DOWNLOAD FUNCTIONS
# ==============================================================================
# SRA data download utilities with multiple download methods
# ==============================================================================
 
#set -euo pipefail

# Guard against double-sourcing
[[ "${DOWNLOAD_SOURCED:-}" == "true" ]] && return 0
export DOWNLOAD_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/a_preprocessing}"
if [[ -z "$SCRIPT_DIR" ]]; then SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."; fi
source "$SCRIPT_DIR/shared_utils_preproc.sh"

# ==============================================================================
# SKIP-IF-EXISTS CHECK (DRY helper — consolidates 4 identical check patterns)
# ==============================================================================
# Returns 0 (skip) if trimmed or raw files already exist for this SRR.
# O(1) glob checks; avoids duplicating the find_trimmed/find_raw pattern.
_srr_already_downloaded() {
	local SRR="$1"
	find_trimmed_fastq "$SRR"
	[[ -n "$trimmed1" ]] && return 0
	find_raw_fastq "$SRR"
	[[ -n "$raw1" ]] && return 0
	return 1
}

# ==============================================================================
# PRIMARY: SRA Toolkit (prefetch + fasterq-dump)
# ==============================================================================

download_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided for download"; return 1; }
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir"
		
		_srr_already_downloaded "$SRR" && { log_info "Files for $SRR exist. Skipping download."; continue; }
		
		log_info "Downloading $SRR..."
		run_with_space_time_log prefetch "$SRR" --output-directory "$raw_dir"
		run_with_space_time_log fasterq-dump --split-files --threads "$THREADS" \
			"$raw_dir/$SRR/$SRR.sra" -O "$raw_dir"
		
		# Compress downloaded files (use shared pigz detection, avoid per-SRR command -v spawn)
		local -a _ccmd=("${_SHARED_GZIP_C:-gzip}")
		[[ "${_ccmd[0]}" == "pigz" ]] && _ccmd=(pigz -p "${THREADS:-4}")
		local _p1 _p2
		if [[ -f "$raw_dir/${SRR}_1.fastq" && -f "$raw_dir/${SRR}_2.fastq" ]]; then
			"${_ccmd[@]}" "$raw_dir/${SRR}_1.fastq" & _p1=$!
			"${_ccmd[@]}" "$raw_dir/${SRR}_2.fastq" & _p2=$!
			wait "$_p1" "$_p2"
		else
			[[ -f "$raw_dir/${SRR}_1.fastq" ]] && "${_ccmd[@]}" "$raw_dir/${SRR}_1.fastq"
			[[ -f "$raw_dir/${SRR}_2.fastq" ]] && "${_ccmd[@]}" "$raw_dir/${SRR}_2.fastq"
		fi
	done
	log_info "All downloads completed."
}


# ==============================================================================
# PARALLEL DOWNLOAD (GNU Parallel wrapper for primary SRA download)
# ==============================================================================

download_srrs_parallel() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }
	
	if ! should_use_parallel; then
		log_info "Running downloads sequentially (USE_GNU_PARALLEL=${USE_GNU_PARALLEL:-FALSE})"
		download_srrs "${SRR_LIST[@]}"
		return $?
	fi
	
	log_info "Running parallel downloads with GNU Parallel (JOBS=${JOBS:-2})"
	export PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE
	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS THREADS_PER_JOB
	# Export _log_impl (core logger) alongside its callers — without it, log_info/log_warn/log_error
	# fail silently in GNU Parallel subshells because they delegate to _log_impl.
	export -f _log_impl timestamp log log_info log_warn log_error find_trimmed_fastq find_raw_fastq _srr_already_downloaded
	
	_download_worker() {
		local SRR="$1"
		
		# Activate conda environment in subshell (uses cached path to avoid dirname subshell)
		if [[ -n "$CONDA_PREFIX" ]]; then
			source "${_CONDA_PROFILE_SCRIPT:-${CONDA_EXE%/*}/../etc/profile.d/conda.sh}" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir"
		
		_srr_already_downloaded "$SRR" && return 0
		
		prefetch "$SRR" --output-directory "$raw_dir" || return 1
		fasterq-dump --split-files --threads "${THREADS_PER_JOB:-2}" "$raw_dir/$SRR/$SRR.sra" -O "$raw_dir" || return 1
		local -a _ccmd=("${_SHARED_GZIP_C:-gzip}")
		[[ "${_ccmd[0]}" == "pigz" ]] && _ccmd=(pigz -p "${THREADS_PER_JOB:-2}")
		# Compress R1/R2 concurrently — paired-end gzip is I/O-bound, ~2x faster with overlap
		[[ -f "$raw_dir/${SRR}_1.fastq" ]] && "${_ccmd[@]}" "$raw_dir/${SRR}_1.fastq" &
		[[ -f "$raw_dir/${SRR}_2.fastq" ]] && "${_ccmd[@]}" "$raw_dir/${SRR}_2.fastq" &
		wait
	}
	export -f _download_worker
	
	printf "%s\n" "${SRR_LIST[@]}" | parallel \
		--env PATH --env CONDA_PREFIX --env _SHARED_GZIP_C \
		-j "${JOBS:-2}" \
		--halt soon,fail,1 \
		--joblog "$RAW_DIR_ROOT/parallel_download.log" \
		_download_worker {}
}
