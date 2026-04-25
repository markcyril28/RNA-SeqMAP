#!/bin/bash
# ==============================================================================
# DOWNLOAD FUNCTIONS
# ==============================================================================
# SRA data download utilities with multiple download methods
# ==============================================================================
 
#set -euo pipefail

# Guard against double-sourcing
[[ "${DOWNLOAD_SOURCED:-}" == "true" ]] && return 0
DOWNLOAD_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/a_preprocessing}"
if [[ -z "$SCRIPT_DIR" ]]; then
	SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
	SCRIPT_DIR="$(cd "$SCRIPT_DIR" 2>/dev/null && pwd)"
fi
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
		run_with_space_time_log prefetch "$SRR" --output-directory "$raw_dir" \
			|| { log_error "prefetch failed for $SRR"; continue; }
		run_with_space_time_log fasterq-dump --split-files --threads "$THREADS" \
			"$raw_dir/$SRR/$SRR.sra" -O "$raw_dir" \
			|| { log_error "fasterq-dump failed for $SRR"; continue; }
		
		# Compress downloaded files (use shared pigz detection, avoid per-SRR command -v spawn)
		local -a _ccmd=("${_SHARED_GZIP_C:-gzip}")
		[[ "${_ccmd[0]}" == "pigz" ]] && _ccmd=(pigz -p "${THREADS:-4}")
		local _p1 _p2
		if [[ -f "$raw_dir/${SRR}_1.fastq" && -f "$raw_dir/${SRR}_2.fastq" ]]; then
			"${_ccmd[@]}" "$raw_dir/${SRR}_1.fastq" & _p1=$!
			"${_ccmd[@]}" "$raw_dir/${SRR}_2.fastq" & _p2=$!
			wait "$_p1" "$_p2"
		elif [[ -f "$raw_dir/${SRR}.fastq" ]]; then
			# Single-end: fasterq-dump produces ${SRR}.fastq (no _1/_2 suffix)
			"${_ccmd[@]}" "$raw_dir/${SRR}.fastq"
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
	# Export orchestrator flag and log paths so parallel workers skip conda activation
	# and write errors to the correct log file (prevents contention in Nextflow/Snakemake)
	export WF_MANAGED_ENV LOG_FILE ERROR_WARN_FILE
	# Export _log_impl (core logger) alongside its callers — without it, log_info/log_warn/log_error
	# fail silently in GNU Parallel subshells because they delegate to _log_impl.
	export -f _log_impl timestamp log log_info log_warn log_error find_trimmed_fastq find_raw_fastq _srr_already_downloaded
	
	_download_worker() {
		local SRR="$1"
		
		# Activate conda environment in subshell (skip when orchestrator manages env)
		if [[ -z "${WF_MANAGED_ENV:-}" && -n "${CONDA_PREFIX:-}" ]]; then
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
		# Compress output files — paired-end R1/R2 concurrently, single-end sequentially
		local _cw1=0 _cw2=0
		[[ -f "$raw_dir/${SRR}_1.fastq" ]] && { "${_ccmd[@]}" "$raw_dir/${SRR}_1.fastq" & _cw1=$!; }
		[[ -f "$raw_dir/${SRR}_2.fastq" ]] && { "${_ccmd[@]}" "$raw_dir/${SRR}_2.fastq" & _cw2=$!; }
		[[ "$_cw1" -ne 0 ]] && wait "$_cw1"
		[[ "$_cw2" -ne 0 ]] && wait "$_cw2"
		# Single-end: fasterq-dump produces ${SRR}.fastq (no _1/_2 suffix)
		[[ -f "$raw_dir/${SRR}.fastq" ]] && "${_ccmd[@]}" "$raw_dir/${SRR}.fastq"
	}
	export -f _download_worker

	# Ensure joblog parent dir exists (created by orchestrator in bash mode,
	# but may be absent when called directly from Nextflow/Snakemake)
	mkdir -p "$RAW_DIR_ROOT"
	parallel \
		--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
		--env _CONDA_PROFILE_SCRIPT --env _SHARED_GZIP_C --env WF_MANAGED_ENV \
		--env RAW_DIR_ROOT --env TRIM_DIR_ROOT --env THREADS_PER_JOB \
		--env LOG_FILE --env ERROR_WARN_FILE \
		-j "${JOBS:-2}" \
		--halt soon,fail,1 \
		--joblog "$RAW_DIR_ROOT/parallel_download_${BASHPID:-$$}.log" \
		_download_worker {} \
		< <(printf "%s\n" "${SRR_LIST[@]}")
}
