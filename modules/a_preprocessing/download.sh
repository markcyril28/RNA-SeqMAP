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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_preproc.sh"

# ==============================================================================
# CORE DOWNLOAD FUNCTION
# ==============================================================================

download_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided for download"; return 1; }
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir"
		
		# Check if trimmed files exist
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed files for $SRR exist. Skipping download."; continue; }
		
		# Check if raw files exist
		find_raw_fastq "$SRR"
		[[ -n "$raw1" ]] && { log_info "Raw files for $SRR exist. Skipping download."; continue; }
		
		log_info "Downloading $SRR..."
		run_with_space_time_log prefetch "$SRR" --output-directory "$raw_dir"
		run_with_space_time_log fasterq-dump --split-files --threads "$THREADS" \
			"$raw_dir/$SRR/$SRR.sra" -O "$raw_dir"
		
		# Compress downloaded files
		[[ -f "$raw_dir/${SRR}_1.fastq" ]] && gzip "$raw_dir/${SRR}_1.fastq" "$raw_dir/${SRR}_2.fastq"
	done
	log_info "All downloads completed."
}

# ==============================================================================
# DOWNLOAD WITH WGET (ENA)
# ==============================================================================

download_srrs_wget() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed $SRR exists. Skipping."; continue; }
		
		find_raw_fastq "$SRR"
		[[ -n "$raw1" ]] && { log_info "Raw $SRR exists. Skipping."; continue; }
		
		# Get ENA links
		local ena_links=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2)
		[[ -z "$ena_links" ]] && { log_warn "No ENA links for $SRR"; continue; }
		
		log_info "Downloading $SRR from ENA..."
		wget -q -c -P "$raw_dir" "ftp://$(echo "$ena_links" | cut -f1 -d';')" "ftp://$(echo "$ena_links" | cut -f2 -d';')" || {
			log_warn "ENA download failed for $SRR, trying SRA..."
			download_srrs "$SRR"
		}
	done
}

# ==============================================================================
# DOWNLOAD WITH KINGFISHER
# ==============================================================================

download_srrs_kingfisher() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }

	# Kingfisher requires download methods; allow override via KINGFISHER_METHODS env, default to prefetch then ENA FTP
	local KF_METHODS="${KINGFISHER_METHODS:-prefetch,ena-ftp}"
	
	if ! command -v kingfisher >/dev/null 2>&1; then
		log_warn "Kingfisher not installed. Falling back to prefetch."
		download_srrs "${SRR_LIST[@]}"
		return $?
	fi
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed $SRR exists. Skipping."; continue; }
		
		find_raw_fastq "$SRR"
		[[ -n "$raw1" ]] && { log_info "Raw $SRR exists. Skipping."; continue; }
		
		log_info "Downloading $SRR with kingfisher..."
		run_with_space_time_log kingfisher get --run-identifiers "$SRR" --output-directory "$raw_dir" \
			--download-threads "$THREADS" --extraction-threads "$THREADS" --gzip --check-md5sums \
			--download-methods "$KF_METHODS" || {
			log_warn "Kingfisher failed for $SRR, trying prefetch..."
			download_srrs "$SRR"
		}
	done
}

# ==============================================================================
# PARALLEL DOWNLOAD FUNCTIONS
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
	export RAW_DIR_ROOT THREADS THREADS_PER_JOB
	export -f log_info log_warn log_error find_trimmed_fastq find_raw_fastq
	
	_download_worker() {
		local SRR="$1"
		
		# Activate conda environment in subshell
		if [[ -n "$CONDA_PREFIX" ]]; then
			source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && return 0
		
		find_raw_fastq "$SRR"
		[[ -n "$raw1" ]] && return 0
		
		prefetch "$SRR" --output-directory "$raw_dir" || return 1
		fasterq-dump --split-files --threads "${THREADS_PER_JOB:-2}" "$raw_dir/$SRR/$SRR.sra" -O "$raw_dir" || return 1
		[[ -f "$raw_dir/${SRR}_1.fastq" ]] && gzip "$raw_dir/${SRR}_1.fastq" "$raw_dir/${SRR}_2.fastq"
	}
	export -f _download_worker
	
	printf "%s\n" "${SRR_LIST[@]}" | parallel --env PATH --env CONDA_PREFIX -j "${JOBS:-2}" _download_worker {}
}
