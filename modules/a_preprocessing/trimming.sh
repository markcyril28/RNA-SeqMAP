#!/bin/bash
# ==============================================================================
# TRIMMING FUNCTIONS
# ==============================================================================
# Read trimming utilities using TrimGalore and Trimmomatic
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${TRIMMING_SOURCED:-}" == "true" ]] && return 0
export TRIMMING_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_preproc.sh"

# ==============================================================================
# TRIMMING CONFIGURATION - IMPORTANT PARAMETERS AT TOP
# ==============================================================================

# Default values (overridden by profile-based settings)
HEADCROP_BASES="${HEADCROP_BASES:-10}"
TAILCROP_BASES="${TAILCROP_BASES:-0}"
MINLEN="${MINLEN:-36}"
SW_SIZE="${SW_SIZE:-4}"
SW_QUAL="${SW_QUAL:-20}"

# ==============================================================================
# CORE TRIMMING LOGIC
# ==============================================================================

# Internal function for trimming a single SRR
_trim_single_srr() {
	local SRR="$1"
	local raw_dir="$RAW_DIR_ROOT/$SRR"
	local trim_dir="$TRIM_DIR_ROOT/$SRR"
	
	mkdir -p "$trim_dir"
	
	# Skip if already trimmed
	find_trimmed_fastq "$SRR"
	[[ -n "${trimmed1:-}" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; return 0; }
	
	# Find raw files
	find_raw_fastq "$SRR"
	[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; return 1; }
	
	# Get trim parameters for this specific SRR from profile
	get_trim_params "$SRR"
	
	log_info "Trimming $SRR with TrimGalore..."
	run_with_space_time_log trim_galore --cores "$THREADS_PER_JOB" \
		--paired "$raw1" "$raw2" --output_dir "$trim_dir"
	
	# Trimmomatic HEADCROP with profile-based parameters
	log_info "Applying HEADCROP:${HEADCROP_BASES} for $SRR..."
	local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
	local tg_r2="$trim_dir/${SRR}_2_val_2.fq"
	[[ -f "${tg_r1}.gz" ]] && gunzip "${tg_r1}.gz"
	[[ -f "${tg_r2}.gz" ]] && gunzip "${tg_r2}.gz"
	
	local tmp_r1="$trim_dir/${SRR}_1_headcrop.fq"
	local tmp_r2="$trim_dir/${SRR}_2_headcrop.fq"
	run_with_space_time_log trimmomatic PE -threads "$THREADS" \
		"$tg_r1" "$tg_r2" "$tmp_r1" /dev/null "$tmp_r2" /dev/null \
		HEADCROP:${HEADCROP_BASES}
	mv "$tmp_r1" "$tg_r1"
	mv "$tmp_r2" "$tg_r2"
	
	# Apply TAILCROP if > 0
	if [[ "${TAILCROP_BASES}" -gt 0 ]]; then
		log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
		run_with_error_capture cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} \
			-o "${tg_r1}.tmp" -p "${tg_r2}.tmp" "$tg_r1" "$tg_r2"
		mv "${tg_r1}.tmp" "$tg_r1"
		mv "${tg_r2}.tmp" "$tg_r2"
	fi
	
	verify_trimming_and_cleanup "$SRR" "$tg_r1" "$tg_r2" "$raw1" "$raw2"
}

# ==============================================================================
# MAIN TRIMMING FUNCTIONS
# ==============================================================================

trim_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided for trimming"; return 1; }
	
	local total=${#SRR_LIST[@]} current=0
	for SRR in "${SRR_LIST[@]}"; do
		((current++))
		log_info "[$current/$total] Processing $SRR..."
		_trim_single_srr "$SRR" || log_warn "Trimming failed for $SRR"
	done
	gzip_trimmed_fastq_files
	log_info "All trimming completed."
}

trim_srrs_trimmomatic() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		local out1="$trim_dir/${SRR}_1_val_1.fq"
		local out2="$trim_dir/${SRR}_2_val_2.fq"
		
		mkdir -p "$trim_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "${trimmed1:-}" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; continue; }
		
		find_raw_fastq "$SRR"
		[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; continue; }
		
		# Get trim parameters for this specific SRR from profile
		get_trim_params "$SRR"
		log_info "Trimming $SRR with Trimmomatic (HEADCROP:$HEADCROP_BASES, TAILCROP:$TAILCROP_BASES, MINLEN:$MINLEN, SW:$SW_SIZE:$SW_QUAL)..."
		
		run_with_space_time_log trimmomatic PE -threads "$THREADS" \
			"$raw1" "$raw2" "$out1" /dev/null "$out2" /dev/null \
			ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
			HEADCROP:${HEADCROP_BASES} SLIDINGWINDOW:${SW_SIZE}:${SW_QUAL} MINLEN:${MINLEN}
		
		# Trim last N bases using cutadapt if TAILCROP > 0
		if [[ "${TAILCROP_BASES}" -gt 0 ]]; then
			log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
			run_with_error_capture cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} \
				-o "${out1}.tmp" -p "${out2}.tmp" "$out1" "$out2"
			mv "${out1}.tmp" "$out1"
			mv "${out2}.tmp" "$out2"
		fi
		
		verify_trimming_and_cleanup "$SRR" "$out1" "$out2" "$raw1" "$raw2"
	done
	gzip_trimmed_fastq_files
}

# ==============================================================================
# COMBINED DOWNLOAD AND TRIM
# ==============================================================================

download_and_trim_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }
	
	# Source download functions
	source "$SCRIPT_DIR/download.sh"
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir" "$trim_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; continue; }
		
		find_raw_fastq "$SRR"
		if [[ -z "$raw1" ]]; then
			log_info "Downloading $SRR..."
			run_with_space_time_log prefetch "$SRR" --output-directory "$raw_dir"
			run_with_space_time_log fasterq-dump --split-files --threads "$THREADS" \
				"$raw_dir/$SRR/$SRR.sra" -O "$raw_dir"
			[[ -f "$raw_dir/${SRR}_1.fastq" ]] && gzip "$raw_dir/${SRR}_1.fastq" "$raw_dir/${SRR}_2.fastq"
			find_raw_fastq "$SRR"
		fi
		
		[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; continue; }
		_trim_single_srr "$SRR"
	done
	gzip_trimmed_fastq_files
}

# ==============================================================================
# PARALLEL TRIMMING
# ==============================================================================

download_and_trim_srrs_parallel() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }
	
	if ! should_use_parallel; then
		log_info "Running download_and_trim sequentially (USE_GNU_PARALLEL=${USE_GNU_PARALLEL:-FALSE})"
		download_and_trim_srrs "${SRR_LIST[@]}"
		return $?
	fi
	
	log_info "Running download_and_trim with GNU Parallel (JOBS=${JOBS:-2})"
	
	# Export PATH and conda environment so tools are available in parallel subshells
	export PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE
	
	# Export all required variables including TIME_DIR and log paths
	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS THREADS_PER_JOB JOBS LOG_FILE
	export TIME_DIR TIME_FILE TIME_TEMP SPACE_TIME_FILE ERROR_WARN_FILE
	export TRIM_PROFILE_DEFAULT DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING
	
	# Serialize SRR_TRIM_PROFILE_MAP to a string for export (associative arrays can't be exported)
	local serialized_profiles=""
	for key in "${!SRR_TRIM_PROFILE_MAP[@]}"; do
		serialized_profiles+="${key}=${SRR_TRIM_PROFILE_MAP[$key]};"
	done
	export SERIALIZED_TRIM_PROFILES="$serialized_profiles"
	
	export -f timestamp log log_info log_warn log_error run_with_space_time_log
	export -f find_trimmed_fastq find_raw_fastq verify_trimming_and_cleanup
	
	_parallel_worker() {
		local SRR="$1"
		
		# Activate conda environment in subshell
		if [[ -n "$CONDA_PREFIX" ]]; then
			source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir" "$trim_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed $SRR exists. Skipping."; return 0; }
		
		find_raw_fastq "$SRR"
		if [[ -z "$raw1" ]]; then
			prefetch "$SRR" --output-directory "$raw_dir" || return 1
			fasterq-dump --split-files --threads "${THREADS_PER_JOB:-4}" "$raw_dir/$SRR/$SRR.sra" -O "$raw_dir" || return 1
			[[ -f "$raw_dir/${SRR}_1.fastq" ]] && gzip "$raw_dir/${SRR}_1.fastq" "$raw_dir/${SRR}_2.fastq"
			find_raw_fastq "$SRR"
		fi
		
		[[ -z "$raw1" ]] && { log_warn "No raw for $SRR"; return 1; }
		
		trim_galore --cores "${THREADS_PER_JOB:-2}" --paired "$raw1" "$raw2" --output_dir "$trim_dir"
		local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
		local tg_r2="$trim_dir/${SRR}_2_val_2.fq"
		[[ -f "${tg_r1}.gz" ]] && gunzip "${tg_r1}.gz" "${tg_r2}.gz"
		
		# Deserialize trim profiles and get HEADCROP for this SRR
		local profile="$TRIM_PROFILE_DEFAULT"
		while IFS='=' read -r key val; do
			[[ "$key" == "$SRR" ]] && { profile="$val"; break; }
		done < <(echo "$SERIALIZED_TRIM_PROFILES" | tr ';' '\n')
		
		local HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL
		IFS=':' read -r HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL <<< "$profile"
		
		trimmomatic PE -threads "${THREADS_PER_JOB:-2}" "$tg_r1" "$tg_r2" \
			"${tg_r1}.tmp" /dev/null "${tg_r2}.tmp" /dev/null HEADCROP:${HEADCROP_BASES}
		mv "${tg_r1}.tmp" "$tg_r1"
		mv "${tg_r2}.tmp" "$tg_r2"
		
		verify_trimming_and_cleanup "$SRR" "$tg_r1" "$tg_r2" "$raw1" "$raw2"
	}
	export -f _parallel_worker
	
	printf "%s\n" "${SRR_LIST[@]}" | parallel --env PATH --env CONDA_PREFIX -j "${JOBS:-2}" _parallel_worker {}
	gzip_trimmed_fastq_files
}

trim_srrs_trimmomatic_parallel() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }

	if ! should_use_parallel; then
		log_info "Running trim_srrs_trimmomatic sequentially (USE_GNU_PARALLEL=${USE_GNU_PARALLEL:-FALSE})"
		trim_srrs_trimmomatic "${SRR_LIST[@]}"
		return $?
	fi

	log_info "Running trim_srrs_trimmomatic with GNU Parallel (JOBS=${JOBS:-2})"
	
	# Export PATH and conda environment so tools are available in parallel subshells
	export PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE
	
	# Export all required variables including TIME_DIR and log paths
	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS THREADS_PER_JOB JOBS LOG_FILE
	export TIME_DIR TIME_FILE TIME_TEMP SPACE_TIME_FILE ERROR_WARN_FILE
	export TRIM_PROFILE_DEFAULT DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING
	
	# Serialize SRR_TRIM_PROFILE_MAP to a string for export (associative arrays can't be exported)
	local serialized_profiles=""
	for key in "${!SRR_TRIM_PROFILE_MAP[@]}"; do
		serialized_profiles+="${key}=${SRR_TRIM_PROFILE_MAP[$key]};"
	done
	export SERIALIZED_TRIM_PROFILES="$serialized_profiles"
	
	export -f timestamp log log_info log_warn log_error run_with_space_time_log run_with_error_capture
	export -f find_trimmed_fastq find_raw_fastq verify_trimming_and_cleanup

	_trimmomatic_parallel_worker() {
		local SRR="$1"
		
		# Activate conda environment in subshell
		if [[ -n "$CONDA_PREFIX" ]]; then
			source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		local out1="$trim_dir/${SRR}_1_val_1.fq"
		local out2="$trim_dir/${SRR}_2_val_2.fq"

		mkdir -p "$trim_dir"

		find_trimmed_fastq "$SRR"
		[[ -n "${trimmed1:-}" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; return 0; }

		find_raw_fastq "$SRR"
		[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; return 1; }

		# Deserialize trim profiles and get params for this SRR
		local profile="$TRIM_PROFILE_DEFAULT"
		while IFS='=' read -r key val; do
			[[ "$key" == "$SRR" ]] && { profile="$val"; break; }
		done < <(echo "$SERIALIZED_TRIM_PROFILES" | tr ';' '\n')
		
		IFS=':' read -r HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL <<< "$profile"
		
		log_info "Trimming $SRR with Trimmomatic (HEADCROP:$HEADCROP_BASES, TAILCROP:$TAILCROP_BASES, MINLEN:$MINLEN, SW:$SW_SIZE:$SW_QUAL)..."

		run_with_space_time_log trimmomatic PE -threads "${THREADS_PER_JOB:-${THREADS:-4}}" \
			"$raw1" "$raw2" "$out1" /dev/null "$out2" /dev/null \
			ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
			HEADCROP:${HEADCROP_BASES} SLIDINGWINDOW:${SW_SIZE}:${SW_QUAL} MINLEN:${MINLEN}

		if [[ "${TAILCROP_BASES:-0}" -gt 0 ]]; then
			log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
			run_with_error_capture cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} \
				-o "${out1}.tmp" -p "${out2}.tmp" "$out1" "$out2"
			mv "${out1}.tmp" "$out1"
			mv "${out2}.tmp" "$out2"
		fi

		verify_trimming_and_cleanup "$SRR" "$out1" "$out2" "$raw1" "$raw2"
	}
	export -f _trimmomatic_parallel_worker

	printf "%s\n" "${SRR_LIST[@]}" | parallel --env PATH --env CONDA_PREFIX -j "${JOBS:-2}" _trimmomatic_parallel_worker {}
	gzip_trimmed_fastq_files
}
