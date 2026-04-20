#!/bin/bash
# ==============================================================================
# TRIMMING FUNCTIONS
# ==============================================================================
# Read trimming utilities using TrimGalore and Trimmomatic
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${TRIMMING_SOURCED:-}" == "true" ]] && return 0
TRIMMING_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/a_preprocessing}"
if [[ -z "$SCRIPT_DIR" ]]; then
	SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
	SCRIPT_DIR="$(cd "$SCRIPT_DIR" 2>/dev/null && pwd)"
fi
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

# Reuse pigz detection from shared_utils_preproc.sh (already sourced above)
# _SHARED_HAS_PIGZ is set at module load in shared_utils_preproc.sh
_TRIMMING_HAS_PIGZ="${_SHARED_HAS_PIGZ:-false}"
export _TRIMMING_HAS_PIGZ

# ==============================================================================
# PRIMARY: TrimGalore → Trimmomatic HEADCROP → optional TAILCROP
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
	
	local -a _decompress_cmd=(gunzip)
	[[ "$_TRIMMING_HAS_PIGZ" == "true" ]] && _decompress_cmd=(pigz -d -p "${THREADS_PER_JOB:-4}")

	if [[ -n "$raw2" && -f "$raw2" ]]; then
		# ── Paired-end ──
		log_info "Trimming $SRR (PE) with TrimGalore..."
		run_with_space_time_log trim_galore --cores "${THREADS_PER_JOB:-2}" \
			--paired "$raw1" "$raw2" --output_dir "$trim_dir"

		log_info "Applying HEADCROP:${HEADCROP_BASES} for $SRR..."
		local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
		local tg_r2="$trim_dir/${SRR}_2_val_2.fq"
		if [[ -f "${tg_r1}.gz" && -f "${tg_r2}.gz" ]]; then
			"${_decompress_cmd[@]}" "${tg_r1}.gz" &
			local _pid1=$!
			"${_decompress_cmd[@]}" "${tg_r2}.gz" &
			local _pid2=$!
			wait "$_pid1" "$_pid2"
		elif [[ -f "${tg_r1}.gz" ]]; then
			"${_decompress_cmd[@]}" "${tg_r1}.gz"
		elif [[ -f "${tg_r2}.gz" ]]; then
			"${_decompress_cmd[@]}" "${tg_r2}.gz"
		fi

		local tmp_r1="$trim_dir/${SRR}_1_headcrop.fq"
		local tmp_r2="$trim_dir/${SRR}_2_headcrop.fq"
		run_with_space_time_log trimmomatic PE -threads "${THREADS_PER_JOB:-2}" \
			"$tg_r1" "$tg_r2" "$tmp_r1" /dev/null "$tmp_r2" /dev/null \
			HEADCROP:${HEADCROP_BASES}
		mv "$tmp_r1" "$tg_r1"
		mv "$tmp_r2" "$tg_r2"

		if [[ "${TAILCROP_BASES}" -gt 0 ]]; then
			log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
			if run_with_error_capture cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} \
				-o "${tg_r1}.tmp" -p "${tg_r2}.tmp" "$tg_r1" "$tg_r2"; then
				mv "${tg_r1}.tmp" "$tg_r1"
				mv "${tg_r2}.tmp" "$tg_r2"
			else
				log_warn "TAILCROP failed for $SRR — continuing with headcropped files"
				rm -f "${tg_r1}.tmp" "${tg_r2}.tmp"
			fi
		fi

		verify_trimming_and_cleanup "$SRR" "$tg_r1" "$tg_r2" "$raw1" "$raw2"
	else
		# ── Single-end ──
		log_info "Trimming $SRR (SE) with TrimGalore..."
		run_with_space_time_log trim_galore --cores "${THREADS_PER_JOB:-2}" \
			"$raw1" --output_dir "$trim_dir"

		log_info "Applying HEADCROP:${HEADCROP_BASES} for $SRR..."
		local tg_r1="$trim_dir/${SRR}_trimmed.fq"
		[[ -f "${tg_r1}.gz" ]] && "${_decompress_cmd[@]}" "${tg_r1}.gz"

		local tmp_r1="$trim_dir/${SRR}_headcrop.fq"
		run_with_space_time_log trimmomatic SE -threads "${THREADS_PER_JOB:-2}" \
			"$tg_r1" "$tmp_r1" HEADCROP:${HEADCROP_BASES}
		mv "$tmp_r1" "$tg_r1"

		if [[ "${TAILCROP_BASES}" -gt 0 ]]; then
			log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
			if run_with_error_capture cutadapt -u -${TAILCROP_BASES} \
				-o "${tg_r1}.tmp" "$tg_r1"; then
				mv "${tg_r1}.tmp" "$tg_r1"
			else
				log_warn "TAILCROP failed for $SRR — continuing with headcropped file"
				rm -f "${tg_r1}.tmp"
			fi
		fi

		verify_trimming_and_cleanup "$SRR" "$tg_r1" "" "$raw1" ""
	fi
}

# ==============================================================================
# SEQUENTIAL TRIMMING (primary entry points)
# ==============================================================================

# trim_srrs() — removed (dead code; superseded by download_and_trim_srrs_parallel)

# ==============================================================================
# ALTERNATIVE: Trimmomatic-only (no TrimGalore pre-step)
# ==============================================================================

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
		
		if [[ -n "$raw2" && -f "$raw2" ]]; then
			# ── Paired-end ──
			run_with_space_time_log trimmomatic PE -threads "${THREADS_PER_JOB:-${THREADS:-4}}" \
				"$raw1" "$raw2" "$out1" /dev/null "$out2" /dev/null \
				ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
				HEADCROP:${HEADCROP_BASES} SLIDINGWINDOW:${SW_SIZE}:${SW_QUAL} MINLEN:${MINLEN}

			# Trim last N bases using cutadapt if TAILCROP > 0
			if [[ "${TAILCROP_BASES}" -gt 0 ]]; then
				log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
				if run_with_error_capture cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} \
					-o "${out1}.tmp" -p "${out2}.tmp" "$out1" "$out2"; then
					mv "${out1}.tmp" "$out1"
					mv "${out2}.tmp" "$out2"
				else
					log_warn "TAILCROP failed for $SRR — continuing with trimmed files"
					rm -f "${out1}.tmp" "${out2}.tmp"
				fi
			fi

			verify_trimming_and_cleanup "$SRR" "$out1" "$out2" "$raw1" "$raw2"
		else
			# ── Single-end ──
			local se_out="$trim_dir/${SRR}_trimmed.fq"
			run_with_space_time_log trimmomatic SE -threads "${THREADS_PER_JOB:-${THREADS:-4}}" \
				"$raw1" "$se_out" \
				ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
				HEADCROP:${HEADCROP_BASES} SLIDINGWINDOW:${SW_SIZE}:${SW_QUAL} MINLEN:${MINLEN}

			if [[ "${TAILCROP_BASES}" -gt 0 ]]; then
				log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
				if run_with_error_capture cutadapt -u -${TAILCROP_BASES} \
					-o "${se_out}.tmp" "$se_out"; then
					mv "${se_out}.tmp" "$se_out"
				else
					log_warn "TAILCROP failed for $SRR — continuing with trimmed file"
					rm -f "${se_out}.tmp"
				fi
			fi

			verify_trimming_and_cleanup "$SRR" "$se_out" "" "$raw1" ""
		fi
	done
	gzip_trimmed_fastq_files
}

# ==============================================================================
# SERIALIZATION HELPER
# ==============================================================================

# Serialize SRR_TRIM_PROFILE_MAP associative array to a semicolon-delimited string
# for export to GNU Parallel subshells (associative arrays can't be exported).
# Usage: export SERIALIZED_TRIM_PROFILES="$(_serialize_trim_profiles)"
_serialize_trim_profiles() {
	local -a _prof_parts=()
	for key in "${!SRR_TRIM_PROFILE_MAP[@]}"; do
		_prof_parts+=("${key}=${SRR_TRIM_PROFILE_MAP[$key]}")
	done
	# Join with IFS directly — avoids printf subprocess fork
	local IFS=';'
	echo "${_prof_parts[*]}"
}

# ==============================================================================
# COMBINED AND PARALLEL VARIANTS
# ==============================================================================
# These functions combine download+trim or run either method in parallel via GNU Parallel.

# Combined download + trim in a single pass (sequential)

download_and_trim_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided"; return 1; }

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
			local -a _ccmd=(gzip)
			[[ "$_TRIMMING_HAS_PIGZ" == "true" ]] && _ccmd=(pigz -p "${THREADS:-4}")
			local _cw1=0 _cw2=0
			[[ -f "$raw_dir/${SRR}_1.fastq" ]] && { "${_ccmd[@]}" "$raw_dir/${SRR}_1.fastq" & _cw1=$!; }
			[[ -f "$raw_dir/${SRR}_2.fastq" ]] && { "${_ccmd[@]}" "$raw_dir/${SRR}_2.fastq" & _cw2=$!; }
			[[ "$_cw1" -ne 0 ]] && wait "$_cw1"
			[[ "$_cw2" -ne 0 ]] && wait "$_cw2"
			find_raw_fastq "$SRR"
		fi

		[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; continue; }
		_trim_single_srr "$SRR"
	done
	gzip_trimmed_fastq_files
}

# Combined download + trim with GNU Parallel

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
	
	export SERIALIZED_TRIM_PROFILES="$(_serialize_trim_profiles)"

	# Export _log_impl (core logger) alongside its callers — without it, log_info/log_warn/log_error
	# fail silently in GNU Parallel subshells because they delegate to _log_impl.
	# Also export strip_ansi_stream and capture_stderr_errors — transitive deps of
	# run_with_space_time_log (pipes to strip_ansi_stream) and run_with_error_capture
	# (pipes to capture_stderr_errors).
	export -f _log_impl timestamp log log_info log_warn log_error run_with_space_time_log run_with_error_capture
	export -f strip_ansi_stream capture_stderr_errors
	export -f find_trimmed_fastq find_raw_fastq verify_trimming_and_cleanup

	_parallel_worker() {
		local SRR="$1"

		# Activate conda environment in subshell (skip when orchestrator manages env)
		if [[ -z "${WF_MANAGED_ENV:-}" && -n "${CONDA_PREFIX:-}" ]]; then
			source "${_CONDA_PROFILE_SCRIPT:-${CONDA_EXE%/*}/../etc/profile.d/conda.sh}" 2>/dev/null || true
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
			local -a _ccmd=(gzip)
			[[ "$_TRIMMING_HAS_PIGZ" == "true" ]] && _ccmd=(pigz -p "${THREADS_PER_JOB:-4}")
			local _cw1=0 _cw2=0
			[[ -f "$raw_dir/${SRR}_1.fastq" ]] && { "${_ccmd[@]}" "$raw_dir/${SRR}_1.fastq" & _cw1=$!; }
			[[ -f "$raw_dir/${SRR}_2.fastq" ]] && { "${_ccmd[@]}" "$raw_dir/${SRR}_2.fastq" & _cw2=$!; }
			[[ "$_cw1" -ne 0 ]] && wait "$_cw1"
			[[ "$_cw2" -ne 0 ]] && wait "$_cw2"
			find_raw_fastq "$SRR"
		fi

		[[ -z "$raw1" ]] && { log_warn "No raw for $SRR"; return 1; }

		# Deserialize trim profiles: O(1) awk lookup replaces O(n) while-read scan
		local profile
		profile=$(awk -F'=' -v srr="$SRR" 'BEGIN{RS=";"} $1==srr{print $2; exit}' <<< "$SERIALIZED_TRIM_PROFILES")
		profile="${profile:-$TRIM_PROFILE_DEFAULT}"

		local HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL
		IFS=':' read -r HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL <<< "$profile"

		local -a _dcmd=(gunzip)
		[[ "$_TRIMMING_HAS_PIGZ" == "true" ]] && _dcmd=(pigz -d -p "${THREADS_PER_JOB:-2}")

		if [[ -n "$raw2" && -f "$raw2" ]]; then
			# ── Paired-end ──
			trim_galore --cores "${THREADS_PER_JOB:-2}" --paired "$raw1" "$raw2" --output_dir "$trim_dir"
			local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
			local tg_r2="$trim_dir/${SRR}_2_val_2.fq"
			if [[ -f "${tg_r1}.gz" && -f "${tg_r2}.gz" ]]; then
				"${_dcmd[@]}" "${tg_r1}.gz" &
				local _p1=$!
				"${_dcmd[@]}" "${tg_r2}.gz" &
				local _p2=$!
				wait "$_p1" "$_p2"
			else
				[[ -f "${tg_r1}.gz" ]] && "${_dcmd[@]}" "${tg_r1}.gz"
				[[ -f "${tg_r2}.gz" ]] && "${_dcmd[@]}" "${tg_r2}.gz"
			fi

			trimmomatic PE -threads "${THREADS_PER_JOB:-2}" "$tg_r1" "$tg_r2" \
				"${tg_r1}.tmp" /dev/null "${tg_r2}.tmp" /dev/null HEADCROP:${HEADCROP_BASES}
			mv "${tg_r1}.tmp" "$tg_r1"
			mv "${tg_r2}.tmp" "$tg_r2"

			if [[ "${TAILCROP_BASES:-0}" -gt 0 ]]; then
				log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
				if cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} \
					-o "${tg_r1}.tmp" -p "${tg_r2}.tmp" "$tg_r1" "$tg_r2"; then
					mv "${tg_r1}.tmp" "$tg_r1"
					mv "${tg_r2}.tmp" "$tg_r2"
				else
					log_warn "cutadapt TAILCROP failed for $SRR (exit $?) — keeping un-tailcropped files"
					rm -f "${tg_r1}.tmp" "${tg_r2}.tmp"
				fi
			fi

			verify_trimming_and_cleanup "$SRR" "$tg_r1" "$tg_r2" "$raw1" "$raw2"
		else
			# ── Single-end ──
			trim_galore --cores "${THREADS_PER_JOB:-2}" "$raw1" --output_dir "$trim_dir"
			local tg_r1="$trim_dir/${SRR}_trimmed.fq"
			[[ -f "${tg_r1}.gz" ]] && "${_dcmd[@]}" "${tg_r1}.gz"

			trimmomatic SE -threads "${THREADS_PER_JOB:-2}" "$tg_r1" \
				"${tg_r1}.tmp" HEADCROP:${HEADCROP_BASES}
			mv "${tg_r1}.tmp" "$tg_r1"

			if [[ "${TAILCROP_BASES:-0}" -gt 0 ]]; then
				log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
				if cutadapt -u -${TAILCROP_BASES} -o "${tg_r1}.tmp" "$tg_r1"; then
					mv "${tg_r1}.tmp" "$tg_r1"
				else
					log_warn "cutadapt TAILCROP failed for $SRR (exit $?) — keeping un-tailcropped file"
					rm -f "${tg_r1}.tmp"
				fi
			fi

			verify_trimming_and_cleanup "$SRR" "$tg_r1" "" "$raw1" ""
		fi
	}
	export -f _parallel_worker

	# Ensure joblog parent dir exists (created by orchestrator in bash mode,
	# but may be absent when called directly from Nextflow/Snakemake)
	mkdir -p "$TRIM_DIR_ROOT"
	parallel \
		--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
		--env _CONDA_PROFILE_SCRIPT --env WF_MANAGED_ENV \
		--env RAW_DIR_ROOT --env TRIM_DIR_ROOT --env THREADS_PER_JOB \
		--env LOG_FILE --env ERROR_WARN_FILE --env _SHARED_GZIP_C --env _SHARED_GZIP_DC \
		--env _TRIMMING_HAS_PIGZ \
		-j "${JOBS:-2}" \
		--halt soon,fail,1 \
		--joblog "$TRIM_DIR_ROOT/parallel_trim_galore_${BASHPID:-$$}.log" \
		_parallel_worker {} \
		< <(printf "%s\n" "${SRR_LIST[@]}")
	gzip_trimmed_fastq_files
}

# Trimmomatic-only with GNU Parallel
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
	
	export SERIALIZED_TRIM_PROFILES="$(_serialize_trim_profiles)"

	# Export _log_impl (core logger) alongside its callers — without it, log_info/log_warn/log_error
	# fail silently in GNU Parallel subshells because they delegate to _log_impl.
	# Also export strip_ansi_stream and capture_stderr_errors — transitive deps of
	# run_with_space_time_log (pipes to strip_ansi_stream) and run_with_error_capture
	# (pipes to capture_stderr_errors).
	export -f _log_impl timestamp log log_info log_warn log_error run_with_space_time_log run_with_error_capture
	export -f strip_ansi_stream capture_stderr_errors
	export -f find_trimmed_fastq find_raw_fastq verify_trimming_and_cleanup

	_trimmomatic_parallel_worker() {
		local SRR="$1"

		# Activate conda environment in subshell (skip when orchestrator manages env)
		if [[ -z "${WF_MANAGED_ENV:-}" && -n "${CONDA_PREFIX:-}" ]]; then
			source "${_CONDA_PROFILE_SCRIPT:-${CONDA_EXE%/*}/../etc/profile.d/conda.sh}" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"

		mkdir -p "$trim_dir"

		find_trimmed_fastq "$SRR"
		[[ -n "${trimmed1:-}" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; return 0; }

		find_raw_fastq "$SRR"
		[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; return 1; }

		# Deserialize trim profiles: O(1) awk lookup replaces O(n) while-read scan
		local profile
		profile=$(awk -F'=' -v srr="$SRR" 'BEGIN{RS=";"} $1==srr{print $2; exit}' <<< "$SERIALIZED_TRIM_PROFILES")
		profile="${profile:-$TRIM_PROFILE_DEFAULT}"

		IFS=':' read -r HEADCROP_BASES TAILCROP_BASES MINLEN SW_SIZE SW_QUAL <<< "$profile"

		log_info "Trimming $SRR with Trimmomatic (HEADCROP:$HEADCROP_BASES, TAILCROP:$TAILCROP_BASES, MINLEN:$MINLEN, SW:$SW_SIZE:$SW_QUAL)..."

		if [[ -n "$raw2" && -f "$raw2" ]]; then
			# ── Paired-end ──
			local out1="$trim_dir/${SRR}_1_val_1.fq"
			local out2="$trim_dir/${SRR}_2_val_2.fq"
			run_with_space_time_log trimmomatic PE -threads "${THREADS_PER_JOB:-${THREADS:-4}}" \
				"$raw1" "$raw2" "$out1" /dev/null "$out2" /dev/null \
				ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
				HEADCROP:${HEADCROP_BASES} SLIDINGWINDOW:${SW_SIZE}:${SW_QUAL} MINLEN:${MINLEN}

			if [[ "${TAILCROP_BASES:-0}" -gt 0 ]]; then
				log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
				if run_with_error_capture cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} \
					-o "${out1}.tmp" -p "${out2}.tmp" "$out1" "$out2"; then
					mv "${out1}.tmp" "$out1"
					mv "${out2}.tmp" "$out2"
				else
					log_warn "cutadapt TAILCROP failed for $SRR (exit $?) — keeping un-tailcropped files"
					rm -f "${out1}.tmp" "${out2}.tmp"
				fi
			fi

			verify_trimming_and_cleanup "$SRR" "$out1" "$out2" "$raw1" "$raw2"
		else
			# ── Single-end: use _trimmed.fq to match find_trimmed_fastq() convention ──
			local se_out="$trim_dir/${SRR}_trimmed.fq"
			run_with_space_time_log trimmomatic SE -threads "${THREADS_PER_JOB:-${THREADS:-4}}" \
				"$raw1" "$se_out" \
				ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
				HEADCROP:${HEADCROP_BASES} SLIDINGWINDOW:${SW_SIZE}:${SW_QUAL} MINLEN:${MINLEN}

			if [[ "${TAILCROP_BASES:-0}" -gt 0 ]]; then
				log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
				if run_with_error_capture cutadapt -u -${TAILCROP_BASES} -o "${se_out}.tmp" "$se_out"; then
					mv "${se_out}.tmp" "$se_out"
				else
					log_warn "cutadapt TAILCROP failed for $SRR (exit $?) — keeping un-tailcropped file"
					rm -f "${se_out}.tmp"
				fi
			fi

			verify_trimming_and_cleanup "$SRR" "$se_out" "" "$raw1" ""
		fi
	}
	export -f _trimmomatic_parallel_worker

	# Ensure joblog parent dir exists (see download.sh comment)
	mkdir -p "$TRIM_DIR_ROOT"
	parallel \
		--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
		--env _CONDA_PROFILE_SCRIPT --env WF_MANAGED_ENV \
		--env RAW_DIR_ROOT --env TRIM_DIR_ROOT --env THREADS_PER_JOB \
		--env LOG_FILE --env ERROR_WARN_FILE --env _SHARED_GZIP_C --env _SHARED_GZIP_DC \
		-j "${JOBS:-2}" \
		--halt soon,fail,1 \
		--joblog "$TRIM_DIR_ROOT/parallel_trimmomatic_${BASHPID:-$$}.log" \
		_trimmomatic_parallel_worker {} \
		< <(printf "%s\n" "${SRR_LIST[@]}")
	gzip_trimmed_fastq_files
}
