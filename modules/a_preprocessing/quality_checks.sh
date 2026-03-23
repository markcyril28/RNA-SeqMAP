#!/bin/bash
# ==============================================================================
# QUALITY CONTROL FUNCTIONS
# ==============================================================================
# FastQC and MultiQC quality control utilities
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${QC_SOURCED:-}" == "true" ]] && return 0
export QC_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/a_preprocessing}"
if [[ -z "$SCRIPT_DIR" ]]; then SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."; fi
source "$SCRIPT_DIR/shared_utils_preproc.sh"

# Cache tool availability once at module load (avoids command -v per SRR)
_QC_HAS_FASTQC=false; command -v fastqc >/dev/null 2>&1 && _QC_HAS_FASTQC=true
_QC_HAS_MULTIQC=false; command -v multiqc >/dev/null 2>&1 && _QC_HAS_MULTIQC=true

# ==============================================================================
# QUALITY CONTROL FUNCTIONS
# ==============================================================================

rename_fastqc_outputs() {
	local outdir="$1"
	local label="$2"
	shift 2

	shopt -s nullglob
	for input in "$@"; do
		[[ -e "$input" ]] || continue

		# Pure bash parameter expansion — avoids basename subshell per file
		local base="${input##*/}"
		base="${base%.*}"
		[[ "$base" == *.fastq ]] && base="${base%.*}"
		[[ "$base" == *.fq ]] && base="${base%.*}"

		local clean_base="$base"
		clean_base="${clean_base/_val_1/}"
		clean_base="${clean_base/_val_2/}"
		clean_base="${clean_base/_val_3/}"

		for ext in html zip; do
			local src="$outdir/${base}_fastqc.${ext}"
			local dest="$outdir/${clean_base}_${label}_fastqc.${ext}"
			[[ -f "$src" ]] && mv "$src" "$dest"
		done

		local src_dir="$outdir/${base}_fastqc"
		local dest_dir="$outdir/${clean_base}_${label}_fastqc"
		[[ -d "$src_dir" ]] && mv "$src_dir" "$dest_dir"
	done
	shopt -u nullglob
}

run_quality_control() {
	local SRR="$1"
	local RAW_DIR="$RAW_DIR_ROOT/$SRR"
	local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
	local srr_outdir="$FASTQC_ROOT/$SRR"

	log_step "QC: $SRR"
	mkdir -p "$srr_outdir"

	if $_QC_HAS_FASTQC; then
		# QC for raw files (single folder per SRR)
		if [[ -d "$RAW_DIR" ]]; then
			local _raw_qc=("$srr_outdir"/*_raw_fastqc.html); if [[ ! -f "${_raw_qc[0]:-}" ]]; then
				log_info "Running FastQC on raw files for $SRR"
				run_with_space_time_log fastqc -t "${THREADS:-2}" -o "$srr_outdir" \
					"$RAW_DIR"/${SRR}*.fastq* 2>/dev/null || log_warn "FastQC failed for raw $SRR"
				rename_fastqc_outputs "$srr_outdir" "raw" "$RAW_DIR"/${SRR}*.fastq*
			fi
		fi

		# QC for trimmed files (same SRR folder, renamed as trimmed)
		if [[ -d "$TrimGalore_DIR" ]]; then
			local _trim_qc=("$srr_outdir"/*_trimmed_fastqc.html); if [[ ! -f "${_trim_qc[0]:-}" ]]; then
				log_info "Running FastQC on trimmed files for $SRR"
				run_with_space_time_log fastqc -t "${THREADS:-2}" -o "$srr_outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || log_warn "FastQC failed for trimmed $SRR"
				rename_fastqc_outputs "$srr_outdir" "trimmed" "$TrimGalore_DIR"/${SRR}*val*.fq*
			fi
		fi
	else
		log_warn "FastQC not found. Skipping QC for $SRR."
	fi
}

run_multiqc() {
	if $_QC_HAS_MULTIQC; then
		run_with_space_time_log multiqc "$FASTQC_ROOT" -o "$FASTQC_ROOT/summary" --force 2>/dev/null || true
	else
		log_warn "MultiQC not found. Skipping aggregation."
	fi
}

run_quality_control_all() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided for QC"; return 1; }
	
	log_step "Running Quality Control for ${#SRR_LIST[@]} samples"
	
	if should_use_parallel; then
		log_info "Running FastQC with GNU Parallel (JOBS=${JOBS:-2})"
		run_quality_control_parallel "${SRR_LIST[@]}"
	else
		log_info "Running FastQC sequentially (USE_GNU_PARALLEL=${USE_GNU_PARALLEL:-FALSE})"
		for SRR in "${SRR_LIST[@]}"; do
			run_quality_control "$SRR"
		done
	fi
	
	# Run MultiQC once at the end to aggregate all results
	run_multiqc
	
	log_info "Quality control completed for all samples."
}

# ==============================================================================
# PARALLEL QC FUNCTIONS
# ==============================================================================

run_quality_control_parallel() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided for parallel QC"; return 1; }
	
	# Export required variables and functions for parallel execution
	export PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE
	export RAW_DIR_ROOT TRIM_DIR_ROOT FASTQC_ROOT THREADS_PER_JOB
	# Export cached paths and error log for parallel worker logging
	export _CONDA_PROFILE_SCRIPT ERROR_WARN_FILE
	# Export _log_impl (core logger) alongside its callers — without it, log_info/log_warn/log_error
	# fail silently in GNU Parallel subshells because they delegate to _log_impl.
	export -f _log_impl timestamp log log_info log_warn log_error log_step rename_fastqc_outputs 2>/dev/null || true

	_qc_worker() {
		local SRR="$1"

		# Activate conda environment in subshell (uses cached path to avoid dirname subshell)
		if [[ -n "$CONDA_PREFIX" ]]; then
			source "${_CONDA_PROFILE_SCRIPT:-${CONDA_EXE%/*}/../etc/profile.d/conda.sh}" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi

		local RAW_DIR="$RAW_DIR_ROOT/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local srr_outdir="$FASTQC_ROOT/$SRR"
		mkdir -p "$srr_outdir"

		# QC for raw files (single folder per SRR)
		if [[ -d "$RAW_DIR" ]]; then
			local _raw_qc=("$srr_outdir"/*_raw_fastqc.html); if [[ ! -f "${_raw_qc[0]:-}" ]]; then
				log_info "Running FastQC on raw files for $SRR"
				fastqc -t "${THREADS_PER_JOB:-2}" -o "$srr_outdir" \
					"$RAW_DIR"/${SRR}*.fastq* 2>/dev/null || log_warn "FastQC failed for raw $SRR"
				rename_fastqc_outputs "$srr_outdir" "raw" "$RAW_DIR"/${SRR}*.fastq*
			fi
		fi

		# QC for trimmed files (same SRR folder, renamed as trimmed)
		if [[ -d "$TrimGalore_DIR" ]]; then
			local _trim_qc=("$srr_outdir"/*_trimmed_fastqc.html); if [[ ! -f "${_trim_qc[0]:-}" ]]; then
				log_info "Running FastQC on trimmed files for $SRR"
				fastqc -t "${THREADS_PER_JOB:-2}" -o "$srr_outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || log_warn "FastQC failed for trimmed $SRR"
				rename_fastqc_outputs "$srr_outdir" "trimmed" "$TrimGalore_DIR"/${SRR}*val*.fq*
			fi
		fi
	}
	export -f _qc_worker
	
	# Run FastQC in parallel for all samples
	parallel \
		--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
		--env RAW_DIR_ROOT --env TRIM_DIR_ROOT --env FASTQC_ROOT --env THREADS_PER_JOB \
		-j "${JOBS:-2}" \
		--halt soon,fail,1 \
		--joblog "$FASTQC_ROOT/parallel_fastqc.log" \
		_qc_worker {} \
		< <(printf "%s\n" "${SRR_LIST[@]}")
}

# ==============================================================================
# QC SUMMARY FUNCTIONS
# ==============================================================================

# generate_qc_summary() — removed (dead code; never called by active pipeline)
