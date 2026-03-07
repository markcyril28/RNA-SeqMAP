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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_preproc.sh"

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

		local base=$(basename "$input")
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

	if command -v fastqc >/dev/null 2>&1; then
		# QC for raw files (single folder per SRR)
		if [[ -d "$RAW_DIR" ]]; then
			if ! compgen -G "$srr_outdir/*_raw_fastqc.html" >/dev/null; then
				log_info "Running FastQC on raw files for $SRR"
				run_with_space_time_log fastqc -t "${THREADS:-2}" -o "$srr_outdir" \
					"$RAW_DIR"/${SRR}*.fastq* 2>/dev/null || log_warn "FastQC failed for raw $SRR"
				rename_fastqc_outputs "$srr_outdir" "raw" "$RAW_DIR"/${SRR}*.fastq*
			fi
		fi

		# QC for trimmed files (same SRR folder, renamed as trimmed)
		if [[ -d "$TrimGalore_DIR" ]]; then
			if ! compgen -G "$srr_outdir/*_trimmed_fastqc.html" >/dev/null; then
				log_info "Running FastQC on trimmed files for $SRR"
				run_with_space_time_log fastqc -t "${THREADS:-2}" -o "$srr_outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || log_warn "FastQC failed for trimmed $SRR"
				rename_fastqc_outputs "$srr_outdir" "trimmed" "$TrimGalore_DIR"/${SRR}*val*.fq*
			fi
		fi
	else
		log_warn "FastQC not found. Skipping QC for $SRR."
	fi

	# Run MultiQC to aggregate results
	run_multiqc
}

run_multiqc() {
	if command -v multiqc >/dev/null 2>&1; then
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
	export -f log_info log_warn log_error log_step run_with_space_time_log rename_fastqc_outputs 2>/dev/null || true
	
	_qc_worker() {
		local SRR="$1"
		
		# Activate conda environment in subshell
		if [[ -n "$CONDA_PREFIX" ]]; then
			source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local RAW_DIR="$RAW_DIR_ROOT/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local srr_outdir="$FASTQC_ROOT/$SRR"
		mkdir -p "$srr_outdir"
		
		# QC for raw files (single folder per SRR)
		if [[ -d "$RAW_DIR" ]]; then
			if ! compgen -G "$srr_outdir/*_raw_fastqc.html" >/dev/null; then
				echo "[INFO] Running FastQC on raw files for $SRR"
				fastqc -t "${THREADS_PER_JOB:-2}" -o "$srr_outdir" \
					"$RAW_DIR"/${SRR}*.fastq* 2>/dev/null || echo "[WARN] FastQC failed for raw $SRR"
				rename_fastqc_outputs "$srr_outdir" "raw" "$RAW_DIR"/${SRR}*.fastq*
			fi
		fi
		
		# QC for trimmed files (same SRR folder, renamed as trimmed)
		if [[ -d "$TrimGalore_DIR" ]]; then
			if ! compgen -G "$srr_outdir/*_trimmed_fastqc.html" >/dev/null; then
				echo "[INFO] Running FastQC on trimmed files for $SRR"
				fastqc -t "${THREADS_PER_JOB:-2}" -o "$srr_outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || echo "[WARN] FastQC failed for trimmed $SRR"
				rename_fastqc_outputs "$srr_outdir" "trimmed" "$TrimGalore_DIR"/${SRR}*val*.fq*
			fi
		fi
	}
	export -f _qc_worker
	
	# Run FastQC in parallel for all samples
	printf "%s\n" "${SRR_LIST[@]}" | parallel --env PATH --env CONDA_PREFIX -j "${JOBS:-2}" _qc_worker {}
}

# ==============================================================================
# QC SUMMARY FUNCTIONS
# ==============================================================================

generate_qc_summary() {
	local output_file="${1:-$FASTQC_ROOT/qc_summary.txt}"
	
	log_step "Generating QC Summary"
	
	{
		echo "=========================================="
		echo "FastQC Summary Report"
		echo "Generated: $(date)"
		echo "=========================================="
		echo ""
		
		# Find all fastqc_data.txt files and extract key metrics
		for summary_file in "$FASTQC_ROOT"/*/*_fastqc/summary.txt; do
			if [[ -f "$summary_file" ]]; then
				local sample_dir=$(dirname "$summary_file")
				local sample_name=$(basename "$sample_dir" | sed 's/_fastqc$//')
				echo "Sample: $sample_name"
				cat "$summary_file"
				echo ""
			fi
		done
	} > "$output_file"
	
	log_info "QC summary saved to: $output_file"
}
