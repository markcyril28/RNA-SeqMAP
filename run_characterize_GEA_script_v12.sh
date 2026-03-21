#!/bin/bash
# ==============================================================================
# GENE EXPRESSION ANALYSIS (GEA) PIPELINE
# RNA-seq analysis pipeline using multiple alignment/quantification methods
# Author: Mark Cyril R. Mercado | Version: v12 | Date: December 2025
# ==============================================================================

set -o pipefail   # -e/-u omitted intentionally (sourced functions use boolean returns)

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT" || exit 1

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================

# Execution mode:
#   skip      - Resume interrupted run; skip steps with existing outputs (default)
#   overwrite - Force clean rerun, overwriting all existing output files
OVERWRITE_MODE="${OVERWRITE_MODE:-skip}"
export OVERWRITE_MODE

# Active configuration file — uncomment as needed:
CONFIG_FILES=(
	# --- Download & Trim ---
	#"config/1_download_and_trim/HPC_download_and_trim.toml"	# Download + trim all SRRs

	# --- Test runs (all M1-M5, 3 SRRs) ---
	#"config/2_alignment/HPC_test_genome_M1_M3.toml"			# M1 + M3 (genome FASTA)
	#"config/2_alignment/HPC_test_transcript_M2_M4_M5.toml"	# M2 + M4 + M5 (transcript FASTA)

	# --- Full runs ---
	"config/2_alignment/HPC_full_ref_guided.toml"				# Reference-guided (M1 + M3)
	"config/2_alignment/HPC_full_non_ref_guided.toml"			# Non-reference-guided (M2 + M4 + M5)

	# --- Local ---
	#"config/2_alignment/local_full_ref_guided.toml"			# Local ref-guided (M1 + M3)
	#"config/2_alignment/local_full_non_ref_guided.toml"		# Local non-ref-guided (M2 + M4 + M5)
)

# ==============================================================================
# PIPELINE FLAG HELPER
# ==============================================================================

# Build associative array from PIPELINE_STAGES for O(1) lookup — zero subshell spawns.
set_pipeline_flags() {
	# Build stage lookup set (single pass)
	local -A _stage_set=()
	local s
	for s in "${PIPELINE_STAGES[@]}"; do _stage_set["$s"]=1; done

	RUN_MAMBA_INSTALLATION="${_stage_set[MAMBA_INSTALLATION]:+TRUE}"
	RUN_MAMBA_INSTALLATION="${RUN_MAMBA_INSTALLATION:-FALSE}"
	RUN_DOWNLOAD_SRR="${_stage_set[DOWNLOAD_SRR]:+TRUE}"
	RUN_DOWNLOAD_SRR="${RUN_DOWNLOAD_SRR:-FALSE}"
	RUN_TRIM_SRR="${_stage_set[TRIM_SRR]:+TRUE}"
	RUN_TRIM_SRR="${RUN_TRIM_SRR:-FALSE}"
	RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR="${_stage_set[DOWNLOAD_TRIM_and_DELETE_RAW_SRR]:+TRUE}"
	RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR="${RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR:-FALSE}"
	RUN_GZIP_TRIMMED_FILES="${_stage_set[GZIP_TRIMMED_FILES]:+TRUE}"
	RUN_GZIP_TRIMMED_FILES="${RUN_GZIP_TRIMMED_FILES:-FALSE}"
	RUN_DELETE_RAW_SRR="${_stage_set[DELETE_RAW_SRR]:+TRUE}"
	RUN_DELETE_RAW_SRR="${RUN_DELETE_RAW_SRR:-FALSE}"
	RUN_QUALITY_CONTROL="${_stage_set[QUALITY_CONTROL]:+TRUE}"
	RUN_QUALITY_CONTROL="${RUN_QUALITY_CONTROL:-FALSE}"
	RUN_METHOD_1_HISAT2_REF_GUIDED="${_stage_set[METHOD_1_HISAT2_REF_GUIDED]:+TRUE}"
	RUN_METHOD_1_HISAT2_REF_GUIDED="${RUN_METHOD_1_HISAT2_REF_GUIDED:-FALSE}"
	RUN_METHOD_2_HISAT2_DE_NOVO="${_stage_set[METHOD_2_HISAT2_DE_NOVO]:+TRUE}"
	RUN_METHOD_2_HISAT2_DE_NOVO="${RUN_METHOD_2_HISAT2_DE_NOVO:-FALSE}"
	RUN_METHOD_3_STAR_ALIGNMENT="${_stage_set[METHOD_3_STAR_ALIGNMENT]:+TRUE}"
	RUN_METHOD_3_STAR_ALIGNMENT="${RUN_METHOD_3_STAR_ALIGNMENT:-FALSE}"
	RUN_METHOD_4_SALMON_SAF="${_stage_set[METHOD_4_SALMON_SAF]:+TRUE}"
	RUN_METHOD_4_SALMON_SAF="${RUN_METHOD_4_SALMON_SAF:-FALSE}"
	RUN_METHOD_5_BOWTIE2_RSEM="${_stage_set[METHOD_5_BOWTIE2_RSEM]:+TRUE}"
	RUN_METHOD_5_BOWTIE2_RSEM="${RUN_METHOD_5_BOWTIE2_RSEM:-FALSE}"
	RUN_DELETE_TRIMMED_FASTQ_FILES="${_stage_set[DELETE_TRIMMED_FASTQ_FILES]:+TRUE}"
	RUN_DELETE_TRIMMED_FASTQ_FILES="${RUN_DELETE_TRIMMED_FASTQ_FILES:-FALSE}"
}

# ==============================================================================
# MAIN PIPELINE FUNCTION
# ==============================================================================

run_all() {
	local fasta="" rnaseq_list=()
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA)      fasta="$2"; shift 2 ;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done ;;
			*) shift ;;
		esac
	done

	local start_time end_time elapsed
	start_time=$(date +%s)

	local fasta_base fasta_tag
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	set_fasta_output_dirs "$fasta_tag"

	setup_logging
	switch_log_stage "1_SRRs"
	# Catalog software once per pipeline execution, not per FASTA
	if [[ "${_SOFTWARE_CATALOGED:-}" != "true" ]]; then
		catalog_all_software
		_SOFTWARE_CATALOGED="true"
	fi
	log_configuration
	log_step "Script started at: $(date -d "@$start_time")"

	log_info "SRR samples to process:"
	for srr in "${rnaseq_list[@]}"; do log_info "$srr"; done

	# Track method failures for aggregate exit code
	local method_failures=0

	# --- Preprocessing ---
	if [[ $RUN_DOWNLOAD_SRR == "TRUE" ]]; then
		if [[ $RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR == "TRUE" ]]; then
			log_warn "DOWNLOAD_SRR skipped: DOWNLOAD_TRIM_and_DELETE_RAW_SRR is enabled"
		else
			log_step "STEP 01a: Download RNA-seq data"
			download_srrs_parallel "${rnaseq_list[@]}"
		fi
	fi

	if [[ $RUN_TRIM_SRR == "TRUE" ]]; then
		if [[ $RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR == "TRUE" ]]; then
			log_warn "TRIM_SRR skipped: DOWNLOAD_TRIM_and_DELETE_RAW_SRR is enabled"
		else
			log_step "STEP 01b: Trim RNA-seq data"
			trim_srrs_trimmomatic_parallel "${rnaseq_list[@]}"
		fi
	fi

	if [[ $RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR == "TRUE" ]]; then
		log_step "STEP 01ab: Download, Trim, and Delete Raw SRR data"
		export DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING="TRUE"
		download_and_trim_srrs_parallel "${rnaseq_list[@]}"
	fi

	if [[ $RUN_DELETE_RAW_SRR == "TRUE" ]]; then
		log_step "STEP 01d: Delete Raw SRR files"
		delete_raw_srr_by_srr_list "${rnaseq_list[@]}"
	fi

	if [[ $RUN_QUALITY_CONTROL == "TRUE" ]]; then
		log_step "STEP 01c: Quality Control analysis"
		run_quality_control_all "${rnaseq_list[@]}"
	fi

	switch_log_stage "2_ALIGNMENT_RESULTs"

	# --- Alignment Methods (parallel when independent) ---
	# M1-M5 produce output in isolated directories and do not depend on each other.
	# When PARALLEL_METHODS is enabled and multiple methods are requested, dispatch
	# them concurrently as background jobs to reduce total wall-clock time.
	local _enabled_methods=()
	local _method_cmds=()

	if [[ $RUN_METHOD_1_HISAT2_REF_GUIDED == "TRUE" ]]; then
		if [[ -z "${gtf_file:-}" || ! -f "${gtf_file:-}" ]]; then
			log_error "GTF file required for reference-guided alignment: ${gtf_file:-<unset>}"
			log_error "Skipping Method 1 — configure gtf_file variable"
			((method_failures++)) || true
		else
			_enabled_methods+=("M1")
			_method_cmds+=("hisat2_ref_guided_pipeline --FASTA \"$fasta\" --GTF \"$gtf_file\" ${HISAT2_STRANDNESS:+--STRANDNESS \"$HISAT2_STRANDNESS\"} --RNASEQ_LIST ${rnaseq_list[*]}")
		fi
	fi
	if [[ $RUN_METHOD_2_HISAT2_DE_NOVO == "TRUE" ]]; then
		_enabled_methods+=("M2")
		_method_cmds+=("hisat2_de_novo_pipeline --FASTA \"$fasta\" --RNASEQ_LIST ${rnaseq_list[*]}")
	fi
	if [[ $RUN_METHOD_3_STAR_ALIGNMENT == "TRUE" ]]; then
		_enabled_methods+=("M3")
		_method_cmds+=("star_alignment_pipeline --FASTA \"$fasta\" --RNASEQ_LIST ${rnaseq_list[*]}")
	fi
	if [[ $RUN_METHOD_4_SALMON_SAF == "TRUE" ]]; then
		if [[ ! -f "${decoy:-}" ]]; then
			log_error "Genome file '${decoy:-<unset>}' not found — skipping Salmon SAF pipeline."
			((method_failures++)) || true
		else
			_enabled_methods+=("M4")
			_method_cmds+=("salmon_saf_pipeline --FASTA \"$fasta\" --GENOME \"$decoy\" --RNASEQ_LIST ${rnaseq_list[*]}")
		fi
	fi
	if [[ $RUN_METHOD_5_BOWTIE2_RSEM == "TRUE" ]]; then
		_enabled_methods+=("M5")
		_method_cmds+=("bowtie2_rsem_pipeline --FASTA \"$fasta\" --RNASEQ_LIST ${rnaseq_list[*]}")
	fi

	if [[ ${#_enabled_methods[@]} -eq 0 ]]; then
		log_info "No alignment methods enabled"
	elif [[ ${#_enabled_methods[@]} -eq 1 ]]; then
		# Single method — run directly (no overhead from background dispatch)
		log_step "Running ${_enabled_methods[0]} (single method)"
		if eval "${_method_cmds[0]}"; then
			log_info "${_enabled_methods[0]} completed successfully"
		else
			log_error "${_enabled_methods[0]} failed (exit code: $?) — continuing"
			((method_failures++)) || true
		fi
	elif [[ "${PARALLEL_METHODS:-TRUE}" == "TRUE" ]] && (( ${#_enabled_methods[@]} > 1 )); then
		# Multiple methods — run concurrently as background jobs.
		# Each method already manages its own GNU Parallel pool for per-sample work,
		# so cross-method parallelism adds no contention beyond shared I/O bandwidth.
		log_step "Running ${#_enabled_methods[@]} methods in parallel: ${_enabled_methods[*]}"
		local -a _method_pids=()
		local -a _method_logs=()
		for _i in "${!_enabled_methods[@]}"; do
			local _m="${_enabled_methods[$_i]}"
			local _mlog="${LOG_DIR}/method_${_m}_${fasta_tag}.log"
			_method_logs+=("$_mlog")
			log_info "Dispatching ${_m} → $_mlog"
			eval "${_method_cmds[$_i]}" > "$_mlog" 2>&1 &
			_method_pids+=($!)
		done

		# Wait for all methods and collect exit codes
		for _i in "${!_enabled_methods[@]}"; do
			local _m="${_enabled_methods[$_i]}"
			local _pid="${_method_pids[$_i]}"
			if wait "$_pid"; then
				log_info "${_m} completed successfully (pid=$_pid)"
			else
				log_error "${_m} failed (exit code: $?, pid=$_pid) — see ${_method_logs[$_i]}"
				((method_failures++)) || true
			fi
		done
		# Batch-append all method logs in a single I/O operation (avoids repeated open/seek/close)
		local _existing_logs=()
		for _mlog in "${_method_logs[@]}"; do
			[[ -f "$_mlog" ]] && _existing_logs+=("$_mlog")
		done
		[[ ${#_existing_logs[@]} -gt 0 ]] && cat "${_existing_logs[@]}" >> "$LOG_FILE"
	else
		# PARALLEL_METHODS=FALSE: sequential fallback
		log_step "Running ${#_enabled_methods[@]} methods sequentially"
		for _i in "${!_enabled_methods[@]}"; do
			local _m="${_enabled_methods[$_i]}"
			log_step "Running ${_m}"
			if eval "${_method_cmds[$_i]}"; then
				log_info "${_m} completed successfully"
			else
				log_error "${_m} failed (exit code: $?) — continuing"
				((method_failures++)) || true
			fi
		done
	fi

	compare_methods_summary "$fasta_tag"

	end_time=$(date +%s)
	elapsed=$((end_time - start_time))
	log_step "Final timing"
	log_info "Script ended at: $(date -d "@$end_time")"
	# Pure bash arithmetic (avoids date subshell spawn)
	log_info "Elapsed time: $(printf '%02d:%02d:%02d' $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60)))"

	if [[ $method_failures -gt 0 ]]; then
		log_error "$method_failures method(s) failed for $fasta_tag"
	fi
	return $method_failures
}

# ==============================================================================
# EXECUTE
# ==============================================================================

[[ ${#CONFIG_FILES[@]} -eq 0 ]] && { echo "[ERROR] No configuration files listed in CONFIG_FILES." >&2; exit 1; }
# NOTE: log_error is not yet available here — logging module is sourced below.

# Source logging early so log_step/log_error/log_info are available before configs.
# The modules_loader double-source guard ensures this is safe when configs re-source it.
source "${PROJECT_ROOT}/modules/logging/logging_utils.sh" 2>/dev/null || {
	# Minimal fallback if logging module cannot be loaded
	log_info()  { echo "[INFO]  $*"; }
	log_warn()  { echo "[WARN]  $*"; }
	log_error() { echo "[ERROR] $*" >&2; }
	log_step()  { echo ""; echo "==> $*"; }
}

# Source TOML parser and shared runtime defaults
source "${PROJECT_ROOT}/config/shared/toml_parser.sh"
source "${PROJECT_ROOT}/config/shared/runtime_defaults.sh"

# Cleanup trap: log summary on exit; clean up STAR temp dirs on signal kill
_pipeline_cleanup() {
	local rc=$?
	# Only search for orphan STAR temp dirs when STAR was actually used (avoids
	# traversing entire PROJECT_ROOT on every exit — saves ~0.5-2s on large trees)
	if [[ "${RUN_METHOD_3_STAR_ALIGNMENT:-}" == "TRUE" ]]; then
		find "${PROJECT_ROOT}/2_ALIGNMENT_RESULTs" -maxdepth 4 -type d -name '_STARtmp*' -exec rm -rf {} + 2>/dev/null || true
	fi
	if [[ $rc -ne 0 ]]; then
		log_error "Pipeline terminated with exit code $rc"
	fi
	log_info "Pipeline finished. See logs under: 1_SRRs/logs/ and 2_ALIGNMENT_RESULTs/logs/"
}
trap _pipeline_cleanup EXIT

# Big O: O(C × R × M × S) where C=configs, R=ref_pairs, M=enabled_methods, S=samples.
# Methods run in parallel when PARALLEL_METHODS=TRUE, reducing M dimension to O(1) wall-clock.
# Per-method inner loops are parallelized via GNU Parallel (S/JOBS threads).
total_failures=0

for config_file in "${CONFIG_FILES[@]}"; do
	log_step "LOADING CONFIGURATION: $config_file"

	[[ -f "$config_file" ]] || { log_error "Configuration file not found: $config_file"; exit 1; }
	GENOME_REF_PAIRS=()
	ALL_FASTA_FILES=()
	PIPELINE_STAGES=()
	SRR_COMBINED_LIST=()
	unset gtf_file STAR_TRANSCRIPTOME_FASTA decoy DECOY KEEP_BAM_GLOBAL STAR_READ_LENGTH HISAT2_STRANDNESS
	load_toml "$config_file" || { log_error "Failed to load config: $config_file"; exit 1; }
	# Export STAR_READ_LENGTH if set by config
	[[ -n "${STAR_READ_LENGTH:-}" ]] && export STAR_READ_LENGTH

	# Load SRR datasets: test configs define SRR_COMBINED_LIST inline;
	# full configs load from shared srr_datasets.toml
	if [[ ${#SRR_COMBINED_LIST[@]} -eq 0 ]]; then
		load_toml_srr_datasets "${PROJECT_ROOT}/config/shared/srr_datasets.toml"
	fi

	# Map TOML uppercase keys to lowercase aliases used by method modules
	[[ -n "${DECOY:-}" ]] && decoy="$DECOY"
	[[ -n "${KEEP_BAM_GLOBAL:-}" ]] && keep_bam_global="$KEEP_BAM_GLOBAL"
	[[ -n "${HISAT2_STRANDNESS:-}" ]] || HISAT2_STRANDNESS=""
	[[ -n "${POST_PROCESSING_ROOT:-}" ]] && export POST_PROCESSING_ROOT

	# Re-derive THREADS_PER_JOB from potentially updated THREADS/JOBS
	if [[ "${USE_GNU_PARALLEL:-FALSE}" == "TRUE" ]]; then
		THREADS_PER_JOB=$((${THREADS:-4} / ${JOBS:-1}))
		[[ $THREADS_PER_JOB -lt 1 ]] && THREADS_PER_JOB=1
	else
		THREADS_PER_JOB="${THREADS:-4}"
	fi
	export THREADS JOBS USE_GNU_PARALLEL THREADS_PER_JOB keep_bam_global

	set_pipeline_flags

	mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT"
	# setup_logging is called inside run_all(); avoid redundant re-initialization per config.
	# Only switch log stage if logging is already initialized (first config handled by run_all).
	[[ "${LOGGING_INITIALIZED:-}" == "true" ]] && switch_log_stage "1_SRRs"

	[[ $RUN_MAMBA_INSTALLATION == "TRUE" ]] && mamba_install
	if [[ $RUN_GZIP_TRIMMED_FILES == "TRUE" ]]; then
		log_step "Gzipping trimmed FASTQ files"
		gzip_trimmed_fastq_files
	fi

	if [[ ${#GENOME_REF_PAIRS[@]} -gt 0 ]]; then
		for _pair in "${GENOME_REF_PAIRS[@]}"; do
			gtf_file="${_pair%%|*}"
			_remainder="${_pair#*|}"
			_fasta="${_remainder%%|*}"
			STAR_TRANSCRIPTOME_FASTA="${_remainder#*|}"
			# If no third field existed, strip fell back to the full string — reset to empty
			[[ "$STAR_TRANSCRIPTOME_FASTA" == "$_fasta" ]] && STAR_TRANSCRIPTOME_FASTA=""
			export gtf_file STAR_TRANSCRIPTOME_FASTA
			_rc=0
			run_all --FASTA "$_fasta" --RNASEQ_LIST "${SRR_COMBINED_LIST[@]}" || _rc=$?
			total_failures=$((total_failures + _rc))
		done
		unset _pair _remainder _fasta gtf_file STAR_TRANSCRIPTOME_FASTA
	elif [[ ${#ALL_FASTA_FILES[@]} -gt 0 ]]; then
		for fasta_input in "${ALL_FASTA_FILES[@]}"; do
			_rc=0
			run_all --FASTA "$fasta_input" --RNASEQ_LIST "${SRR_COMBINED_LIST[@]}" || _rc=$?
			total_failures=$((total_failures + _rc))
		done
	else
		log_warn "No GENOME_REF_PAIRS or ALL_FASTA_FILES defined in $config_file — skipping alignment"
	fi

	# --- Cleanup ---
	switch_log_stage "1_SRRs"
	if [[ $RUN_DELETE_TRIMMED_FASTQ_FILES == "TRUE" ]]; then
		log_step "Deleting trimmed FASTQ files"
		delete_trimmed_fastq_by_srr_list "${SRR_COMBINED_LIST[@]}"
	fi

	log_step "FINISHED CONFIG: $config_file"
done

if [[ $total_failures -gt 0 ]]; then
	log_error "PIPELINE COMPLETED WITH $total_failures METHOD FAILURE(S)"
	exit 1
fi
log_info "Pipeline execution completed"
