#!/bin/bash
# ==============================================================================
# GENE EXPRESSION ANALYSIS (GEA) PIPELINE
# RNA-seq analysis pipeline using multiple alignment/quantification methods
# Author: Mark Cyril R. Mercado | Version: v12 | Date: December 2025
# ==============================================================================

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
	#"config/1_download_and_trim/HPC_download_and_trim.sh"		# Download + trim all SRRs

	# --- Test runs (all M1-M5, 3 SRRs) ---
	#"config/2_alignment/HPC_test_genome_M1_M3.sh"				# M1 + M3 (genome FASTA)
	#"config/2_alignment/HPC_test_transcript_M2_M4_M5.sh"		# M2 + M4 + M5 (transcript FASTA)

	# --- Test runs (individual methods) ---
	#"config/2_alignment/HPC_test_hisat2.sh"
	#"config/2_alignment/HPC_test_star.sh"
	#"config/2_alignment/HPC_test_salmon_bowtie2.sh"
	#"config/2_alignment/HPC_test_genome.sh"

	# --- Full runs ---
	"config/2_alignment/HPC_full_ref_guided.sh"				# Reference-guided (M1 + M3)
	"config/2_alignment/HPC_full_non_ref_guided.sh"			# Non-reference-guided (M2 + M4 + M5)

	# --- Local ---
	#"config/2_alignment/local_test.sh"							# Local testing
)

# ==============================================================================
# PIPELINE FLAG HELPER
# ==============================================================================

_has_stage() { printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^$1$" && echo "TRUE" || echo "FALSE"; }

set_pipeline_flags() {
	RUN_MAMBA_INSTALLATION=$(_has_stage "MAMBA_INSTALLATION")
	RUN_DOWNLOAD_SRR=$(_has_stage "DOWNLOAD_SRR")
	RUN_TRIM_SRR=$(_has_stage "TRIM_SRR")
	RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR=$(_has_stage "DOWNLOAD_TRIM_and_DELETE_RAW_SRR")
	RUN_GZIP_TRIMMED_FILES=$(_has_stage "GZIP_TRIMMED_FILES")
	RUN_DELETE_RAW_SRR=$(_has_stage "DELETE_RAW_SRR")
	RUN_QUALITY_CONTROL=$(_has_stage "QUALITY_CONTROL")
	RUN_METHOD_1_HISAT2_REF_GUIDED=$(_has_stage "METHOD_1_HISAT2_REF_GUIDED")
	RUN_METHOD_2_HISAT2_DE_NOVO=$(_has_stage "METHOD_2_HISAT2_DE_NOVO")
	RUN_METHOD_3_STAR_ALIGNMENT=$(_has_stage "METHOD_3_STAR_ALIGNMENT")
	RUN_METHOD_4_SALMON_SAF=$(_has_stage "METHOD_4_SALMON_SAF")
	RUN_METHOD_5_BOWTIE2_RSEM=$(_has_stage "METHOD_5_BOWTIE2_RSEM")
	RUN_DELETE_TRIMMED_FASTQ_FILES=$(_has_stage "DELETE_TRIMMED_FASTQ_FILES")
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
	catalog_all_software
	log_configuration
	log_step "Script started at: $(date -d @$start_time)"

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

	# --- Alignment Methods ---
	if [[ $RUN_METHOD_1_HISAT2_REF_GUIDED == "TRUE" ]]; then
		log_step "STEP 02a: HISAT2 Reference-Guided Pipeline"
		if [[ -z "${gtf_file:-}" || ! -f "${gtf_file:-}" ]]; then
			log_error "GTF file required for reference-guided alignment: ${gtf_file:-<unset>}"
			log_error "Skipping Method 1 — configure gtf_file variable"
			((method_failures++)) || true
		elif hisat2_ref_guided_pipeline --FASTA "$fasta" --GTF "$gtf_file" \
			${HISAT2_STRANDNESS:+--STRANDNESS "$HISAT2_STRANDNESS"} \
			--RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 1 completed successfully"
		else
			log_error "Method 1 failed (exit code: $?) — continuing"
			((method_failures++)) || true
		fi
	fi

	if [[ $RUN_METHOD_2_HISAT2_DE_NOVO == "TRUE" ]]; then
		log_step "STEP 02b: HISAT2 De Novo Pipeline"
		if hisat2_de_novo_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 2 completed successfully"
		else
			log_error "Method 2 failed (exit code: $?) — continuing"
			((method_failures++)) || true
		fi
	fi

	if [[ $RUN_METHOD_3_STAR_ALIGNMENT == "TRUE" ]]; then
		log_step "STEP 03: STAR Splice-Aware Alignment"
		if star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 3 completed successfully"
		else
			log_error "Method 3 failed (exit code: $?) — continuing"
			((method_failures++)) || true
		fi
	fi

	if [[ $RUN_METHOD_4_SALMON_SAF == "TRUE" ]]; then
		log_step "STEP 04: Salmon SAF Quantification"
		if [[ ! -f "${decoy:-}" ]]; then
			log_warn "Genome file '${decoy:-<unset>}' not found — skipping Salmon SAF pipeline."
		elif salmon_saf_pipeline --FASTA "$fasta" --GENOME "$decoy" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 4 completed successfully"
		else
			log_error "Method 4 failed (exit code: $?) — continuing"
			((method_failures++)) || true
		fi
	fi

	if [[ $RUN_METHOD_5_BOWTIE2_RSEM == "TRUE" ]]; then
		log_step "STEP 05: Bowtie2 + RSEM Quantification"
		if bowtie2_rsem_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 5 completed successfully"
		else
			log_error "Method 5 failed (exit code: $?) — continuing"
			((method_failures++)) || true
		fi
	fi

	compare_methods_summary "$fasta_tag"
	catalog_all_software

	end_time=$(date +%s)
	elapsed=$((end_time - start_time))
	log_step "Final timing"
	log_info "Script ended at: $(date -d @$end_time)"
	log_info "Elapsed time: $(date -u -d @${elapsed} +%H:%M:%S)"

	if [[ $method_failures -gt 0 ]]; then
		log_error "$method_failures method(s) failed for $fasta_tag"
	fi
	return $method_failures
}

# ==============================================================================
# EXECUTE
# ==============================================================================

[[ ${#CONFIG_FILES[@]} -eq 0 ]] && { echo "ERROR: No configuration files listed in CONFIG_FILES."; exit 1; }

# Cleanup trap: log summary on exit; clean up STAR temp dirs on signal kill
_pipeline_cleanup() {
	local rc=$?
	# Remove any orphan STAR temp directories left by interrupted runs
	find "${PROJECT_ROOT}" -maxdepth 4 -type d -name '_STARtmp*' -exec rm -rf {} + 2>/dev/null || true
	if [[ $rc -ne 0 ]]; then
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] Pipeline terminated with exit code $rc" >&2
	fi
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] Pipeline finished. See logs under: 1_SRRs/logs/ and 2_ALIGNMENT_RESULTs/logs/"
}
trap _pipeline_cleanup EXIT

total_failures=0

for config_file in "${CONFIG_FILES[@]}"; do
	echo ""
	echo "=============================================================================="
	echo "  LOADING CONFIGURATION: $config_file"
	echo "=============================================================================="

	[[ -f "$config_file" ]] || { echo "ERROR: Configuration file not found: $config_file"; exit 1; }
	GENOME_REF_PAIRS=()
	ALL_FASTA_FILES=()
	unset gtf_file STAR_TRANSCRIPTOME_FASTA decoy 2>/dev/null || true
	source "$config_file"
	set_pipeline_flags

	mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT"
	setup_logging
	switch_log_stage "1_SRRs"

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
			export gtf_file STAR_TRANSCRIPTOME_FASTA
			local _rc=0
			run_all --FASTA "$_fasta" --RNASEQ_LIST "${SRR_COMBINED_LIST[@]}" || _rc=$?
			total_failures=$((total_failures + _rc))
		done
		unset _pair _remainder _fasta gtf_file STAR_TRANSCRIPTOME_FASTA
	elif [[ ${#ALL_FASTA_FILES[@]} -gt 0 ]]; then
		for fasta_input in "${ALL_FASTA_FILES[@]}"; do
			local _rc=0
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

	echo ""
	echo "=============================================================================="
	echo "  FINISHED CONFIG: $config_file"
	echo "=============================================================================="
done

if [[ $total_failures -gt 0 ]]; then
	echo "PIPELINE COMPLETED WITH $total_failures METHOD FAILURE(S)"
	exit 1
fi
echo "END OF SCRIPT"
