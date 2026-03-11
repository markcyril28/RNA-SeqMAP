#!/bin/bash
# ==============================================================================
# METHOD 5: BOWTIE2 + RSEM QUANTIFICATION PIPELINE
# ==============================================================================
# Quantify expression using Bowtie2 alignment + RSEM
# Reviewer-preferred method for publication
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M5_RSEM_SOURCED:-}" == "true" ]] && return 0
export M5_RSEM_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# RSEM CONFIGURATION - IMPORTANT PARAMETERS (tweak here)
# ==============================================================================

# Number of samples to process in parallel (requires USE_GNU_PARALLEL=TRUE)
MAX_PARALLEL_SAMPLES="${MAX_PARALLEL_SAMPLES:-2}"

# Threads allocated per RSEM job (auto-calculated from THREADS / MAX_PARALLEL_SAMPLES)
THREADS_PER_RSEM_JOB="${THREADS_PER_RSEM_JOB:-$((THREADS / MAX_PARALLEL_SAMPLES))}"
[[ $THREADS_PER_RSEM_JOB -lt 1 ]] && THREADS_PER_RSEM_JOB=1

# Library strandedness: none (unstranded), forward (sense), reverse (antisense/dUTP)
# Set to "reverse" for dUTP-based stranded libraries (most modern Illumina RNA-seq)
# Set to "auto" to auto-detect using Salmon --libType A (recommended)
RSEM_STRANDEDNESS="${RSEM_STRANDEDNESS:-auto}"

# Random seed for RSEM's EM algorithm (reproducibility across runs)
RSEM_SEED="${RSEM_SEED:-42}"

# ==============================================================================
# STRANDEDNESS AUTO-DETECTION (via Salmon --libType A)
# ==============================================================================

# Detect library strandedness by running Salmon on a small read subset.
# Caches the result per FASTA tag in the index directory.
# Usage: _rsem_detect_strandedness <fasta> <first_srr> <index_root>
# Sets: RSEM_STRANDEDNESS (global)
_rsem_detect_strandedness() {
	local fasta="$1"
	local first_srr="$2"
	local index_root="$3"
	local cache_file="$index_root/.detected_strandedness"

	# Return cached result if available
	if [[ -f "$cache_file" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		RSEM_STRANDEDNESS=$(cat "$cache_file")
		log_info "[STRANDEDNESS] Using cached result: $RSEM_STRANDEDNESS (from $cache_file)"
		return 0
	fi

	# Check that Salmon is available
	if ! command -v salmon >/dev/null 2>&1; then
		log_warn "[STRANDEDNESS] Salmon not found — falling back to 'none' (unstranded)"
		RSEM_STRANDEDNESS="none"
		return 0
	fi

	# Locate trimmed reads for the first sample
	find_trimmed_fastq "$first_srr"
	if [[ -z "$trimmed1" ]]; then
		log_warn "[STRANDEDNESS] No trimmed reads for $first_srr — falling back to 'none'"
		RSEM_STRANDEDNESS="none"
		return 0
	fi

	log_step "[STRANDEDNESS] Auto-detecting library strandedness using Salmon (sample: $first_srr)"

	local tmp_dir
	tmp_dir=$(mktemp -d "${index_root}/strandedness_detect_XXXXXX")
	local salmon_idx="$tmp_dir/salmon_idx"
	local salmon_quant="$tmp_dir/salmon_quant"

	# Build a lightweight Salmon index (using the same transcriptome FASTA)
	salmon index -t "$fasta" -i "$salmon_idx" --threads 4 -k 23 2>"$tmp_dir/salmon_index.log" || {
		log_warn "[STRANDEDNESS] Salmon index failed — falling back to 'none'"
		rm -rf "$tmp_dir"
		RSEM_STRANDEDNESS="none"
		return 0
	}

	# Quantify with --libType A (auto-detect) using a read subset (--numBootstraps 0, --validateMappings off for speed)
	local salmon_exit
	if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
		salmon quant -i "$salmon_idx" -l A \
			-1 "$trimmed1" -2 "$trimmed2" \
			-o "$salmon_quant" --threads 4 \
			--skipQuant 2>"$tmp_dir/salmon_quant.log"
		salmon_exit=$?
	else
		salmon quant -i "$salmon_idx" -l A \
			-r "$trimmed1" \
			-o "$salmon_quant" --threads 4 \
			--skipQuant 2>"$tmp_dir/salmon_quant.log"
		salmon_exit=$?
	fi

	if [[ $salmon_exit -ne 0 ]]; then
		log_warn "[STRANDEDNESS] Salmon quant failed — falling back to 'none'"
		rm -rf "$tmp_dir"
		RSEM_STRANDEDNESS="none"
		return 0
	fi

	# Parse the inferred library type from lib_format_counts.json
	local lib_format="$salmon_quant/lib_format_counts.json"
	if [[ ! -f "$lib_format" ]]; then
		log_warn "[STRANDEDNESS] lib_format_counts.json not found — falling back to 'none'"
		rm -rf "$tmp_dir"
		RSEM_STRANDEDNESS="none"
		return 0
	fi

	# Extract the expected_format field (Salmon's inferred library type)
	local inferred_type
	inferred_type=$(python3 -c "
import json, sys
with open('$lib_format') as f:
    data = json.load(f)
print(data.get('expected_format', 'U'))
" 2>/dev/null)

	# Map Salmon library type codes to RSEM strandedness
	# Salmon paired-end: IU=unstranded, ISF=forward(sense), ISR=reverse(antisense)
	# Salmon single-end: U=unstranded, SF=forward, SR=reverse
	case "$inferred_type" in
		IU|U)   RSEM_STRANDEDNESS="none" ;;
		ISF|SF) RSEM_STRANDEDNESS="forward" ;;
		ISR|SR) RSEM_STRANDEDNESS="reverse" ;;
		*)
			log_warn "[STRANDEDNESS] Unrecognized Salmon library type '$inferred_type' — falling back to 'none'"
			RSEM_STRANDEDNESS="none"
			;;
	esac

	# Also log the read counts per orientation for transparency
	local num_compat
	num_compat=$(python3 -c "
import json
with open('$lib_format') as f:
    d = json.load(f)
for k in ['compatible_fragment_ratio', 'num_compatible_fragments',
           'num_assigned_fragments']:
    if k in d: print(f'  {k}: {d[k]}')
for k in sorted(d.keys()):
    if 'strand' in k.lower() or k.startswith('read'): print(f'  {k}: {d[k]}')
" 2>/dev/null)
	[[ -n "$num_compat" ]] && log_info "[STRANDEDNESS] Salmon fragment stats:\n$num_compat"

	log_info "[STRANDEDNESS] Detected: $RSEM_STRANDEDNESS (Salmon inferred: $inferred_type)"

	# Cache the result
	mkdir -p "$index_root"
	echo "$RSEM_STRANDEDNESS" > "$cache_file"

	# Cleanup temp files
	rm -rf "$tmp_dir"

	return 0
}

# ==============================================================================
# MAIN PIPELINE: BOWTIE2 + RSEM
# ==============================================================================

bowtie2_rsem_pipeline() {
	local fasta="" rnaseq_list=()
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do rnaseq_list+=("$1"); shift; done;;
			*) log_error "Unknown arg: $1"; return 1;;
		esac
	done

	[[ -z "$fasta" ]] && { log_error "Usage: --FASTA genes.fa"; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")

	# Convert line endings if dos2unix is available
	command -v dos2unix >/dev/null 2>&1 && dos2unix "$fasta" 2>/dev/null || true

	local tag="$(basename "${fasta%.*}")"
	set_fasta_output_dirs "$tag"
	local rsem_idx="$RSEM_INDEX_ROOT/rsem_ref"
	local quant_root="$RSEM_QUANT_ROOT"
	local matrix_dir="$RSEM_MATRIX_ROOT"

	mkdir -p "$RSEM_INDEX_ROOT" "$quant_root" "$matrix_dir"

	# AUTO-DETECT STRANDEDNESS (runs once, caches result)
	if [[ "$RSEM_STRANDEDNESS" == "auto" ]]; then
		_rsem_detect_strandedness "$fasta" "${rnaseq_list[0]}" "$RSEM_INDEX_ROOT"
	fi

	# BUILD RSEM REFERENCE
	# Check all required index files: .grp (RSEM) + .rev.2.bt2 (written last by Bowtie2 build)
	if [[ -f "${rsem_idx}.grp" && -f "${rsem_idx}.rev.2.bt2" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[RSEM INDEX] RSEM reference already exists. Skipping."
	else
		log_step "Building RSEM reference for $tag"
		log_file_size "$fasta" "Input FASTA for RSEM index - $tag"
		run_with_space_time_log --input "$fasta" --output "$RSEM_INDEX_ROOT" \
			rsem-prepare-reference --bowtie2 "$fasta" "$rsem_idx"
		log_file_size "$RSEM_INDEX_ROOT" "RSEM index output - $tag"
	fi

	# QUANTIFY SAMPLES - parallel or sequential
	if _rsem_should_use_parallel && [[ ${#rnaseq_list[@]} -gt 1 ]]; then
		_rsem_quantify_parallel "$rsem_idx" "$quant_root" rnaseq_list[@] || \
			log_error "[RSEM] Some parallel quantification jobs failed — check logs before relying on matrix output"
	else
		_rsem_quantify_sequential "$rsem_idx" "$quant_root" rnaseq_list[@]
	fi

	# GENERATE MATRICES
	_create_rsem_matrices "$fasta" "$tag" "$quant_root" "$matrix_dir" rnaseq_list[@] || \
		log_warn "[RSEM] Matrix generation failed for $tag"

	log_step "COMPLETED: Bowtie2-RSEM pipeline for $tag"
}

# ==============================================================================
# QUANTIFICATION: SEQUENTIAL AND PARALLEL
# ==============================================================================

# Sequential sample quantification (default)
_rsem_quantify_sequential() {
	local rsem_idx="$1"
	local quant_root="$2"
	local arr_name="${3:-}"

	# Expand array from indirect reference
	local samples=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && samples+=("$s")
		done
	fi

	log_info "[RSEM QUANT] Running SEQUENTIAL quantification for ${#samples[@]} samples ($THREADS threads)"

	local failed_samples=0
	for SRR in "${samples[@]}"; do
		_rsem_process_single_sample "$SRR" "$rsem_idx" "$quant_root" "$THREADS" || {
			log_warn "[RSEM QUANT] Sample $SRR failed — continuing with remaining samples"
			failed_samples=$((failed_samples + 1))
		}
	done
	[[ $failed_samples -gt 0 ]] && log_warn "[RSEM QUANT] $failed_samples/${#samples[@]} sample(s) failed"
}

# ==============================================================================
# PARALLEL WORKER (top-level for reliable export -f across bash versions)
# ==============================================================================
# All required variables must be exported by _rsem_quantify_parallel before
# GNU Parallel dispatches this function:
#   rsem_idx, quant_root, threads_per_job, abs_trim_dir_root,
#   abs_error_warn_file, keep_bam_global, OVERWRITE_MODE,
#   BOWTIE2_MODE, RSEM_STRANDEDNESS
_rsem_parallel_worker() {
	local SRR="$1"

	_plog() {
		local level="$1"; shift
		local ts
		ts=$(date '+%Y-%m-%d %H:%M:%S')
		echo "[$ts] [$level] [RSEM-$SRR] $*"
		[[ "$level" != "INFO" && -n "${abs_error_warn_file:-}" ]] && echo "[$ts] [$level] [RSEM-$SRR] $*" >> "$abs_error_warn_file"
	}

	[[ -z "$SRR" ]] && { _plog "ERROR" "Empty SRR ID - skipping"; return 1; }

	# Reactivate conda in subshell if needed
	if [[ -n "${CONDA_PREFIX:-}" && -n "${CONDA_EXE:-}" ]]; then
		source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
		conda activate "${CONDA_DEFAULT_ENV:-base}" 2>/dev/null || true
	fi

	local out_dir="$quant_root/$SRR"
	mkdir -p "$out_dir"

	# Skip if already processed
	if [[ -f "$out_dir/${SRR}.genes.results" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		_plog "INFO" "Results already exist. Skipping."
		return 0
	fi

	# Locate trimmed reads (inlined to avoid function-export issues across bash versions)
	local trim_dir="$abs_trim_dir_root/$SRR"
	local trimmed1="" trimmed2=""

	if [[ -f "$trim_dir/${SRR}_1_val_1.fq.gz" && -f "$trim_dir/${SRR}_2_val_2.fq.gz" ]]; then
		trimmed1="$trim_dir/${SRR}_1_val_1.fq.gz"
		trimmed2="$trim_dir/${SRR}_2_val_2.fq.gz"
	elif [[ -f "$trim_dir/${SRR}_1_val_1.fq" && -f "$trim_dir/${SRR}_2_val_2.fq" ]]; then
		trimmed1="$trim_dir/${SRR}_1_val_1.fq"
		trimmed2="$trim_dir/${SRR}_2_val_2.fq"
	elif [[ -f "$trim_dir/${SRR}_trimmed.fq.gz" ]]; then
		trimmed1="$trim_dir/${SRR}_trimmed.fq.gz"
	elif [[ -f "$trim_dir/${SRR}_trimmed.fq" ]]; then
		trimmed1="$trim_dir/${SRR}_trimmed.fq"
	else
		for f in "$trim_dir"/${SRR}*val_1*.fq* "$trim_dir"/${SRR}*val_1*.gz; do
			[[ -f "$f" ]] && { trimmed1="$f"; break; }
		done
		for f in "$trim_dir"/${SRR}*val_2*.fq* "$trim_dir"/${SRR}*val_2*.gz; do
			[[ -f "$f" ]] && { trimmed2="$f"; break; }
		done
	fi

	if [[ -z "$trimmed1" ]]; then
		_plog "ERROR" "Missing trimmed reads in: $trim_dir"
		_plog "ERROR" "Contents: $(ls -la "$trim_dir" 2>&1 || echo 'Directory does not exist')"
		return 1
	fi

	if [[ ! -f "${rsem_idx}.grp" || ! -f "${rsem_idx}.rev.2.bt2" ]]; then
		_plog "ERROR" "RSEM index incomplete or not found: ${rsem_idx}.grp / .rev.2.bt2"
		return 1
	fi

	_plog "INFO" "Processing with $threads_per_job threads (strandedness: ${RSEM_STRANDEDNESS:-none})"

	local rsem_log="$out_dir/${SRR}.rsem.log"
	local rsem_exit_code
	if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
		rsem-calculate-expression \
			--paired-end \
			--bowtie2 \
			--bowtie2-sensitivity-level "${BOWTIE2_MODE:-sensitive}" \
			--strandedness "${RSEM_STRANDEDNESS:-none}" \
			--seed "$RSEM_SEED" \
			--num-threads "$threads_per_job" \
			"$trimmed1" "$trimmed2" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$rsem_log"
		rsem_exit_code=${PIPESTATUS[0]}
	else
		rsem-calculate-expression \
			--bowtie2 \
			--bowtie2-sensitivity-level "${BOWTIE2_MODE:-sensitive}" \
			--strandedness "${RSEM_STRANDEDNESS:-none}" \
			--seed "$RSEM_SEED" \
			--num-threads "$threads_per_job" \
			"$trimmed1" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$rsem_log"
		rsem_exit_code=${PIPESTATUS[0]}
	fi

	if [[ $rsem_exit_code -ne 0 ]]; then
		_plog "ERROR" "RSEM failed (exit code: $rsem_exit_code)"
		[[ -f "$rsem_log" ]] && tail -20 "$rsem_log" 2>/dev/null | \
			while IFS= read -r line; do _plog "ERROR" "  $line"; done
		return $rsem_exit_code
	fi

	if [[ ! -f "$out_dir/${SRR}.genes.results" ]]; then
		_plog "ERROR" "RSEM completed but output missing: $out_dir/${SRR}.genes.results"
		ls -la "$out_dir" 2>&1 | while IFS= read -r line; do _plog "ERROR" "  $line"; done
		return 1
	fi

	# Cleanup BAM files unless explicitly kept
	if [[ "${keep_bam_global:-n}" != "y" ]]; then
		rm -f "$out_dir/${SRR}.transcript.bam" "$out_dir/${SRR}.genome.bam" \
			  "$out_dir/${SRR}.transcript.sorted.bam" "$out_dir/${SRR}.transcript.sorted.bam.bai"
	fi

	_plog "INFO" "Completed successfully"
	return 0
}
export -f _rsem_parallel_worker

# Parallel sample quantification (requires USE_GNU_PARALLEL=TRUE)
_rsem_quantify_parallel() {
	local rsem_idx="$1"
	local quant_root="$2"
	local arr_name="${3:-}"

	# Expand array from indirect reference
	local valid_samples=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && valid_samples+=("$s")
		done
	fi

	local num_samples=${#valid_samples[@]}
	if [[ $num_samples -eq 0 ]]; then
		log_error "[RSEM QUANT] No valid samples to process!"
		return 1
	fi

	local parallel_jobs="${MAX_PARALLEL_SAMPLES:-2}"
	local threads_per_job="${THREADS_PER_RSEM_JOB:-$((THREADS / parallel_jobs))}"
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	log_step "[RSEM QUANT] Running PARALLEL quantification: $num_samples samples, $parallel_jobs concurrent jobs, $threads_per_job threads/job"
	log_info "[RSEM QUANT] Sample list: ${valid_samples[*]}"

	# Set up environment for parallel subshells using shared utility
	_prepare_parallel_env
	export rsem_idx quant_root threads_per_job OVERWRITE_MODE BOWTIE2_MODE RSEM_STRANDEDNESS RSEM_SEED

	printf "%s\n" "${valid_samples[@]}" | parallel \
		--env PATH \
		--env CONDA_PREFIX \
		--env CONDA_DEFAULT_ENV \
		--env CONDA_EXE \
		--env abs_trim_dir_root \
		--env abs_error_warn_file \
		--env keep_bam_global \
		--env rsem_idx \
		--env quant_root \
		--env threads_per_job \
		--env OVERWRITE_MODE \
		--env BOWTIE2_MODE \
		--env RSEM_STRANDEDNESS \
		--env RSEM_SEED \
		-j "$parallel_jobs" \
		--halt soon,fail=1 \
		--joblog "$quant_root/parallel_rsem.log" \
		--progress \
		_rsem_parallel_worker {}

	local parallel_exit=$?

	# Report results — count only files for the current batch (avoids inflation from stale runs)
	local successful=0
	for s in "${valid_samples[@]}"; do
		[[ -f "$quant_root/$s/${s}.genes.results" ]] && successful=$((successful + 1))
	done
	log_info "[RSEM QUANT] Parallel quantification complete: $successful/$num_samples samples succeeded"

	if [[ -f "$quant_root/parallel_rsem.log" ]]; then
		local failed
		failed=$(awk 'NR>1 && $7!=0' "$quant_root/parallel_rsem.log" | wc -l)
		[[ $failed -gt 0 ]] && log_warn "[RSEM QUANT] $failed sample(s) failed - check $quant_root/parallel_rsem.log"
	fi

	return $parallel_exit
}

# ==============================================================================
# MATRIX GENERATION
# ==============================================================================

_create_rsem_matrices() {
	local fasta="$1"
	local tag="$2"
	local quant_root="$3"
	local matrix_dir="$4"
	local arr_name="${5:-}"

	# Expand array from indirect reference
	local samples=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && samples+=("$s")
		done
	fi

	log_step "Generating gene and transcript matrices (RSEM)"

	# Check or create gene_trans_map
	local gene_trans_map="${fasta}.gene_trans_map"
	if [[ ! -s "$gene_trans_map" ]]; then
		log_info "[RSEM MATRIX] Creating gene-transcript mapping file..."
		if ! create_gene_trans_map "$fasta" "$gene_trans_map"; then
			log_warn "[RSEM MATRIX] create_gene_trans_map returned non-zero — checking output file"
		fi
		if [[ ! -s "$gene_trans_map" ]]; then
			log_warn "[RSEM MATRIX] gene_trans_map is empty or creation failed — skipping abundance_estimates_to_matrix.pl, using manual fallback"
			_create_manual_rsem_matrix "$quant_root" "$matrix_dir" samples[@]
			_prepare_rsem_deseq2_output "$tag" "$quant_root" "$matrix_dir" samples[@] || return 1
			return 0
		fi
	fi

	# Generate matrices
	if command -v abundance_estimates_to_matrix.pl >/dev/null 2>&1; then
		# Build explicit file list from samples array to avoid picking up stale results from other runs
		local rsem_result_files=()
		for s in "${samples[@]}"; do
			[[ -f "$quant_root/$s/${s}.genes.results" ]] && rsem_result_files+=("$quant_root/$s/${s}.genes.results")
		done
		if [[ ${#rsem_result_files[@]} -eq 0 ]]; then
			log_warn "No RSEM result files found — skipping abundance_estimates_to_matrix.pl, using manual fallback"
			_create_manual_rsem_matrix "$quant_root" "$matrix_dir" samples[@]
		else
			run_with_space_time_log abundance_estimates_to_matrix.pl \
				--est_method RSEM \
				--gene_trans_map "$gene_trans_map" \
				--out_prefix "$matrix_dir/genes" \
				--name_sample_by_basedir "${rsem_result_files[@]}" || {
				log_warn "abundance_estimates_to_matrix.pl failed. Creating manual count matrix..."
				_create_manual_rsem_matrix "$quant_root" "$matrix_dir" samples[@]
			}
		fi
	else
		log_warn "abundance_estimates_to_matrix.pl not found. Creating manual count matrix..."
		_create_manual_rsem_matrix "$quant_root" "$matrix_dir" samples[@]
	fi

	# Prepare DESeq2 outputs
	_prepare_rsem_deseq2_output "$tag" "$quant_root" "$matrix_dir" samples[@] || return 1
}


_create_manual_rsem_matrix() {
	local quant_root="$1"
	local matrix_dir="$2"
	local arr_name="${3:-}"

	# Expand array from indirect reference
	local srr_list=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && srr_list+=("$s")
		done
	fi

	local temp_gene_ids="$matrix_dir/temp_gene_ids.txt"

	local first_sample=""
	for SRR in "${srr_list[@]}"; do
		if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
			first_sample="$SRR"
			awk 'NR>1 {print $1}' "$quant_root/$SRR/${SRR}.genes.results" > "$temp_gene_ids"
			break
		fi
	done

	if [[ -z "$first_sample" ]]; then
		log_warn "[RSEM MATRIX] No samples have RSEM results — cannot create manual matrix"
		return 1
	fi

	local num_genes
	num_genes=$(wc -l < "$temp_gene_ids")

	for SRR in "${srr_list[@]}"; do
		if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
			# Preserve fractional expected counts — tximport requires non-integer values
			awk 'NR>1 {print $5}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_counts.tmp"
			awk 'NR>1 {print $6}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_tpm.tmp"
			awk 'NR>1 {print $7}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_fpkm.tmp"
		else
			log_warn "[RSEM MATRIX] Missing results for $SRR — filling with zeros in count matrix"
			yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_counts.tmp"
			yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_tpm.tmp"
			yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_fpkm.tmp"
		fi
	done

	# Build ordered file lists matching srr_list order (avoids glob sort mismatch)
	local count_files=() tpm_files=() fpkm_files=()
	for SRR in "${srr_list[@]}"; do
		count_files+=("$matrix_dir/${SRR}_counts.tmp")
		tpm_files+=("$matrix_dir/${SRR}_tpm.tmp")
		fpkm_files+=("$matrix_dir/${SRR}_fpkm.tmp")
	done

	# Create count matrix
	echo -n "gene_id" > "$matrix_dir/genes.counts.matrix"
	for SRR in "${srr_list[@]}"; do echo -ne "\t$SRR" >> "$matrix_dir/genes.counts.matrix"; done
	echo "" >> "$matrix_dir/genes.counts.matrix"
	paste "$temp_gene_ids" "${count_files[@]}" >> "$matrix_dir/genes.counts.matrix"

	# Create TPM matrix
	echo -n "gene_id" > "$matrix_dir/genes.TPM.not_cross_norm"
	for SRR in "${srr_list[@]}"; do echo -ne "\t$SRR" >> "$matrix_dir/genes.TPM.not_cross_norm"; done
	echo "" >> "$matrix_dir/genes.TPM.not_cross_norm"
	paste "$temp_gene_ids" "${tpm_files[@]}" >> "$matrix_dir/genes.TPM.not_cross_norm"

	# Create FPKM matrix
	echo -n "gene_id" > "$matrix_dir/genes.FPKM.not_cross_norm"
	for SRR in "${srr_list[@]}"; do echo -ne "\t$SRR" >> "$matrix_dir/genes.FPKM.not_cross_norm"; done
	echo "" >> "$matrix_dir/genes.FPKM.not_cross_norm"
	paste "$temp_gene_ids" "${fpkm_files[@]}" >> "$matrix_dir/genes.FPKM.not_cross_norm"

	rm -f "$temp_gene_ids" "$matrix_dir"/*_counts.tmp "$matrix_dir"/*_tpm.tmp "$matrix_dir"/*_fpkm.tmp
}

_prepare_rsem_deseq2_output() {
	local tag="$1"
	local quant_root="$2"
	local matrix_dir="$3"
	local arr_name="${4:-}"

	# Expand array from indirect reference
	local srr_list=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && srr_list+=("$s")
		done
	fi

	log_step "Preparing DESeq2-compatible count matrix for RSEM pipeline"
	log_info "[NOTE] RSEM reports expected counts; for DESeq2, prefer tximport (script will be generated)."

	local deseq2_dir="$matrix_dir/deseq2_input"
	local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
	local sample_metadata="$deseq2_dir/sample_metadata.csv"
	mkdir -p "$deseq2_dir"

	# Verify quantifications
	local quant_count=0
	for SRR in "${srr_list[@]}"; do
		[[ -f "$quant_root/$SRR/${SRR}.genes.results" ]] && quant_count=$((quant_count + 1))
	done

	[[ $quant_count -lt 2 ]] && { log_error "Insufficient RSEM quantifications (found: $quant_count, need: ≥2)"; return 1; }
	log_info "[RSEM] Found $quant_count samples with successful quantifications"

	# Convert to CSV (always regenerate if source matrix is newer)
	if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
		if [[ ! -f "$gene_count_matrix" || "$matrix_dir/genes.counts.matrix" -nt "$gene_count_matrix" ]]; then
			log_info "[RSEM MATRIX] Converting count matrix to CSV format..."
			sed 's/\t/,/g' "$matrix_dir/genes.counts.matrix" | sed '1s/gene_id/Gene_ID/' > "$gene_count_matrix"
		fi
	fi

	# Create sample metadata (regenerate if missing or samples changed)
	create_sample_metadata "$sample_metadata" "${srr_list[@]}"

	# Generate tximport script (regenerate if missing)
	local tximport_script="$deseq2_dir/run_tximport_rsem.R"
	if [[ ! -f "$tximport_script" || "${OVERWRITE_MODE:-skip}" == "overwrite" ]]; then
		generate_tximport_script "rsem" "$quant_root" "$tximport_script" "$sample_metadata"
	fi

	# Create TPM and FPKM matrices (always regenerate if source is newer)
	if [[ -f "$matrix_dir/genes.TPM.not_cross_norm" ]]; then
		local tpm_matrix="$deseq2_dir/gene_tpm_matrix.csv"
		if [[ ! -f "$tpm_matrix" || "$matrix_dir/genes.TPM.not_cross_norm" -nt "$tpm_matrix" ]]; then
			sed 's/\t/,/g' "$matrix_dir/genes.TPM.not_cross_norm" | sed '1s/gene_id/Gene_ID/' > "$tpm_matrix"
		fi
	fi

	if [[ -f "$matrix_dir/genes.FPKM.not_cross_norm" ]]; then
		local fpkm_matrix="$deseq2_dir/gene_fpkm_matrix.csv"
		if [[ ! -f "$fpkm_matrix" || "$matrix_dir/genes.FPKM.not_cross_norm" -nt "$fpkm_matrix" ]]; then
			sed 's/\t/,/g' "$matrix_dir/genes.FPKM.not_cross_norm" | sed '1s/gene_id/Gene_ID/' > "$fpkm_matrix"
		fi
	fi

	# Create summary
	_create_rsem_summary "$tag" "$quant_root" "$deseq2_dir" srr_list[@]

	# Validate
	[[ -f "$gene_count_matrix" ]] && validate_count_matrix "$gene_count_matrix" "gene" 2

	log_info "DESeq2 input files:"
	log_info "  - Gene count matrix: $gene_count_matrix"
	log_info "  - Sample metadata: $sample_metadata"
}

_create_rsem_summary() {
	local tag="$1"
	local quant_root="$2"
	local deseq2_dir="$3"
	local arr_name="${4:-}"

	# Expand array from indirect reference
	local srr_list=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && srr_list+=("$s")
		done
	fi

	local summary_file="$deseq2_dir/rsem_summary.txt"

	local LOW_ALIGN_THRESHOLD=50  # Flag samples below this alignment rate (%)
	local outlier_samples=()

	{
		echo "==================================================================="
		echo "RSEM Quantification Summary for $tag"
		echo "==================================================================="
		echo "Date: $(date)"
		echo "Samples processed: ${#srr_list[@]}"
		echo "Method: RSEM with Bowtie2 alignment (strandedness: ${RSEM_STRANDEDNESS:-none}, sensitivity: ${BOWTIE2_MODE:-sensitive})"
		echo ""
		echo "Per-sample statistics:"
		echo "-------------------------------------------------------------------"

		for SRR in "${srr_list[@]}"; do
			if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
				local total=$(awk 'NR>1' "$quant_root/$SRR/${SRR}.genes.results" | wc -l)
				local expressed=$(awk 'NR>1 && $5>0' "$quant_root/$SRR/${SRR}.genes.results" | wc -l)
				local counts=$(awk 'NR>1 {sum+=$5} END {print int(sum)}' "$quant_root/$SRR/${SRR}.genes.results")

				# Extract Bowtie2 alignment rate from RSEM log
				local align_rate="N/A"
				local rsem_log="$quant_root/$SRR/${SRR}.rsem.log"
				if [[ -f "$rsem_log" ]]; then
					# Bowtie2 reports "XX.XX% overall alignment rate" in its summary
					align_rate=$(grep -oP '[0-9]+\.[0-9]+(?=% overall alignment rate)' "$rsem_log" | tail -1)
					if [[ -z "$align_rate" ]]; then
						align_rate="N/A"
					elif (( $(echo "$align_rate < $LOW_ALIGN_THRESHOLD" | bc -l) )); then
						outlier_samples+=("$SRR ($align_rate%)")
					fi
				fi

				echo "$SRR: $expressed/$total expressed genes, $counts expected counts, alignment rate: ${align_rate}%"
			fi
		done

		if [[ ${#outlier_samples[@]} -gt 0 ]]; then
			echo ""
			echo "WARNING: Low alignment rate samples (<${LOW_ALIGN_THRESHOLD}%):"
			for s in "${outlier_samples[@]}"; do
				echo "  - $s"
			done
			echo "These samples may have contamination, wrong reference, or quality issues."
		fi
	} > "$summary_file"

	# Also log outlier warnings to the main pipeline log
	if [[ ${#outlier_samples[@]} -gt 0 ]]; then
		log_warn "[RSEM SUMMARY] ${#outlier_samples[@]} sample(s) with low Bowtie2 alignment rate (<${LOW_ALIGN_THRESHOLD}%):"
		for s in "${outlier_samples[@]}"; do
			log_warn "  - $s"
		done
	fi

	log_info "[RSEM SUMMARY] Summary saved to: $summary_file"
}

# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

# Check if GNU Parallel should be used for sample processing
_rsem_should_use_parallel() {
	[[ "${USE_GNU_PARALLEL:-FALSE}" != "TRUE" ]] && return 1
	command -v parallel >/dev/null 2>&1 || return 1
	return 0
}

# Process a single sample (shared by sequential and parallel modes)
_rsem_process_single_sample() {
	local SRR="$1"
	local rsem_idx="$2"
	local quant_root="$3"
	local threads_to_use="$4"

	local out_dir="$quant_root/$SRR"
	mkdir -p "$out_dir"

	# Skip if already processed
	if [[ -f "$out_dir/${SRR}.genes.results" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[RSEM QUANT] RSEM results for $SRR already exist. Skipping."
		return 0
	fi

	# Find trimmed reads
	find_trimmed_fastq "$SRR"
	if [[ -z "$trimmed1" ]]; then
		log_warn "Missing trimmed reads for $SRR in $TRIM_DIR_ROOT/$SRR"
		ls -la "$TRIM_DIR_ROOT/$SRR" 2>&1 | while IFS= read -r line; do log_warn "  $line"; done
		return 1
	fi

	local rsem_log="$out_dir/${SRR}.rsem.log"
	local rsem_exit_code
	if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
		log_step "Running Bowtie2 + RSEM (paired-end) for $SRR (threads: $threads_to_use, strandedness: ${RSEM_STRANDEDNESS:-none})"
		rsem-calculate-expression \
			--paired-end \
			--bowtie2 \
			--bowtie2-sensitivity-level "${BOWTIE2_MODE:-sensitive}" \
			--strandedness "${RSEM_STRANDEDNESS:-none}" \
			--seed "$RSEM_SEED" \
			--num-threads "$threads_to_use" \
			"$trimmed1" "$trimmed2" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$rsem_log"
		rsem_exit_code=${PIPESTATUS[0]}
	else
		log_step "Running Bowtie2 + RSEM (single-end) for $SRR (threads: $threads_to_use, strandedness: ${RSEM_STRANDEDNESS:-none})"
		rsem-calculate-expression \
			--bowtie2 \
			--bowtie2-sensitivity-level "${BOWTIE2_MODE:-sensitive}" \
			--strandedness "${RSEM_STRANDEDNESS:-none}" \
			--seed "$RSEM_SEED" \
			--num-threads "$threads_to_use" \
			"$trimmed1" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$rsem_log"
		rsem_exit_code=${PIPESTATUS[0]}
	fi

	if [[ $rsem_exit_code -ne 0 ]]; then
		log_error "[RSEM QUANT] RSEM failed for $SRR (exit code: $rsem_exit_code)"
		log_error "[RSEM QUANT] Check log: $rsem_log"
		[[ -f "$rsem_log" ]] && tail -20 "$rsem_log" | while IFS= read -r line; do log_error "  $line"; done
		return $rsem_exit_code
	fi

	# Verify output exists even when RSEM exits 0 (edge case: disk full, interrupted write)
	if [[ ! -f "$out_dir/${SRR}.genes.results" ]]; then
		log_error "[RSEM QUANT] RSEM completed but output missing: $out_dir/${SRR}.genes.results"
		ls -la "$out_dir" 2>&1 | while IFS= read -r line; do log_error "  $line"; done
		return 1
	fi

	log_file_size "$out_dir/${SRR}.genes.results" "RSEM gene results - $SRR"

	# Cleanup BAM files
	if [[ "${keep_bam_global:-n}" != "y" ]]; then
		log_info "[CLEANUP] Deleting RSEM BAM files for $SRR to save disk space"
		rm -f "$out_dir/${SRR}.transcript.bam" "$out_dir/${SRR}.genome.bam" \
			  "$out_dir/${SRR}.transcript.sorted.bam" "$out_dir/${SRR}.transcript.sorted.bam.bai"
	fi

	return 0
}
