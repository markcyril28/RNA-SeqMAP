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
# RSEM CONFIGURATION - IMPORTANT PARAMETERS AT TOP
# ==============================================================================

# Bowtie2 alignment mode for RSEM (options: very_sensitive, sensitive, fast, very_fast)
# Note: RSEM uses end-to-end mode by default, not local mode
BOWTIE2_MODE="${BOWTIE2_MODE:-sensitive}"

# Parallel sample processing configuration
# MAX_PARALLEL_SAMPLES: Number of samples to process in parallel (default: 2)
# THREADS_PER_RSEM_JOB: Threads allocated per RSEM job (auto-calculated if not set)
MAX_PARALLEL_SAMPLES="${MAX_PARALLEL_SAMPLES:-2}"
THREADS_PER_RSEM_JOB="${THREADS_PER_RSEM_JOB:-$((THREADS / MAX_PARALLEL_SAMPLES))}"
[[ $THREADS_PER_RSEM_JOB -lt 1 ]] && THREADS_PER_RSEM_JOB=1

# ==============================================================================
# HELPER: CHECK IF GNU PARALLEL SHOULD BE USED
# ==============================================================================

_rsem_should_use_parallel() {
	# Check if parallel processing is enabled and available
	[[ "${USE_GNU_PARALLEL:-FALSE}" != "TRUE" ]] && return 1
	command -v parallel >/dev/null 2>&1 || return 1
	return 0
}

# ==============================================================================
# HELPER: PROCESS SINGLE SAMPLE (used by both sequential and parallel modes)
# ==============================================================================

_rsem_process_single_sample() {
	local SRR="$1"
	local rsem_idx="$2"
	local quant_root="$3"
	local bowtie2_mode="$4"
	local threads_to_use="$5"
	
	local out_dir="$quant_root/$SRR"
	mkdir -p "$out_dir"
	
	# Skip if already processed
	if [[ -f "$out_dir/${SRR}.genes.results" ]]; then
		log_info "[RSEM QUANT] RSEM results for $SRR already exist. Skipping."
		return 0
	fi
	
	# Find trimmed reads
	log_info "[RSEM DEBUG] Looking for trimmed reads in: $TRIM_DIR_ROOT/$SRR"
	find_trimmed_fastq "$SRR"
	if [[ -z "$trimmed1" ]]; then
		log_warn "Missing trimmed reads for $SRR in $TRIM_DIR_ROOT/$SRR"
		log_warn "[DEBUG] Directory contents:"
		ls -la "$TRIM_DIR_ROOT/$SRR" 2>&1 | while read line; do log_warn "  $line"; done
		return 1
	fi
	
	log_step "Running Bowtie2 ($bowtie2_mode) + RSEM for $SRR (threads: $threads_to_use)"
	log_info "[RSEM DEBUG] trimmed1=$trimmed1"
	[[ -n "$trimmed2" ]] && log_info "[RSEM DEBUG] trimmed2=$trimmed2"
	log_info "[RSEM DEBUG] rsem_idx=$rsem_idx"
	
	local rsem_exit_code=0
	if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
		# Paired-end reads
		log_info "[RSEM QUANT] Using paired-end reads for $SRR"
		rsem-calculate-expression \
			--paired-end \
			--bowtie2 \
			--bowtie2-sensitivity-level "$bowtie2_mode" \
			--num-threads "$threads_to_use" \
			"$trimmed1" "$trimmed2" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$out_dir/${SRR}.rsem.log" || rsem_exit_code=$?
	else
		# Single-end reads
		log_info "[RSEM QUANT] Using single-end reads for $SRR"
		rsem-calculate-expression \
			--bowtie2 \
			--bowtie2-sensitivity-level "$bowtie2_mode" \
			--num-threads "$threads_to_use" \
			"$trimmed1" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$out_dir/${SRR}.rsem.log" || rsem_exit_code=$?
	fi
	
	if [[ $rsem_exit_code -ne 0 ]]; then
		log_error "[RSEM QUANT] RSEM failed for $SRR (exit code: $rsem_exit_code)"
		log_error "[RSEM QUANT] Check log: $out_dir/${SRR}.rsem.log"
		# Show last 20 lines of RSEM log for quick diagnosis
		if [[ -f "$out_dir/${SRR}.rsem.log" ]]; then
			log_error "[RSEM QUANT] Last 20 lines of RSEM output:"
			tail -20 "$out_dir/${SRR}.rsem.log" | while read line; do log_error "  $line"; done
		fi
		return $rsem_exit_code
	fi
	
	log_file_size "$out_dir/${SRR}.genes.results" "RSEM gene results - $SRR"
	
	# Cleanup BAM files
	if [[ "$keep_bam_global" != "y" ]]; then
		log_info "[CLEANUP] Deleting RSEM BAM files for $SRR to save disk space"
		rm -f "$out_dir/${SRR}.transcript.bam" "$out_dir/${SRR}.genome.bam" \
			  "$out_dir/${SRR}.transcript.sorted.bam" "$out_dir/${SRR}.transcript.sorted.bam.bai"
	fi
	
	return 0
}

# ==============================================================================
# BOWTIE2 + RSEM PIPELINE
# ==============================================================================

bowtie2_rsem_pipeline() {
	local fasta="" rnaseq_list=() bowtie2_mode="${BOWTIE2_MODE}"
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--BOWTIE2_MODE) bowtie2_mode="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do rnaseq_list+=("$1"); shift; done;;
			*) log_error "Unknown arg: $1"; return 1;;
		esac
	done

	[[ -z "$fasta" ]] && { log_error "Usage: --FASTA genes.fa [--BOWTIE2_MODE sensitive]"; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	
	# Validate bowtie2 mode - RSEM accepts: very_sensitive, sensitive, fast, very_fast
	# Convert hyphenated/local modes to underscore format for RSEM compatibility
	case "$bowtie2_mode" in
		very-sensitive-local|very-sensitive|very_sensitive)
			bowtie2_mode="very_sensitive"
			;;
		sensitive-local|sensitive)
			bowtie2_mode="sensitive"
			;;
		fast-local|fast)
			bowtie2_mode="fast"
			;;
		very-fast-local|very-fast|very_fast)
			bowtie2_mode="very_fast"
			;;
		*)
			log_warn "[BOWTIE2] Unknown mode '$bowtie2_mode', defaulting to 'sensitive'"
			bowtie2_mode="sensitive"
			;;
	esac
	log_info "[BOWTIE2] Using RSEM alignment mode: $bowtie2_mode"
	
	# Convert line endings if dos2unix is available
	command -v dos2unix >/dev/null 2>&1 && dos2unix "$fasta" 2>/dev/null || true

	local tag="$(basename "${fasta%.*}")"
	local rsem_idx="$RSEM_INDEX_ROOT/${tag}_rsem"
	local quant_root="$RSEM_QUANT_ROOT/$tag"
	local matrix_dir="$RSEM_MATRIX_ROOT/$tag"

	mkdir -p "$RSEM_INDEX_ROOT" "$quant_root" "$matrix_dir"

	# BUILD RSEM REFERENCE
	if [[ -f "${rsem_idx}.grp" ]]; then
		log_info "[RSEM INDEX] RSEM reference already exists. Skipping."
	else
		log_step "Building RSEM reference for $tag"
		log_file_size "$fasta" "Input FASTA for RSEM index - $tag"
		run_with_space_time_log --input "$fasta" --output "$RSEM_INDEX_ROOT" \
			rsem-prepare-reference --bowtie2 "$fasta" "$rsem_idx"
		log_file_size "$RSEM_INDEX_ROOT" "RSEM index output - $tag"
	fi

	# QUANTIFY SAMPLES - PARALLEL OR SEQUENTIAL
	if _rsem_should_use_parallel && [[ ${#rnaseq_list[@]} -gt 1 ]]; then
		_rsem_quantify_parallel "$rsem_idx" "$quant_root" "$bowtie2_mode" rnaseq_list[@]
	else
		_rsem_quantify_sequential "$rsem_idx" "$quant_root" "$bowtie2_mode" rnaseq_list[@]
	fi

	# GENERATE MATRICES
	_create_rsem_matrices "$fasta" "$tag" "$quant_root" "$matrix_dir" "$bowtie2_mode" rnaseq_list[@]
	
	log_step "COMPLETED: Bowtie2-RSEM pipeline for $tag (mode: $bowtie2_mode)"
}

# ==============================================================================
# SEQUENTIAL SAMPLE QUANTIFICATION
# ==============================================================================

_rsem_quantify_sequential() {
	local rsem_idx="$1"
	local quant_root="$2"
	local bowtie2_mode="$3"
	local arr_name="${4:-}"
	
	# Expand array from indirect reference
	local samples=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && samples+=("$s")
		done
	fi
	
	log_info "[RSEM QUANT] Running SEQUENTIAL quantification for ${#samples[@]} samples (threads per job: $THREADS)"
	
	for SRR in "${samples[@]}"; do
		_rsem_process_single_sample "$SRR" "$rsem_idx" "$quant_root" "$bowtie2_mode" "$THREADS"
	done
}

# ==============================================================================
# PARALLEL SAMPLE QUANTIFICATION (GNU Parallel)
# ==============================================================================

_rsem_quantify_parallel() {
	local rsem_idx="$1"
	local quant_root="$2"
	local bowtie2_mode="$3"
	local arr_name="${4:-}"
	
	# Debug: Log received parameters
	log_info "[RSEM DEBUG] rsem_idx=$rsem_idx"
	log_info "[RSEM DEBUG] quant_root=$quant_root"
	log_info "[RSEM DEBUG] bowtie2_mode=$bowtie2_mode"
	log_info "[RSEM DEBUG] arr_name='$arr_name'"
	
	# Expand array from indirect reference
	local valid_samples=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && valid_samples+=("$s")
		done
	fi
	
	# Debug: Log expanded samples
	log_info "[RSEM DEBUG] valid_samples count: ${#valid_samples[@]}"
	log_info "[RSEM DEBUG] valid_samples content: ${valid_samples[*]}"
	
	local num_samples=${#valid_samples[@]}
	if [[ $num_samples -eq 0 ]]; then
		log_error "[RSEM QUANT] No valid samples to process!"
		log_error "[RSEM QUANT] arr_name='$arr_name'"
		log_error "[RSEM QUANT] Check that --RNASEQ_LIST samples are being passed correctly"
		return 1
	fi
	
	local parallel_jobs="${MAX_PARALLEL_SAMPLES:-2}"
	local threads_per_job="${THREADS_PER_RSEM_JOB:-$((THREADS / parallel_jobs))}"
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1
	
	log_step "[RSEM QUANT] Running PARALLEL quantification: $num_samples samples, $parallel_jobs concurrent jobs, $threads_per_job threads/job"
	log_info "[RSEM QUANT] Sample list: ${valid_samples[*]}"
	
	# Export required variables and functions for parallel subshells
	export PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE
	export TRIM_DIR_ROOT LOG_FILE TIME_DIR TIME_FILE ERROR_WARN_FILE
	export keep_bam_global
	
	# Convert TRIM_DIR_ROOT to absolute path for parallel subshells
	local abs_trim_dir_root="$TRIM_DIR_ROOT"
	[[ "$abs_trim_dir_root" != /* ]] && abs_trim_dir_root="$(pwd)/$abs_trim_dir_root"
	
	# Convert ERROR_WARN_FILE to absolute path for parallel subshells
	local abs_error_warn_file="${ERROR_WARN_FILE:-}"
	[[ -n "$abs_error_warn_file" && "$abs_error_warn_file" != /* ]] && abs_error_warn_file="$(pwd)/$abs_error_warn_file"
	
	export rsem_idx quant_root bowtie2_mode threads_per_job abs_trim_dir_root abs_error_warn_file
	
	# Note: We inline find_trimmed_fastq logic in the worker to avoid function export issues
	
	# Define parallel worker function
	_rsem_parallel_worker() {
		local SRR="$1"
		local ts=$(date '+%Y-%m-%d %H:%M:%S')
		
		# Helper to log to error file
		_log_err() {
			local level="$1"; shift
			echo "[$ts] [$level] [RSEM-$SRR] $*"
			[[ -n "$abs_error_warn_file" ]] && echo "[$ts] [$level] [RSEM-$SRR] $*" >> "$abs_error_warn_file"
		}
		
		# Validate SRR ID
		if [[ -z "$SRR" ]]; then
			_log_err "ERROR" "Empty SRR ID received - skipping"
			return 1
		fi
		
		# Reactivate conda environment in subshell if needed
		if [[ -n "$CONDA_PREFIX" && -n "$CONDA_EXE" ]]; then
			source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local out_dir="$quant_root/$SRR"
		mkdir -p "$out_dir"
		
		# Skip if already processed
		if [[ -f "$out_dir/${SRR}.genes.results" ]]; then
			echo "[RSEM QUANT] Results for $SRR already exist. Skipping."
			return 0
		fi
		
		# Inline find_trimmed_fastq logic (avoids function export issues)
		local trim_dir="$abs_trim_dir_root/$SRR"
		local trimmed1="" trimmed2=""
		
		_log_err "DEBUG" "Looking for trimmed reads in: $trim_dir"
		
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
			# Fallback: try to find any val_1/val_2 files
			for f in "$trim_dir"/${SRR}*val_1*.fq* "$trim_dir"/${SRR}*val_1*.gz; do
				[[ -f "$f" ]] && { trimmed1="$f"; break; }
			done
			for f in "$trim_dir"/${SRR}*val_2*.fq* "$trim_dir"/${SRR}*val_2*.gz; do
				[[ -f "$f" ]] && { trimmed2="$f"; break; }
			done
		fi
		
		if [[ -z "$trimmed1" ]]; then
			_log_err "ERROR" "Missing trimmed reads for $SRR in $trim_dir"
			_log_err "ERROR" "Directory contents: $(ls -la "$trim_dir" 2>&1 || echo 'Directory does not exist')"
			return 1
		fi
		
		_log_err "INFO" "Processing $SRR with $threads_per_job threads"
		_log_err "DEBUG" "trimmed1=$trimmed1"
		[[ -n "$trimmed2" ]] && _log_err "DEBUG" "trimmed2=$trimmed2"
		_log_err "DEBUG" "rsem_idx=$rsem_idx"
		
		# Verify RSEM index exists
		if [[ ! -f "${rsem_idx}.grp" ]]; then
			_log_err "ERROR" "RSEM index not found: ${rsem_idx}.grp"
			_log_err "ERROR" "Index directory contents: $(ls -la "$(dirname "$rsem_idx")" 2>&1 || echo 'Directory does not exist')"
			return 1
		fi
		
		local rsem_log="$out_dir/${SRR}.rsem.log"
		local rsem_exit_code=0
		if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
			rsem-calculate-expression \
				--paired-end \
				--bowtie2 \
				--bowtie2-sensitivity-level "$bowtie2_mode" \
				--num-threads "$threads_per_job" \
				"$trimmed1" "$trimmed2" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$rsem_log" || rsem_exit_code=$?
		else
			rsem-calculate-expression \
				--bowtie2 \
				--bowtie2-sensitivity-level "$bowtie2_mode" \
				--num-threads "$threads_per_job" \
				"$trimmed1" "$rsem_idx" "$out_dir/$SRR" 2>&1 | tee "$rsem_log" || rsem_exit_code=$?
		fi
		
		if [[ $rsem_exit_code -ne 0 ]]; then
			_log_err "ERROR" "RSEM failed (exit code: $rsem_exit_code)"
			_log_err "ERROR" "RSEM log: $rsem_log"
			# Log last 20 lines of RSEM output to error log
			if [[ -f "$rsem_log" ]]; then
				_log_err "ERROR" "Last 20 lines of RSEM output:"
				tail -20 "$rsem_log" 2>/dev/null | while IFS= read -r line; do
					_log_err "ERROR" "  $line"
				done
			fi
			return $rsem_exit_code
		fi
		
		# Verify output was created
		if [[ ! -f "$out_dir/${SRR}.genes.results" ]]; then
			_log_err "ERROR" "RSEM completed but output file not found: $out_dir/${SRR}.genes.results"
			_log_err "ERROR" "Output directory contents: $(ls -la "$out_dir" 2>&1)"
			if [[ -f "$rsem_log" ]]; then
				_log_err "ERROR" "Last 20 lines of RSEM output:"
				tail -20 "$rsem_log" 2>/dev/null | while IFS= read -r line; do
					_log_err "ERROR" "  $line"
				done
			fi
			return 1
		fi
		
		# Cleanup BAM files
		if [[ "$keep_bam_global" != "y" ]]; then
			rm -f "$out_dir/${SRR}.transcript.bam" "$out_dir/${SRR}.genome.bam" \
				  "$out_dir/${SRR}.transcript.sorted.bam" "$out_dir/${SRR}.transcript.sorted.bam.bai"
		fi
		
		_log_err "INFO" "Completed successfully"
		return 0
	}
	export -f _rsem_parallel_worker
	
	# Run parallel quantification with validated sample list
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
		--env bowtie2_mode \
		--env threads_per_job \
		-j "$parallel_jobs" \
		--joblog "$quant_root/parallel_rsem.log" \
		--progress \
		_rsem_parallel_worker {}
	
	local parallel_exit=$?
	
	# Report results
	local successful=$(find "$quant_root" -name "*.genes.results" 2>/dev/null | wc -l)
	log_info "[RSEM QUANT] Parallel quantification complete: $successful/$num_samples samples succeeded"
	
	if [[ -f "$quant_root/parallel_rsem.log" ]]; then
		local failed=$(awk 'NR>1 && $7!=0' "$quant_root/parallel_rsem.log" | wc -l)
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
	local bowtie2_mode="$5"
	local arr_name="${6:-}"

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
	if [[ ! -f "$gene_trans_map" ]]; then
		log_info "[RSEM MATRIX] Creating gene-transcript mapping file..."
		_create_gene_trans_map_rsem "$fasta" "$gene_trans_map"
	fi
	
	# Generate matrices
	if command -v abundance_estimates_to_matrix.pl >/dev/null 2>&1; then
		run_with_space_time_log abundance_estimates_to_matrix.pl \
			--est_method RSEM \
			--gene_trans_map "$gene_trans_map" \
			--out_prefix "$matrix_dir/genes" \
			--name_sample_by_basedir "$quant_root"/*/*.genes.results || {
			log_warn "abundance_estimates_to_matrix.pl failed. Creating manual count matrix..."
			_create_manual_rsem_matrix "$quant_root" "$matrix_dir" samples[@]
		}
	else
		log_warn "abundance_estimates_to_matrix.pl not found. Creating manual count matrix..."
		_create_manual_rsem_matrix "$quant_root" "$matrix_dir" samples[@]
	fi
	
	# Prepare DESeq2 outputs
	_prepare_rsem_deseq2_output "$tag" "$quant_root" "$matrix_dir" "$bowtie2_mode" samples[@]
}

_create_gene_trans_map_rsem() {
	local fasta="$1"
	local output="$2"
	
	local is_trinity=false
	grep -q "^>TRINITY_" "$fasta" 2>/dev/null && is_trinity=true
	
	if [[ "$is_trinity" == "true" ]]; then
		log_info "[RSEM MATRIX] Detected Trinity assembly format"
		awk '/^>/ {
			trans = $1; gsub(/^>/, "", trans)
			gene = trans
			if (match(gene, /^(.+)_i[0-9]+$/, arr)) { gene = arr[1] }
			print gene "\t" trans
		}' "$fasta" > "$output"
	else
		grep "^>" "$fasta" | sed 's/^>//' | awk '{
			trans=$1
			if (match($0, /gene=([^ ]+)/, arr)) { gene=arr[1] }
			else if (match($0, /gene_id[=:]([^ ]+)/, arr)) { gene=arr[1] }
			else if (match(trans, /^([^|]+)\|/, arr)) { gene=arr[1] }
			else { gene=trans }
			print gene "\t" trans
		}' > "$output"
	fi
	
	local unique_genes=$(cut -f1 "$output" | sort -u | wc -l)
	local total_transcripts=$(wc -l < "$output")
	log_info "[RSEM MATRIX] Created gene-transcript map: $unique_genes genes, $total_transcripts transcripts"
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
	
	if [[ -n "$first_sample" ]]; then
		for SRR in "${srr_list[@]}"; do
			if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
				awk 'NR>1 {print int($5 + 0.5)}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_counts.tmp"
				awk 'NR>1 {print $6}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_tpm.tmp"
				awk 'NR>1 {print $7}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_fpkm.tmp"
			else
				local num_genes=$(wc -l < "$temp_gene_ids")
				yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_counts.tmp"
				yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_tpm.tmp"
				yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_fpkm.tmp"
			fi
		done
		
		# Create count matrix
		echo -n "gene_id" > "$matrix_dir/genes.counts.matrix"
		for SRR in "${srr_list[@]}"; do echo -ne "\t$SRR" >> "$matrix_dir/genes.counts.matrix"; done
		echo "" >> "$matrix_dir/genes.counts.matrix"
		paste "$temp_gene_ids" "$matrix_dir"/*_counts.tmp >> "$matrix_dir/genes.counts.matrix"
		
		# Create TPM matrix
		echo -n "gene_id" > "$matrix_dir/genes.TPM.not_cross_norm"
		for SRR in "${srr_list[@]}"; do echo -ne "\t$SRR" >> "$matrix_dir/genes.TPM.not_cross_norm"; done
		echo "" >> "$matrix_dir/genes.TPM.not_cross_norm"
		paste "$temp_gene_ids" "$matrix_dir"/*_tpm.tmp >> "$matrix_dir/genes.TPM.not_cross_norm"
		
		# Create FPKM matrix
		echo -n "gene_id" > "$matrix_dir/genes.FPKM.not_cross_norm"
		for SRR in "${srr_list[@]}"; do echo -ne "\t$SRR" >> "$matrix_dir/genes.FPKM.not_cross_norm"; done
		echo "" >> "$matrix_dir/genes.FPKM.not_cross_norm"
		paste "$temp_gene_ids" "$matrix_dir"/*_fpkm.tmp >> "$matrix_dir/genes.FPKM.not_cross_norm"
		
		rm -f "$temp_gene_ids" "$matrix_dir"/*_counts.tmp "$matrix_dir"/*_tpm.tmp "$matrix_dir"/*_fpkm.tmp
	fi
}

_prepare_rsem_deseq2_output() {
	local tag="$1"
	local quant_root="$2"
	local matrix_dir="$3"
	local bowtie2_mode="$4"
	local arr_name="${5:-}"

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
		[[ -f "$quant_root/$SRR/${SRR}.genes.results" ]] && ((quant_count++))
	done
	
	[[ $quant_count -lt 2 ]] && { log_error "Insufficient RSEM quantifications (found: $quant_count, need: ≥2)"; return 1; }
	log_info "[RSEM] Found $quant_count samples with successful quantifications"
	
	# Convert to CSV
	if [[ -f "$matrix_dir/genes.counts.matrix" && ! -f "$gene_count_matrix" ]]; then
		log_info "[RSEM MATRIX] Converting count matrix to CSV format..."
		sed 's/\t/,/g' "$matrix_dir/genes.counts.matrix" | sed '1s/gene_id/Gene_ID/' > "$gene_count_matrix"
	fi
	
	# Create sample metadata
	[[ ! -f "$sample_metadata" ]] && create_sample_metadata "$sample_metadata" srr_list[@]
	
	# Generate tximport script
	local tximport_script="$deseq2_dir/run_tximport_rsem.R"
	[[ ! -f "$tximport_script" ]] && generate_tximport_script "rsem" "$quant_root" "$tximport_script" "$sample_metadata"
	
	# Create TPM and FPKM matrices
	if [[ -f "$matrix_dir/genes.TPM.not_cross_norm" ]]; then
		local tpm_matrix="$deseq2_dir/gene_tpm_matrix.csv"
		[[ ! -f "$tpm_matrix" ]] && sed 's/\t/,/g' "$matrix_dir/genes.TPM.not_cross_norm" | sed '1s/gene_id/Gene_ID/' > "$tpm_matrix"
	fi
	
	if [[ -f "$matrix_dir/genes.FPKM.not_cross_norm" ]]; then
		local fpkm_matrix="$deseq2_dir/gene_fpkm_matrix.csv"
		[[ ! -f "$fpkm_matrix" ]] && sed 's/\t/,/g' "$matrix_dir/genes.FPKM.not_cross_norm" | sed '1s/gene_id/Gene_ID/' > "$fpkm_matrix"
	fi
	
	# Create summary
	_create_rsem_summary "$tag" "$quant_root" "$deseq2_dir" "$bowtie2_mode" srr_list[@]
	
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
	local bowtie2_mode="$4"
	local arr_name="${5:-}"
	
	# Expand array from indirect reference
	local srr_list=()
	if [[ -n "$arr_name" ]]; then
		local tmp_arr=("${!arr_name}")
		for s in "${tmp_arr[@]}"; do
			[[ -n "$s" ]] && srr_list+=("$s")
		done
	fi
	
	local summary_file="$deseq2_dir/rsem_summary.txt"
	[[ -f "$summary_file" ]] && return 0
	
	{
		echo "==================================================================="
		echo "RSEM Quantification Summary for $tag"
		echo "==================================================================="
		echo "Date: $(date)"
		echo "Samples processed: ${#srr_list[@]}"
		echo "Method: RSEM with Bowtie2 alignment ($bowtie2_mode)"
		echo ""
		echo "Bowtie2 mode options (RSEM uses end-to-end mode):"
		echo "  very_sensitive: Most thorough (slower, best accuracy)"
		echo "  sensitive: Balanced (default)"
		echo "  fast: Faster"
		echo "  very_fast: Fastest (less accurate)"
		echo ""
		echo "Per-sample statistics:"
		echo "-------------------------------------------------------------------"
		
		for SRR in "${srr_list[@]}"; do
			if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
				local total=$(awk 'NR>1' "$quant_root/$SRR/${SRR}.genes.results" | wc -l)
				local expressed=$(awk 'NR>1 && $5>0' "$quant_root/$SRR/${SRR}.genes.results" | wc -l)
				local counts=$(awk 'NR>1 {sum+=$5} END {print int(sum)}' "$quant_root/$SRR/${SRR}.genes.results")
				echo "$SRR: $expressed/$total expressed genes, $counts expected counts"
			fi
		done
	} > "$summary_file"
	
	log_info "[RSEM SUMMARY] Summary saved to: $summary_file"
}
