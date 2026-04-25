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
M5_RSEM_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/b_main_methods}"
if [[ -z "$SCRIPT_DIR" ]]; then
	SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
	SCRIPT_DIR="$(cd "$SCRIPT_DIR" 2>/dev/null && pwd)"
fi
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# RSEM CONFIGURATION - IMPORTANT PARAMETERS
# ==============================================================================

# Number of samples to process in parallel (requires USE_GNU_PARALLEL=TRUE)
# RSEM/Bowtie2 optimal: 8 threads per job (memory-constrained, ~2.5GB/thread)
MAX_PARALLEL_SAMPLES="${MAX_PARALLEL_SAMPLES:-2}"
if [[ "${_JOBS_MODE:-}" == "auto" || "${_JOBS_MODE:-}" == "AUTO" ]]; then
	MAX_PARALLEL_SAMPLES=$(( THREADS / 8 ))
	(( MAX_PARALLEL_SAMPLES < 1 )) && MAX_PARALLEL_SAMPLES=1
fi

# Threads allocated per RSEM job (auto-calculated from THREADS / MAX_PARALLEL_SAMPLES)
# Memory-aware: RSEM uses ~2-3GB per thread; cap to prevent OOM on constrained systems.
# Deferred RSEM thread calculation — only computed when first needed via _ensure_rsem_threads()
_RSEM_THREADS_COMPUTED=false
_ensure_rsem_threads() {
	[[ "$_RSEM_THREADS_COMPUTED" == "true" ]] && return 0
	_RSEM_THREADS_COMPUTED=true
	if [[ -z "${THREADS_PER_RSEM_JOB:-}" ]]; then
		THREADS_PER_RSEM_JOB=$((THREADS / MAX_PARALLEL_SAMPLES))
		# Memory guard: estimate available RAM and cap threads so total < 75% of RAM
		# RSEM uses ~2.5GB per thread on average
		local _rsem_avail_mb _rsem_usable_mb _rsem_max_threads_per_job
		_rsem_avail_mb=$(_get_available_ram_mb 2>/dev/null || echo 16384)
		_rsem_usable_mb=$(( _rsem_avail_mb * 75 / 100 ))
		_rsem_max_threads_per_job=$(( _rsem_usable_mb / 2560 / MAX_PARALLEL_SAMPLES ))
		[[ $_rsem_max_threads_per_job -lt 1 ]] && _rsem_max_threads_per_job=1
		[[ $THREADS_PER_RSEM_JOB -gt $_rsem_max_threads_per_job ]] && THREADS_PER_RSEM_JOB=$_rsem_max_threads_per_job
	fi
	[[ $THREADS_PER_RSEM_JOB -lt 1 ]] && THREADS_PER_RSEM_JOB=1
	# Safety cap: never exceed total THREADS (edge case with MAX_PARALLEL_SAMPLES=1 and low RAM)
	[[ $THREADS_PER_RSEM_JOB -gt $THREADS ]] && THREADS_PER_RSEM_JOB=$THREADS
}

# Cache tool availability at module load — avoids command -v subprocess per call
_M5_HAS_SALMON=false; command -v salmon &>/dev/null && _M5_HAS_SALMON=true
_M5_HAS_DOS2UNIX=false; command -v dos2unix &>/dev/null && _M5_HAS_DOS2UNIX=true

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

	# Return cached result if available (must be non-empty and a valid strandedness value)
	if [[ -s "$cache_file" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		RSEM_STRANDEDNESS=$(<"$cache_file")
		case "$RSEM_STRANDEDNESS" in
			none|forward|reverse)
				log_info "[STRANDEDNESS] Using cached result: $RSEM_STRANDEDNESS (from $cache_file)"
				return 0
				;;
			*)
				log_warn "[STRANDEDNESS] Cache file contains invalid value '$RSEM_STRANDEDNESS', re-detecting"
				RSEM_STRANDEDNESS="auto"
				;;
		esac
	fi

	# Check that Salmon is available (uses module-level cache to avoid per-call subprocess)
	if [[ "${_M5_HAS_SALMON:-false}" != "true" ]]; then
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
	tmp_dir=$(mktemp -d "${TMPDIR:-${index_root}}/strandedness_detect_XXXXXX") || {
		log_warn "[STRANDEDNESS] Failed to create temp directory — falling back to 'none'"
		RSEM_STRANDEDNESS="none"
		return 0
	}
	local salmon_idx="$tmp_dir/salmon_idx"
	local salmon_quant="$tmp_dir/salmon_quant"

	# Build a lightweight Salmon index (smaller k-mer for speed)
	local _detect_threads=$(( ${THREADS:-8} / 4 ))
	[[ $_detect_threads -lt 2 ]] && _detect_threads=2
	salmon index -t "$fasta" -i "$salmon_idx" --threads "$_detect_threads" -k 23 --keepDuplicates 2>"$tmp_dir/salmon_index.log" || {
		log_warn "[STRANDEDNESS] Salmon index failed — falling back to 'none'"
		rm -rf "$tmp_dir"
		RSEM_STRANDEDNESS="none"
		return 0
	}

	# Subsample reads for faster strandedness detection (~200k reads is sufficient)
	# This avoids mapping the entire FASTQ just to detect library type
	# Write uncompressed temp FASTQs — eliminates decompress→head→recompress→decompress cycle
	# (Salmon accepts uncompressed FASTQ, saves 2 compression rounds per read file)
	local _sub1="$tmp_dir/sub_R1.fq" _sub2=""
	local _subsample_lines=800000  # 200k reads × 4 lines per FASTQ record

	# Launch R1 subsampling (always needed)
	if [[ "$trimmed1" == *.gz ]]; then
		${_SHARED_GZIP_DC:-gzip -dc} "$trimmed1" | head -n $_subsample_lines > "$_sub1" &
	else
		head -n $_subsample_lines "$trimmed1" > "$_sub1" &
	fi
	local _sub1_pid=$!

	# Launch R2 subsampling in parallel (paired-end only) — 2x speedup vs sequential
	if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
		_sub2="$tmp_dir/sub_R2.fq"
		if [[ "$trimmed2" == *.gz ]]; then
			${_SHARED_GZIP_DC:-gzip -dc} "$trimmed2" | head -n $_subsample_lines > "$_sub2" &
		else
			head -n $_subsample_lines "$trimmed2" > "$_sub2" &
		fi
		local _sub2_pid=$!
	fi

	# Wait for both subsampling tasks to complete
	wait "$_sub1_pid" || true
	[[ -n "${_sub2_pid:-}" ]] && { wait "$_sub2_pid" || true; }

	# Quantify with --libType A (auto-detect) using subsampled reads + --skipQuant for speed
	local salmon_exit
	if [[ -n "$_sub2" && -f "$_sub2" ]]; then
		salmon quant -i "$salmon_idx" -l A \
			-1 "$_sub1" -2 "$_sub2" \
			-o "$salmon_quant" --threads "$_detect_threads" \
			--skipQuant 2>"$tmp_dir/salmon_quant.log"
		salmon_exit=$?
	else
		salmon quant -i "$salmon_idx" -l A \
			-r "$_sub1" \
			-o "$salmon_quant" --threads "$_detect_threads" \
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

	# Single-pass awk: extract both expected_format and fragment stats (saves 1 awk fork)
	local inferred_type num_compat _awk_out
	_awk_out=$(awk '
		/"expected_format"/ { split($0, a, "\""); fmt = a[4] }
		/compatible_fragment_ratio|num_compatible_fragments|num_assigned_fragments|strand|"read/ {
			gsub(/[{}":,]/, " "); gsub(/^[ \t]+|[ \t]+$/, "")
			n = split($0, f, /[ \t]+/)
			if (n >= 2) stats = stats "  " f[1] ": " f[2] "\n"
		}
		END { printf "%s\n---\n%s", (fmt ? fmt : "U"), stats }
	' "$lib_format" 2>/dev/null)
	inferred_type="${_awk_out%%$'\n'*}"
	num_compat="${_awk_out#*---$'\n'}"
	[[ "$num_compat" == "$_awk_out" ]] && num_compat=""

	# Map Salmon library type codes to RSEM strandedness
	# Salmon paired-end: IU=unstranded, ISF=forward(sense), ISR=reverse(antisense)
	# Salmon single-end: U=unstranded, SF=forward, SR=reverse
	case "$inferred_type" in
		IU|U)   RSEM_STRANDEDNESS="none" ;;
		ISF|SF) RSEM_STRANDEDNESS="forward" ;;
		ISR|SR) RSEM_STRANDEDNESS="reverse" ;;
		*)
			log_warn "[STRANDEDNESS] Unrecognized Salmon library type '${inferred_type:-}' — falling back to 'none'"
			RSEM_STRANDEDNESS="none"
			;;
	esac

	if [[ -n "${num_compat:-}" ]]; then
		log_info "[STRANDEDNESS] Salmon fragment stats:"
		log_info "$num_compat"
	fi

	log_info "[STRANDEDNESS] Detected: $RSEM_STRANDEDNESS (Salmon inferred: $inferred_type)"

	# Cache the result (atomic write via tmp+mv to avoid partial reads under concurrent runs)
	mkdir -p "$index_root"
	local _cache_tmp="${cache_file}.${BASHPID:-$$}"
	echo "$RSEM_STRANDEDNESS" > "$_cache_tmp" && mv -f "$_cache_tmp" "$cache_file"

	# Cleanup temp files
	rm -rf "$tmp_dir"

	return 0
}

# ==============================================================================
# MAIN PIPELINE: BOWTIE2 + RSEM
# ==============================================================================

bowtie2_rsem_pipeline() {
	_ensure_rsem_threads
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
	[[ ! -f "$fasta" ]] && { log_error "FASTA file not found: $fasta"; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	[[ ${#rnaseq_list[@]} -eq 0 ]] && { log_error "No RNA-seq samples provided."; return 1; }

	# Derive tag from original FASTA name for consistent output dir naming
	local _fn="${fasta##*/}"; local tag="${_fn%.*}"
	set_fasta_output_dirs "$tag"
	local rsem_idx="$RSEM_INDEX_ROOT/rsem_ref"
	local quant_root="$RSEM_QUANT_ROOT"
	local matrix_dir="$RSEM_MATRIX_ROOT"

	mkdir -p "$RSEM_INDEX_ROOT" "$quant_root" "$matrix_dir"

	# Convert line endings only if CRLF detected in first 100KB (avoids full-file scan
	# on multi-GB FASTA references; CRLF is always present in the header if present at all).
	# Write to a method-local copy to avoid modifying the shared reference FASTA
	# (which could race with concurrent M3/M4 runs or invalidate cache sentinels).
	local fasta_use="$fasta"
	if [[ "$_M5_HAS_DOS2UNIX" == "true" ]] && head -c 102400 "$fasta" 2>/dev/null | grep -q $'\r'; then
		local _m5_clean="${RSEM_INDEX_ROOT}/${tag}_nodoscr.fa"
		if [[ ! -f "$_m5_clean" || "$fasta" -nt "$_m5_clean" ]]; then
			dos2unix < "$fasta" > "$_m5_clean" 2>/dev/null || { log_warn "dos2unix copy failed, using original FASTA"; _m5_clean="$fasta"; }
		fi
		fasta_use="$_m5_clean"
	fi

	# Create gene-transcript mapping BEFORE building reference (needed for proper
	# gene-level aggregation when using transcript FASTAs with multiple isoforms)
	local gene_trans_map="${fasta_use}.gene_trans_map"
	if [[ ! -s "$gene_trans_map" ]]; then
		create_gene_trans_map "$fasta_use" "$gene_trans_map"
	fi

	# AUTO-DETECT STRANDEDNESS (runs once, caches result)
	if [[ "$RSEM_STRANDEDNESS" == "auto" ]]; then
		_rsem_detect_strandedness "$fasta_use" "${rnaseq_list[0]}" "$RSEM_INDEX_ROOT"
	fi

	# BUILD RSEM REFERENCE
	# Pass --transcript-to-gene-map when available so RSEM performs proper gene-level
	# aggregation in .genes.results (critical for transcript FASTAs with multiple isoforms)
	local _gtm_args=()
	if [[ -s "$gene_trans_map" ]]; then
		_gtm_args=(--transcript-to-gene-map "$gene_trans_map")
	fi
	# Check all required index files: .grp (RSEM) + .rev.2.bt2 (written last by Bowtie2 build)
	if [[ -f "${rsem_idx}.grp" && ( -f "${rsem_idx}.rev.2.bt2" || -f "${rsem_idx}.rev.2.bt2l" ) && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[RSEM INDEX] RSEM reference already exists. Skipping."
		if [[ ${#_gtm_args[@]} -gt 0 && ! -f "$RSEM_INDEX_ROOT/.has_gene_trans_map" ]]; then
			log_warn "[RSEM INDEX] Existing reference may lack gene-transcript mapping."
			log_warn "[RSEM INDEX] Set OVERWRITE_MODE=overwrite to rebuild with proper gene-level aggregation."
		fi
	else
		log_step "Building RSEM reference for $tag"
		log_file_size "$fasta_use" "Input FASTA for RSEM index - $tag"
		run_with_space_time_log --input "$fasta_use" --output "$RSEM_INDEX_ROOT" \
			rsem-prepare-reference --bowtie2 -p "$THREADS" "${_gtm_args[@]}" "$fasta_use" "$rsem_idx"
		log_file_size "$RSEM_INDEX_ROOT" "RSEM index output - $tag"
		# Mark that this reference was built with gene-transcript mapping
		[[ ${#_gtm_args[@]} -gt 0 ]] && touch "$RSEM_INDEX_ROOT/.has_gene_trans_map"
	fi

	# QUANTIFY SAMPLES - parallel or sequential
	if _rsem_should_use_parallel && [[ ${#rnaseq_list[@]} -gt 1 ]]; then
		_rsem_quantify_parallel "$rsem_idx" "$quant_root" rnaseq_list[@] || \
			log_error "[RSEM] Some parallel quantification jobs failed — check logs before relying on matrix output"
	else
		_rsem_quantify_sequential "$rsem_idx" "$quant_root" rnaseq_list[@] || \
			log_error "[RSEM] Some sequential quantification jobs failed — check logs before relying on matrix output"
	fi

	# GENERATE MATRICES
	_create_rsem_matrices "$fasta_use" "$tag" "$quant_root" "$matrix_dir" rnaseq_list[@] || \
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
	if [[ $failed_samples -gt 0 ]]; then
		log_warn "[RSEM QUANT] $failed_samples/${#samples[@]} sample(s) failed"
		return 1
	fi
	return 0
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
		printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || ts=$(date '+%Y-%m-%d %H:%M:%S')
		echo "[$ts] [$level] [RSEM-$SRR] $*"
		if [[ "$level" != "INFO" && -n "${abs_error_warn_file:-}" ]]; then
			# Use per-method error log when _PARALLEL_METHOD_ID is set (avoids write contention)
			local _ew="${abs_error_warn_file}"
			[[ -n "${_PARALLEL_METHOD_ID:-}" ]] && _ew="${abs_error_warn_file%.log}_${_PARALLEL_METHOD_ID}.log"
			echo "[$ts] [$level] [RSEM-$SRR] $*" >> "$_ew"
		fi
	}

	[[ -z "$SRR" ]] && { _plog "ERROR" "Empty SRR ID - skipping"; return 1; }

	# Reactivate conda in subshell if needed (skip when orchestrator manages env)
	if [[ -z "${WF_MANAGED_ENV:-}" && -n "${CONDA_PREFIX:-}" && -n "${CONDA_EXE:-}" ]]; then
		source "${_CONDA_PROFILE_SCRIPT:-${CONDA_EXE%/*}/../etc/profile.d/conda.sh}" 2>/dev/null || true
		conda activate "${CONDA_DEFAULT_ENV:-base}" 2>/dev/null || true
	fi

	local out_dir="$quant_root/$SRR"
	mkdir -p "$out_dir"

	# Skip if already processed
	if [[ -f "$out_dir/${SRR}.genes.results" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		_plog "INFO" "Results already exist. Skipping."
		return 0
	fi

	# Locate trimmed reads via shared helper (exported by shared_utils_method.sh;
	# sets trimmed1 and trimmed2 — also handles conda reactivation, so the earlier
	# conda block above is a no-op safety net for edge cases)
	_init_parallel_worker "$SRR"

	if [[ -z "$trimmed1" ]]; then
		local _td="$abs_trim_dir_root/$SRR"
		_plog "ERROR" "Missing trimmed reads in: $_td"
		# Guard ls subshell — only spawn if directory exists (avoids unconditional fork)
		if [[ -d "$_td" ]]; then
			_plog "ERROR" "Contents: $(ls -la "$_td" 2>&1)"
		else
			_plog "ERROR" "Directory does not exist: $_td"
		fi
		return 1
	fi

	if [[ ! -f "${rsem_idx}.grp" || ( ! -f "${rsem_idx}.rev.2.bt2" && ! -f "${rsem_idx}.rev.2.bt2l" ) ]]; then
		_plog "ERROR" "RSEM index incomplete or not found: ${rsem_idx}.grp / .rev.2.bt2[l]"
		return 1
	fi

	_plog "INFO" "Processing with $threads_per_job threads (strandedness: ${RSEM_STRANDEDNESS:-none})"

	# Skip BAM output when BAMs won't be kept (saves significant disk I/O)
	local _no_bam_flag=""
	[[ "${keep_bam_global:-n}" != "y" ]] && _no_bam_flag="--no-bam-output"

	local rsem_log="$out_dir/${SRR}.rsem.log"
	local rsem_exit_code
	# Redirect directly to log file — parallel already captures stdout from this worker,
	# so tee would double the I/O. Direct redirect also simplifies exit code capture.
	if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
		rsem-calculate-expression \
			--paired-end \
			--bowtie2 \
			--bowtie2-sensitivity-level "${BOWTIE2_MODE:-sensitive}" \
			--strandedness "${RSEM_STRANDEDNESS:-none}" \
			--seed "$RSEM_SEED" \
			--num-threads "$threads_per_job" \
			${_no_bam_flag:+"$_no_bam_flag"} \
			"$trimmed1" "$trimmed2" "$rsem_idx" "$out_dir/$SRR" > "$rsem_log" 2>&1
		rsem_exit_code=$?
	else
		rsem-calculate-expression \
			--bowtie2 \
			--bowtie2-sensitivity-level "${BOWTIE2_MODE:-sensitive}" \
			--strandedness "${RSEM_STRANDEDNESS:-none}" \
			--seed "$RSEM_SEED" \
			--num-threads "$threads_per_job" \
			${_no_bam_flag:+"$_no_bam_flag"} \
			"$trimmed1" "$rsem_idx" "$out_dir/$SRR" > "$rsem_log" 2>&1
		rsem_exit_code=$?
	fi

	if [[ $rsem_exit_code -ne 0 ]]; then
		_plog "ERROR" "RSEM failed (exit code: $rsem_exit_code)"
		# Batch log output: capture tail once, log as single block.
		# Avoids O(L) subshell forks from while-read pipeline (L = lines).
		if [[ -f "$rsem_log" ]]; then
			local _tail_buf; _tail_buf=$(tail -20 "$rsem_log" 2>/dev/null) || true
			[[ -n "$_tail_buf" ]] && _plog "ERROR" "$_tail_buf"
		fi
		return $rsem_exit_code
	fi

	if [[ ! -f "$out_dir/${SRR}.genes.results" ]]; then
		_plog "ERROR" "RSEM completed but output missing: $out_dir/${SRR}.genes.results"
		if [[ -d "$out_dir" ]]; then
			_plog "ERROR" "  Contents: $(ls -la "$out_dir" 2>&1)"
		else
			_plog "ERROR" "  Output directory does not exist: $out_dir"
		fi
		return 1
	fi

	# Cleanup any residual BAM files (safety net — --no-bam-output should prevent creation)
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

	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	# Adaptive: cap parallel jobs at sample count to maximize per-job thread allocation
	(( parallel_jobs > num_samples )) && parallel_jobs=$num_samples
	(( parallel_jobs < 1 )) && parallel_jobs=1
	local threads_per_job="${THREADS_PER_RSEM_JOB:-$((THREADS / parallel_jobs))}"
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	log_step "[RSEM QUANT] Running PARALLEL quantification: $num_samples samples, $parallel_jobs concurrent jobs, $threads_per_job threads/job"
	log_info "[RSEM QUANT] Sample list: ${valid_samples[*]}"

	# Set up environment for parallel subshells using shared utility
	_prepare_parallel_env "M5"
	export rsem_idx quant_root threads_per_job OVERWRITE_MODE BOWTIE2_MODE RSEM_STRANDEDNESS RSEM_SEED

	# O(S/parallel_jobs × (N log N + T)) — S samples batched across parallel_jobs slots;
	# each worker runs Bowtie2 O(N log N) + RSEM EM O(T × iterations)
	parallel \
		--env PATH \
		--env CONDA_PREFIX \
		--env CONDA_DEFAULT_ENV \
		--env CONDA_EXE \
		--env _CONDA_PROFILE_SCRIPT \
		--env WF_MANAGED_ENV \
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
		--halt soon,fail,1 \
		--joblog "$quant_root/parallel_rsem_${BASHPID:-$$}.log" \
		--progress \
		_rsem_parallel_worker {} \
		< <(printf "%s\n" "${valid_samples[@]}")

	local parallel_exit=$?

	# Report results — count only files for the current batch (avoids inflation from stale runs)
	local successful=0
	for s in "${valid_samples[@]}"; do
		[[ -f "$quant_root/$s/${s}.genes.results" ]] && successful=$((successful + 1))
	done
	log_info "[RSEM QUANT] Parallel quantification complete: $successful/$num_samples samples succeeded"

	if [[ -f "$quant_root/parallel_rsem_${BASHPID:-$$}.log" ]]; then
		local failed
		failed=$(awk 'NR>1 && $7!=0 {n++} END{print n+0}' "$quant_root/parallel_rsem_${BASHPID:-$$}.log")
		[[ $failed -gt 0 ]] && log_warn "[RSEM QUANT] $failed sample(s) failed - check $quant_root/parallel_rsem_${BASHPID:-$$}.log"
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
	# Cache availability probe (avoid per-call PATH scan, consistent with _HAS_PREPDE in m1)
	if [[ -z "${_M5_HAS_ABUND_MATRIX:-}" ]]; then
		_M5_HAS_ABUND_MATRIX=false; command -v abundance_estimates_to_matrix.pl >/dev/null 2>&1 && _M5_HAS_ABUND_MATRIX=true
	fi
	if [[ "$_M5_HAS_ABUND_MATRIX" == "true" ]]; then
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
			awk -F'\t' 'NR>1 {print $1}' "$quant_root/$SRR/${SRR}.genes.results" > "$temp_gene_ids"
			break
		fi
	done

	if [[ -z "$first_sample" ]]; then
		log_warn "[RSEM MATRIX] No samples have RSEM results — cannot create manual matrix"
		return 1
	fi

	# O(1) fork — wc -l is a single C-level scan, faster than O(G) bash read loop
	local num_genes; num_genes=$(wc -l < "$temp_gene_ids")

	for SRR in "${srr_list[@]}"; do
		if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
			# Single awk pass: extract counts/TPM/FPKM and zero-fill if gene count mismatches.
			# Handles both normal case and mismatch in one process (was 2 awk invocations on mismatch).
			# O(G) per sample, single file traversal.
			local sample_genes
			sample_genes=$(awk -F'\t' -v n="$num_genes" \
				-v c="$matrix_dir/${SRR}_counts.tmp" \
				-v t="$matrix_dir/${SRR}_tpm.tmp" \
				-v f="$matrix_dir/${SRR}_fpkm.tmp" '
				NR>1 { print $5 > c; print $6 > t; print $7 > f; g++ }
				END  {
					if (g < n) { for (i = g; i < n; i++) { print 0 > c; print 0 > t; print 0 > f } }
					print g+0
				}
			' "$quant_root/$SRR/${SRR}.genes.results")
			if [[ "$sample_genes" -ne "$num_genes" ]]; then
				log_warn "[RSEM MATRIX] $SRR has $sample_genes genes (expected $num_genes) — zero-filled to match"
			fi
		else
			log_warn "[RSEM MATRIX] Missing results for $SRR — filling with zeros in count matrix"
			# Single awk generates all three zero-fill files (replaces 3 yes|head pipelines = 6 processes)
			awk -v n="$num_genes" -v c="$matrix_dir/${SRR}_counts.tmp" \
				-v t="$matrix_dir/${SRR}_tpm.tmp" -v f="$matrix_dir/${SRR}_fpkm.tmp" \
				'BEGIN{for(i=0;i<n;i++){print 0>c; print 0>t; print 0>f}}'
		fi
	done

	# Build ordered file lists matching srr_list order (avoids glob sort mismatch)
	local count_files=() tpm_files=() fpkm_files=()
	for SRR in "${srr_list[@]}"; do
		count_files+=("$matrix_dir/${SRR}_counts.tmp")
		tpm_files+=("$matrix_dir/${SRR}_tpm.tmp")
		fpkm_files+=("$matrix_dir/${SRR}_fpkm.tmp")
	done

	# Build headers with single printf (was N+1 echo calls per matrix, 3 matrices)
	local header; printf -v header '\t%s' "${srr_list[@]}"

	# Build all 3 matrices concurrently — each paste reads temp_gene_ids once in parallel
	# (background jobs share OS page cache so temp_gene_ids is only loaded from disk once)
	# O(G × S) per matrix — paste joins G gene-id rows across S per-sample column files; 3 matrices run concurrently
	{ printf 'gene_id%s\n' "$header"; paste "$temp_gene_ids" "${count_files[@]}"; } > "$matrix_dir/genes.counts.matrix" &
	local _pid_counts=$!
	{ printf 'gene_id%s\n' "$header"; paste "$temp_gene_ids" "${tpm_files[@]}"; } > "$matrix_dir/genes.TPM.not_cross_norm" &
	local _pid_tpm=$!
	{ printf 'gene_id%s\n' "$header"; paste "$temp_gene_ids" "${fpkm_files[@]}"; } > "$matrix_dir/genes.FPKM.not_cross_norm" &
	local _pid_fpkm=$!
	local _matrix_fail=0
	wait "$_pid_counts" || { log_warn "[RSEM MATRIX] Count matrix assembly failed"; _matrix_fail=1; }
	wait "$_pid_tpm"   || { log_warn "[RSEM MATRIX] TPM matrix assembly failed"; _matrix_fail=1; }
	wait "$_pid_fpkm"  || { log_warn "[RSEM MATRIX] FPKM matrix assembly failed"; _matrix_fail=1; }

	rm -f "$temp_gene_ids" "$matrix_dir"/*_counts.tmp "$matrix_dir"/*_tpm.tmp "$matrix_dir"/*_fpkm.tmp

	return $_matrix_fail
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
	# All three sed conversions are independent — run concurrently as background jobs
	# to reduce wall-clock from O(3×file_size) to O(max(file_size)) for the conversion step.
	local _csv_pids=()
	if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
		if [[ ! -f "$gene_count_matrix" || "$matrix_dir/genes.counts.matrix" -nt "$gene_count_matrix" ]]; then
			log_info "[RSEM MATRIX] Converting count matrix to CSV format..."
			sed '1{s/^gene_id\t/Gene_ID\t/; s/^\t/Gene_ID\t/}; s/\t/,/g' "$matrix_dir/genes.counts.matrix" > "$gene_count_matrix" &
			_csv_pids+=($!)
		fi
	fi

	# Create TPM and FPKM matrices (always regenerate if source is newer)
	if [[ -f "$matrix_dir/genes.TPM.not_cross_norm" ]]; then
		local tpm_matrix="$deseq2_dir/gene_tpm_matrix.csv"
		if [[ ! -f "$tpm_matrix" || "$matrix_dir/genes.TPM.not_cross_norm" -nt "$tpm_matrix" ]]; then
			sed '1{s/^gene_id\t/Gene_ID\t/; s/^\t/Gene_ID\t/}; s/\t/,/g' "$matrix_dir/genes.TPM.not_cross_norm" > "$tpm_matrix" &
			_csv_pids+=($!)
		fi
	fi

	if [[ -f "$matrix_dir/genes.FPKM.not_cross_norm" ]]; then
		local fpkm_matrix="$deseq2_dir/gene_fpkm_matrix.csv"
		if [[ ! -f "$fpkm_matrix" || "$matrix_dir/genes.FPKM.not_cross_norm" -nt "$fpkm_matrix" ]]; then
			sed '1{s/^gene_id\t/Gene_ID\t/; s/^\t/Gene_ID\t/}; s/\t/,/g' "$matrix_dir/genes.FPKM.not_cross_norm" > "$fpkm_matrix" &
			_csv_pids+=($!)
		fi
	fi

	# Wait for all CSV conversions to complete before proceeding
	local _csv_pid _csv_fail=0
	for _csv_pid in "${_csv_pids[@]+${_csv_pids[@]}}"; do
		wait "$_csv_pid" 2>/dev/null || { log_warn "[RSEM MATRIX] CSV conversion job (PID=$_csv_pid) failed"; _csv_fail=1; }
	done
	if (( _csv_fail )); then
		log_error "[RSEM MATRIX] One or more CSV conversions failed in $deseq2_dir"
	fi

	# Create sample metadata (regenerate if missing or samples changed)
	create_sample_metadata "$sample_metadata" "${srr_list[@]}"

	# Generate tximport script (regenerate if missing)
	local tximport_script="$deseq2_dir/run_tximport_rsem.R"
	if [[ ! -f "$tximport_script" || "${OVERWRITE_MODE:-skip}" == "overwrite" ]]; then
		generate_tximport_script "rsem" "$quant_root" "$tximport_script" "$sample_metadata"
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
		# O(1) bash builtin — avoids $(date) subprocess fork
		printf -v _rsem_ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _rsem_ts=$(date '+%Y-%m-%d %H:%M:%S')
		echo "Date: $_rsem_ts"
		echo "Samples processed: ${#srr_list[@]}"
		echo "Method: RSEM with Bowtie2 alignment (strandedness: ${RSEM_STRANDEDNESS:-none}, sensitivity: ${BOWTIE2_MODE:-sensitive})"
		echo ""
		echo "Per-sample statistics:"
		echo "-------------------------------------------------------------------"

		for SRR in "${srr_list[@]}"; do
			if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
				# Single AWK pass: count total genes, expressed genes, and sum counts
				local total expressed counts
				read -r total expressed counts < <(awk -F'\t' 'NR>1 { total++; if ($5>0) expr++; s+=$5 }
					END { printf "%d %d %d", total+0, expr+0, int(s) }
				' "$quant_root/$SRR/${SRR}.genes.results")

				# Extract Bowtie2 alignment rate from RSEM log
				local align_rate="N/A"
				local rsem_log="$quant_root/$SRR/${SRR}.rsem.log"
				if [[ -f "$rsem_log" ]]; then
					align_rate=$(awk '/% overall alignment rate/{match($0,/[0-9]+\.[0-9]+/);r=substr($0,RSTART,RLENGTH)} END{print r}' "$rsem_log")
					if [[ -z "$align_rate" ]]; then
						align_rate="N/A"
					# Bash integer comparison avoids awk fork per sample. O(S) forks saved.
					elif (( ${align_rate%%.*} < LOW_ALIGN_THRESHOLD )); then
						outlier_samples+=("$SRR ($align_rate%)")
					fi
				fi

				local _rate_display="${align_rate}"
				[[ "$align_rate" != "N/A" ]] && _rate_display="${align_rate}%"
				echo "$SRR: $expressed/$total expressed genes, $counts expected counts, alignment rate: ${_rate_display}"
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

# Check if GNU Parallel should be used for sample processing.
# Consistent with M1-M4 inline guards: default TRUE (opt-out, not opt-in),
# parallel binary available, AND JOBS > 1 (single-job parallel is wasteful).
_rsem_should_use_parallel() {
	[[ "${USE_GNU_PARALLEL:-TRUE}" == "FALSE" ]] && return 1
	$_SHARED_HAS_PARALLEL || return 1
	[[ "${JOBS:-1}" -gt 1 ]] || return 1
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
		if [[ -d "$TRIM_DIR_ROOT/$SRR" ]]; then
			log_warn "  Contents: $(ls -la "$TRIM_DIR_ROOT/$SRR" 2>&1)"
		else
			log_warn "  Directory does not exist: $TRIM_DIR_ROOT/$SRR"
		fi
		return 1
	fi

	# Skip BAM output when BAMs won't be kept (saves significant disk I/O)
	local _no_bam_flag=""
	[[ "${keep_bam_global:-n}" != "y" ]] && _no_bam_flag="--no-bam-output"

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
			${_no_bam_flag:+"$_no_bam_flag"} \
			"$trimmed1" "$trimmed2" "$rsem_idx" "$out_dir/$SRR" > "$rsem_log" 2>&1
		rsem_exit_code=$?
	else
		log_step "Running Bowtie2 + RSEM (single-end) for $SRR (threads: $threads_to_use, strandedness: ${RSEM_STRANDEDNESS:-none})"
		rsem-calculate-expression \
			--bowtie2 \
			--bowtie2-sensitivity-level "${BOWTIE2_MODE:-sensitive}" \
			--strandedness "${RSEM_STRANDEDNESS:-none}" \
			--seed "$RSEM_SEED" \
			--num-threads "$threads_to_use" \
			${_no_bam_flag:+"$_no_bam_flag"} \
			"$trimmed1" "$rsem_idx" "$out_dir/$SRR" > "$rsem_log" 2>&1
		rsem_exit_code=$?
	fi

	if [[ $rsem_exit_code -ne 0 ]]; then
		log_error "[RSEM QUANT] RSEM failed for $SRR (exit code: $rsem_exit_code)"
		log_error "[RSEM QUANT] Check log: $rsem_log"
		# Single log call with tail output (avoids per-line subshell)
		[[ -f "$rsem_log" ]] && log_error "$(tail -20 "$rsem_log" 2>/dev/null)"
		return $rsem_exit_code
	fi

	# Verify output exists even when RSEM exits 0 (edge case: disk full, interrupted write)
	if [[ ! -f "$out_dir/${SRR}.genes.results" ]]; then
		log_error "[RSEM QUANT] RSEM completed but output missing: $out_dir/${SRR}.genes.results"
		if [[ -d "$out_dir" ]]; then
			log_error "  Contents: $(ls -la "$out_dir" 2>&1)"
		else
			log_error "  Output directory does not exist: $out_dir"
		fi
		return 1
	fi

	log_file_size "$out_dir/${SRR}.genes.results" "RSEM gene results - $SRR"

	# Cleanup any residual BAM files (safety net — --no-bam-output should prevent creation)
	if [[ "${keep_bam_global:-n}" != "y" ]]; then
		log_info "[CLEANUP] Removing any residual RSEM BAM files for $SRR"
		rm -f "$out_dir/${SRR}.transcript.bam" "$out_dir/${SRR}.genome.bam" \
			  "$out_dir/${SRR}.transcript.sorted.bam" "$out_dir/${SRR}.transcript.sorted.bam.bai"
	fi

	return 0
}
