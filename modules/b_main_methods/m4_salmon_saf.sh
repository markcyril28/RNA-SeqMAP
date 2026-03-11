#!/bin/bash
# ==============================================================================
# METHOD 4: SALMON SAF (SELECTIVE ALIGNMENT WITH DECOYS) PIPELINE
# ==============================================================================
# Quantify expression using decoy-aware Salmon (Selective Alignment)
# Fast and accurate pseudo-alignment
# ==============================================================================

# NOTE: set -e/-u/-o pipefail are intentionally NOT set here — this file is sourced
# as a library and must not alter the parent shell's error-exit behaviour.
# Pipe failures in salmon quant are checked manually via PIPESTATUS (lines ~119, ~127).

# Guard against double-sourcing
[[ "${M4_SALMON_SOURCED:-}" == "true" ]] && return 0
export M4_SALMON_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# SALMON CONFIGURATION - IMPORTANT PARAMETERS (tweak here)
# ==============================================================================

# k-mer size for Salmon index (31 is standard; reduce to 21 for very short reads)
SALMON_KMER_SIZE="${SALMON_KMER_SIZE:-31}"
# Number of bootstraps for uncertainty estimation (0 = off; 100 for sleuth/Swish)
SALMON_NUM_BOOTSTRAPS="${SALMON_NUM_BOOTSTRAPS:-0}"

# ==============================================================================
# SALMON SAF PIPELINE
# ==============================================================================

salmon_saf_pipeline() {
	local fasta="" genome="" rnaseq_list=()
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--GENOME) genome="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do rnaseq_list+=("$1"); shift; done;;
			*) log_error "Unknown arg: $1"; return 1;;
		esac
	done

	[[ -z "$fasta" || -z "$genome" ]] && { log_error "Usage: --FASTA genes.fa --GENOME genome.fa"; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")

	local tag="$(basename "${fasta%.*}")"
	set_fasta_output_dirs "$tag"
	# Use an absolute path so $work is unambiguous regardless of the caller's CWD
	local work="$SALMON_SAF_ROOT/tmp_${tag}_gentrome"
	local idx_dir="$SALMON_INDEX_ROOT/decoySAF"
	local quant_root="$SALMON_QUANT_ROOT"
	local matrix_dir="$SALMON_MATRIX_ROOT"

	mkdir -p "$SALMON_INDEX_ROOT" "$quant_root" "$matrix_dir"

	# BUILD DECOY-AWARE INDEX
	if [[ -f "$idx_dir/versionInfo.json" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[SALMON INDEX] Decoy index already exists. Skipping."
	else
		# $work is created here only — not unconditionally — so it isn't left as an empty
		# orphan directory on runs where the index already exists.
		mkdir -p "$work"
		log_step "Building decoy-aware Salmon index for $tag"
		grep "^>" "$genome" | awk '{print substr($1,2)}' > "$work/decoys.txt"
		cat "$fasta" "$genome" > "$work/gentrome.fa"
		log_file_size "$work/gentrome.fa" "Gentrome FASTA for Salmon - $tag"
		log_file_size "$work/decoys.txt" "Decoy list for Salmon - $tag"
		run_with_space_time_log --input "$work" --output "$idx_dir" salmon index \
			-t "$work/gentrome.fa" \
			-d "$work/decoys.txt" \
			-i "$idx_dir" \
			-k "$SALMON_KMER_SIZE" -p "$THREADS"
		log_file_size "$idx_dir" "Salmon index output - $tag"
		# Only clean up the temporary gentrome directory if the index was built successfully.
		# If salmon index failed, versionInfo.json will be absent; keep $work for debugging.
		if [[ -f "$idx_dir/versionInfo.json" ]]; then
			log_info "[CLEANUP] Removing temporary gentrome work directory"
			rm -rf "$work"
		else
			log_warn "[INDEX] salmon index may have failed — keeping $work for inspection"
		fi
	fi

	# QUANTIFICATION PER SRR
	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]]; then
		log_step "[PARALLEL] Salmon quantification: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export idx_dir quant_root threads_per_job
		local salmon_num_bootstraps="$SALMON_NUM_BOOTSTRAPS"
		export salmon_num_bootstraps

		_m4_salmon_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log SALMON "$SRR" WARN "Missing trimmed reads - skipping"; return 0; }

			local out_dir="$quant_root/$SRR"
			mkdir -p "$out_dir"
			[[ -f "$out_dir/quant.sf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]] && { _parallel_log SALMON "$SRR" INFO "Already quantified - skipping"; return 0; }

			_parallel_log SALMON "$SRR" INFO "Quantifying with $threads_per_job threads"
			local salmon_exit=0
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				salmon quant -i "$idx_dir" -l A \
					-1 "$trimmed1" -2 "$trimmed2" \
					-p "$threads_per_job" \
					--numBootstraps "$salmon_num_bootstraps" \
					--gcBias --seqBias \
					-o "$out_dir" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				salmon_exit=${PIPESTATUS[0]}
			else
				salmon quant -i "$idx_dir" -l A \
					-r "$trimmed1" \
					-p "$threads_per_job" \
					--numBootstraps "$salmon_num_bootstraps" \
					--seqBias \
					-o "$out_dir" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				salmon_exit=${PIPESTATUS[0]}
			fi

			[[ $salmon_exit -ne 0 ]] && { _parallel_log SALMON "$SRR" ERROR "Salmon failed (exit=$salmon_exit)"; return $salmon_exit; }
			[[ ! -f "$out_dir/quant.sf" ]] && { _parallel_log SALMON "$SRR" ERROR "quant.sf not created"; return 1; }
			_parallel_log SALMON "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m4_salmon_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env idx_dir --env quant_root --env threads_per_job --env salmon_num_bootstraps \
			--env OVERWRITE_MODE \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$quant_root/parallel_salmon_saf.log" \
			_m4_salmon_parallel_worker {}

		local par_exit=$?
		# Count only samples from this run's rnaseq_list to avoid inflating the count
		# with stale quant.sf files from prior runs under quant_root.
		local successful=0
		for _s in "${rnaseq_list[@]}"; do
			[[ -f "$quant_root/$_s/quant.sf" ]] && successful=$((successful + 1))
		done
		log_info "[PARALLEL] Salmon SAF complete: $successful/${#rnaseq_list[@]} samples succeeded"
		[[ $par_exit -ne 0 ]] && log_warn "[PARALLEL] Some jobs failed - check $quant_root/parallel_salmon_saf.log"
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local out_dir="$quant_root/$SRR"
			mkdir -p "$out_dir"

			[[ -f "$out_dir/quant.sf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]] && { log_info "[SALMON QUANT] Quantification for $SRR already exists. Skipping."; continue; }

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "Missing trimmed reads for $SRR. Skipping."; continue; }

			log_step "Quantifying expression for $SRR with Salmon"

			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				log_info "[SALMON QUANT] Using paired-end reads for $SRR"
				run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$out_dir" salmon quant \
					-i "$idx_dir" -l A \
					-1 "$trimmed1" -2 "$trimmed2" \
					-p "$THREADS" \
					--numBootstraps "$SALMON_NUM_BOOTSTRAPS" \
					--gcBias --seqBias \
					-o "$out_dir"
			else
				log_info "[SALMON QUANT] Using single-end reads for $SRR"
				run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$out_dir" salmon quant \
					-i "$idx_dir" -l A \
					-r "$trimmed1" \
					-p "$THREADS" \
					--numBootstraps "$SALMON_NUM_BOOTSTRAPS" \
					--seqBias \
					-o "$out_dir"
			fi
			log_file_size "$out_dir/quant.sf" "Salmon quantification output - $SRR"
		done
	fi

	# MERGE MATRICES
	_create_salmon_matrices "$fasta" "$tag" "$quant_root" "$matrix_dir" rnaseq_list[@]

	log_step "COMPLETED: Salmon-SAF pipeline for $tag"
}

# ==============================================================================
# MATRIX GENERATION
# ==============================================================================

_create_salmon_matrices() {
	local fasta="$1"
	local tag="$2"
	local quant_root="$3"
	local matrix_dir="$4"
	local _arr_name="${5:-}"
	local srr_list=()
	if [[ -n "$_arr_name" ]]; then
		local _tmp=("${!_arr_name}")
		for _s in "${_tmp[@]}"; do [[ -n "$_s" ]] && srr_list+=("$_s"); done
	fi

	log_step "Generating gene and transcript matrices (Salmon)"

	# Check or create gene_trans_map
	local gene_trans_map="${fasta}.gene_trans_map"
	if [[ ! -f "$gene_trans_map" ]]; then
		log_info "[SALMON MATRIX] Creating gene-transcript mapping file..."
		create_gene_trans_map "$fasta" "$gene_trans_map"
	fi
	# Export so tximport_salmon_to_matrices.R can locate it via GENE_TRANS_MAP_FILE
	# without needing to do a recursive glob search across inputs/
	export GENE_TRANS_MAP_FILE="$gene_trans_map"

	# abundance_estimates_to_matrix.pl (Trinity) only supports RSEM|eXpress|kallisto,
	# not salmon — skip it and use the manual fallback directly.
	log_info "[SALMON MATRIX] abundance_estimates_to_matrix.pl does not support --est_method salmon; using manual count matrix builder."
	_create_manual_salmon_matrix "$quant_root" "$matrix_dir" srr_list[@]
}

_create_manual_salmon_matrix() {
	local quant_root="$1"
	local matrix_dir="$2"
	local _arr_name="${3:-}"
	local srr_list=()
	if [[ -n "$_arr_name" ]]; then
		local _tmp=("${!_arr_name}")
		for _s in "${_tmp[@]}"; do [[ -n "$_s" ]] && srr_list+=("$_s"); done
	fi

	# Fallback: discover samples from quant.sf files if SRR list is empty
	if [[ ${#srr_list[@]} -eq 0 ]]; then
		while IFS= read -r _qf; do
			srr_list+=("$(basename "$(dirname "$_qf")")");
		done < <(find "$quant_root" -name "quant.sf" 2>/dev/null)
	fi

	local temp_gene_ids="$matrix_dir/temp_gene_ids.txt"
	local temp_counts="$matrix_dir/temp_counts.txt"

	local first_sample=""
	for SRR in "${srr_list[@]}"; do
		if [[ -f "$quant_root/$SRR/quant.sf" ]]; then
			first_sample="$SRR"
			awk 'NR>1 {print $1}' "$quant_root/$SRR/quant.sf" > "$temp_gene_ids"
			break
		fi
	done

	if [[ -n "$first_sample" ]]; then
		for SRR in "${srr_list[@]}"; do
			if [[ -f "$quant_root/$SRR/quant.sf" ]]; then
				awk 'NR>1 {print int($5 + 0.5)}' "$quant_root/$SRR/quant.sf" > "$matrix_dir/${SRR}_counts.tmp"
			else
				# Separate declaration from assignment so wc errors are not masked by local
				local num_genes
				num_genes=$(wc -l < "$temp_gene_ids")
				yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_counts.tmp"
			fi
		done

		local _paste_args=("$temp_gene_ids")
		for _s in "${srr_list[@]}"; do
			[[ -f "$matrix_dir/${_s}_counts.tmp" ]] && _paste_args+=("$matrix_dir/${_s}_counts.tmp")
		done
		paste "${_paste_args[@]}" > "$temp_counts"

		# NOTE: Column 1 of quant.sf is the transcript Name, so this fallback matrix
		# is transcript-level despite the "genes" filename.  The authoritative gene-level
		# matrices are produced by tximport_salmon_to_matrices.R (uses tximport aggregation).
		# When abundance_estimates_to_matrix.pl is available it applies --gene_trans_map
		# to produce true gene-level output; this manual path does not.
		echo -n "transcript_id" > "$matrix_dir/genes.counts.matrix"
		for SRR in "${srr_list[@]}"; do
			echo -ne "\t$SRR" >> "$matrix_dir/genes.counts.matrix"
		done
		echo "" >> "$matrix_dir/genes.counts.matrix"
		cat "$temp_counts" >> "$matrix_dir/genes.counts.matrix"

		rm -f "$temp_gene_ids" "$temp_counts" "$matrix_dir"/*_counts.tmp
	fi
}
