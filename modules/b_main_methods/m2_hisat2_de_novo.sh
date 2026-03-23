#!/bin/bash
# ==============================================================================
# METHOD 2: HISAT2 DE NOVO PIPELINE
# ==============================================================================
# HISAT2 De Novo alignment with StringTie assembly
# Does not use reference GTF - discovers transcripts de novo
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M2_HISAT2_DN_SOURCED:-}" == "true" ]] && return 0
export M2_HISAT2_DN_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/b_main_methods}"
if [[ -z "$SCRIPT_DIR" ]]; then SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."; fi
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# HISAT2 DE NOVO PIPELINE
# ==============================================================================

hisat2_de_novo_pipeline() {
	local fasta="" rnaseq_list=()
	
	# Parse arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done
	
	# Validate inputs
	[[ -z "$fasta" ]] && { log_error "No FASTA file specified. Use --FASTA <fasta_file>."; return 1; }
	[[ ! -f "$fasta" ]] && { log_error "FASTA file '$fasta' not found."; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")

	local fasta_base fasta_tag index_prefix
	fasta_base="${fasta##*/}"
	fasta_tag="${fasta_base%.*}"
	set_fasta_output_dirs "$fasta_tag"
	index_prefix="$HISAT2_DE_NOVO_INDEX_DIR/${fasta_tag}_index"

	# BUILD HISAT2 INDEX
	mkdir -p "$HISAT2_DE_NOVO_INDEX_DIR"
	# Check for existing index via specific file (avoids ls glob subprocess + edge cases)
	if [[ -f "${index_prefix}.1.ht2" ]] && [[ "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[HISAT2 INDEX] De novo index exists - skipping build"
	else
		log_step "Building HISAT2 de novo index from $fasta"
		log_file_size "$fasta" "Input FASTA for HISAT2 de novo index"
		run_with_space_time_log --input "$fasta" --output "$HISAT2_DE_NOVO_INDEX_DIR" \
			hisat2-build -p "${THREADS}" "$fasta" "$index_prefix" \
			|| { log_error "HISAT2 de novo index build failed for $fasta_tag"; rm -f "${index_prefix}".*.ht2; return 1; }
		log_file_size "$HISAT2_DE_NOVO_INDEX_DIR" "HISAT2 de novo index output"
	fi

	# STRANDNESS AUTO-DETECTION
	# M2 de novo aligns to a transcriptome FASTA where all reference seqs are on
	# the + strand.  By aligning the first sample and checking read1 strand bias
	# we can infer the library prep protocol:
	#   RF (dUTP / TruSeq):  read1 maps predominantly reverse
	#   FR (ligation):       read1 maps predominantly forward
	#   Unstranded:          ~50/50 split
	# The result is cached in the index dir so resume runs skip re-detection.
	local hisat2_strand_opts="" stringtie_strand_opt="" _detected_strand=""
	local _strand_cache="$HISAT2_DE_NOVO_INDEX_DIR/detected_strandness.txt"

	if [[ -f "$_strand_cache" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		source "$_strand_cache"
		log_info "[STRANDNESS] Cached: ${_detected_strand:-unstranded} | HISAT2: ${hisat2_strand_opts:-none} | StringTie: ${stringtie_strand_opt:-none}"
	else
		# Find first sample with trimmed FASTQs
		local _det_srr="" _det_t1="" _det_t2=""
		for _det_srr in "${rnaseq_list[@]}"; do
			find_trimmed_fastq "$_det_srr"
			if [[ -n "$trimmed1" ]]; then
				_det_t1="$trimmed1"; _det_t2="${trimmed2:-}"
				break
			fi
		done

		if [[ -n "$_det_t1" ]]; then
			# Align first sample for strandness detection.  If a stranded
			# library is detected, the BAM is removed so the main loop
			# re-aligns with proper --rna-strandness; otherwise it is kept.
			local _det_dir="$HISAT2_DE_NOVO_ROOT/$_det_srr"
			local _det_bam="$_det_dir/${_det_srr}_${fasta_tag}_trimmed_mapped_sorted.bam"

			if [[ ! -f "$_det_bam" || "${OVERWRITE_MODE:-skip}" == "overwrite" ]]; then
				mkdir -p "$_det_dir"
				log_step "Aligning $_det_srr for strandness auto-detection"
				# Pipe directly into samtools sort — eliminates 10-50GB SAM intermediate
				local _det_summary="$_det_dir/${_det_srr}_${fasta_tag}_detection_summary.txt"
				# Cap sort threads: min(THREADS, 4); at least 1
				local _det_sort_threads=$(( THREADS < 4 ? THREADS : 4 ))
				(( _det_sort_threads < 1 )) && _det_sort_threads=1
				local _det_sort_mem _det_wi_flag=""
				_det_sort_mem=$(_samtools_sort_mem "$_det_sort_threads" 1)
				_samtools_has_write_index && _det_wi_flag="--write-index"

				if [[ -n "$_det_t2" && -f "$_det_t2" ]]; then
					hisat2 -p "$THREADS" --dta -x "$index_prefix" \
						-1 "$_det_t1" -2 "$_det_t2" 2>"$_det_summary" \
						| samtools sort -@ "$_det_sort_threads" -m "$_det_sort_mem" $_det_wi_flag -o "$_det_bam"
				else
					hisat2 -p "$THREADS" --dta -x "$index_prefix" \
						-U "$_det_t1" 2>"$_det_summary" \
						| samtools sort -@ "$_det_sort_threads" -m "$_det_sort_mem" $_det_wi_flag -o "$_det_bam"
				fi
				local _det_ps=("${PIPESTATUS[@]}")
				[[ -s "$_det_summary" ]] && sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g' "$_det_summary"
				if [[ ${_det_ps[0]} -ne 0 || ${_det_ps[1]} -ne 0 ]]; then
					log_warn "[STRANDNESS] Detection alignment failed — running unstranded"
					rm -f "$_det_bam"
				elif [[ -z "$_det_wi_flag" ]]; then
					local _idx_t=$THREADS; (( _idx_t > 4 )) && _idx_t=4
					samtools index -@ "$_idx_t" "$_det_bam" 2>/dev/null
				fi
				# Symlink detection summary to standard alignment summary name so QC can find it
				local _std_summary="$_det_dir/${_det_srr}_${fasta_tag}_alignment_summary.txt"
				[[ -s "$_det_summary" && ! -f "$_std_summary" ]] && ln -sf "${_det_summary##*/}" "$_std_summary"
			fi

			if [[ -f "$_det_bam" ]]; then
				_m2_infer_strand_from_bam "$_det_bam"
				# If stranded library detected, remove detection BAM so the main
				# loop re-aligns this sample with proper --rna-strandness options
				if [[ -n "$hisat2_strand_opts" ]]; then
					log_info "[STRANDNESS] Removing detection BAM — will re-align $_det_srr with $hisat2_strand_opts"
					rm -f "$_det_bam" "${_det_bam}.bai" "${_det_bam}.csi"
				fi
			fi
		else
			log_warn "[STRANDNESS] No trimmed FASTQs found — running unstranded"
		fi

		# Cache for resume runs
		printf '_detected_strand="%s"\nhisat2_strand_opts="%s"\nstringtie_strand_opt="%s"\n' \
			"${_detected_strand:-unstranded}" "$hisat2_strand_opts" "$stringtie_strand_opt" \
			> "$_strand_cache"
	fi

	# ALIGNMENT AND STRINGTIE ASSEMBLY
	# HISAT2+StringTie optimal: 16 threads per job (scales well up to ~16)
	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	if [[ "${_JOBS_MODE:-}" == "auto" || "${_JOBS_MODE:-}" == "AUTO" ]]; then
		parallel_jobs=$(( THREADS / 16 ))
		(( parallel_jobs < 1 )) && parallel_jobs=1
	fi
	# Adaptive: cap parallel jobs at sample count to maximize per-job thread allocation
	local _n_samples=${#rnaseq_list[@]}
	(( parallel_jobs > _n_samples )) && parallel_jobs=$_n_samples
	(( parallel_jobs < 1 )) && parallel_jobs=1
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if $_SHARED_HAS_PARALLEL && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
		log_step "[PARALLEL] HISAT2 De Novo: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export fasta_tag index_prefix threads_per_job hisat2_strand_opts stringtie_strand_opt
		local abs_hisat2_dn_root="$HISAT2_DE_NOVO_ROOT"
		[[ "$abs_hisat2_dn_root" != /* ]] && abs_hisat2_dn_root="$PWD/$abs_hisat2_dn_root"
		local abs_stringtie_dn_root="$STRINGTIE_HISAT2_DE_NOVO_ROOT"
		[[ "$abs_stringtie_dn_root" != /* ]] && abs_stringtie_dn_root="$PWD/$abs_stringtie_dn_root"
		export abs_hisat2_dn_root abs_stringtie_dn_root

		_m2_align_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log HISAT2_DN "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			# Check StringTie output first — if final results exist, skip entirely
			# (BAM may have been cleaned up by keep_bam_global="n")
			local out_dir="$abs_stringtie_dn_root/$SRR"
			local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
			local out_abund="$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo.tsv"

			if [[ -f "$out_gtf" && -f "$out_abund" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_DN "$SRR" INFO "StringTie output exists - skipping"
				return 0
			fi

			local HISAT2_DIR="$abs_hisat2_dn_root/$SRR"
			mkdir -p "$HISAT2_DIR"
			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"

			if [[ -f "$bam" && ( -f "${bam}.bai" || -f "${bam}.csi" ) && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_DN "$SRR" INFO "BAM exists - skipping alignment"
			else
				_parallel_log HISAT2_DN "$SRR" INFO "Aligning with $threads_per_job threads"
				# Pipe hisat2 directly into samtools sort — eliminates 10-50GB SAM intermediate per sample
				# Cap sort threads at 4: I/O-bound beyond that, and sort memory is per-thread
				# Cap sort threads: min(threads_per_job, 4); at least 1
				local sort_threads=$(( threads_per_job < 4 ? threads_per_job : 4 ))
				(( sort_threads < 1 )) && sort_threads=1
				local sort_mem _sort_wi_flag=""
				sort_mem=$(_samtools_sort_mem "$sort_threads" "$parallel_jobs")
				_samtools_has_write_index && _sort_wi_flag="--write-index"
				local _align_log="$HISAT2_DIR/${SRR}_${fasta_tag}_alignment_summary.txt"
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "$threads_per_job" --dta $hisat2_strand_opts -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" 2>"$_align_log" \
						| samtools sort -@ "$sort_threads" -m "$sort_mem" $_sort_wi_flag -o "$bam"
				else
					hisat2 -p "$threads_per_job" --dta $hisat2_strand_opts -x "$index_prefix" \
						-U "$trimmed1" 2>"$_align_log" \
						| samtools sort -@ "$sort_threads" -m "$sort_mem" $_sort_wi_flag -o "$bam"
				fi
				local _ps=("${PIPESTATUS[@]}")
				# Display alignment summary (strip ANSI codes for clean log output)
				[[ -s "$_align_log" ]] && sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g' "$_align_log"
				[[ ${_ps[0]} -ne 0 ]] && { _parallel_log HISAT2_DN "$SRR" ERROR "HISAT2 failed (exit=${_ps[0]})"; rm -f "$bam"; return ${_ps[0]}; }
				[[ ${_ps[1]} -ne 0 ]] && { _parallel_log HISAT2_DN "$SRR" ERROR "samtools sort failed"; rm -f "$bam"; return 1; }
				if [[ -z "$_sort_wi_flag" ]]; then
					# Cap index threads at 4 — samtools index is I/O-bound
					local _idx_t=$threads_per_job; (( _idx_t > 4 )) && _idx_t=4
					samtools index -@ "$_idx_t" "$bam" 2>/dev/null
					[[ $? -ne 0 ]] && { _parallel_log HISAT2_DN "$SRR" ERROR "samtools index failed"; rm -f "$bam"; return 1; }
				fi
			fi

			# StringTie assembly (de novo)
			mkdir -p "$out_dir"

			if [[ -f "$out_gtf" && -f "$out_abund" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_DN "$SRR" INFO "De novo assembly exists - skipping"
			else
				_parallel_log HISAT2_DN "$SRR" INFO "Assembling transcripts (de novo)"
				stringtie -p "$threads_per_job" $stringtie_strand_opt "$bam" -o "$out_gtf" \
					-A "$out_abund" 2>&1 || \
					{ _parallel_log HISAT2_DN "$SRR" ERROR "StringTie failed"; rm -f "$out_gtf" "$out_abund"; return 1; }
			fi

			# Only delete BAM after verifying StringTie produced valid output
			if [[ "$keep_bam_global" != "y" ]]; then
				if [[ -f "$out_gtf" && -s "$out_gtf" ]]; then
					rm -f "$bam" "${bam}.bai" "${bam}.csi"
				else
					_parallel_log HISAT2_DN "$SRR" WARN "Retaining BAM — StringTie GTF missing or empty"
				fi
			fi
			_parallel_log HISAT2_DN "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m2_align_parallel_worker
		# Pre-detect capabilities so parallel workers don't each test independently
		_samtools_has_write_index || true
		_get_available_ram_mb > /dev/null

		parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env index_prefix --env threads_per_job \
			--env hisat2_strand_opts --env stringtie_strand_opt \
			--env abs_hisat2_dn_root --env abs_stringtie_dn_root \
			--env OVERWRITE_MODE --env _SAMTOOLS_HAS_WRITE_INDEX --env _CACHED_AVAIL_MB \
			-j "$parallel_jobs" \
			--halt soon,fail,1 \
			--joblog "$HISAT2_DE_NOVO_ROOT/parallel_hisat2_denovo.log" \
			_m2_align_parallel_worker {} \
			< <(printf '%s\n' "${rnaseq_list[@]}")

		local par_exit=$?
		log_info "[PARALLEL] HISAT2 De Novo complete (exit=$par_exit)"
		if [[ $par_exit -ne 0 ]]; then
			log_warn "[PARALLEL] Some jobs failed - check $HISAT2_DE_NOVO_ROOT/parallel_hisat2_denovo.log"
			return $par_exit
		fi
	else
		# Sequential fallback
		# Pre-compute sort params once (invariant across samples in sequential mode)
		local _seq_sort_threads=$(( THREADS < 4 ? THREADS : 4 ))
		(( _seq_sort_threads < 1 )) && _seq_sort_threads=1
		local _seq_sort_mem _seq_sort_wi_flag=""
		_seq_sort_mem=$(_samtools_sort_mem "$_seq_sort_threads" 1)
		_samtools_has_write_index && _seq_sort_wi_flag="--write-index"
		local _seq_failures=0
		for SRR in "${rnaseq_list[@]}"; do
			# Check StringTie output first — if final results exist, skip entirely
			# (BAM may have been cleaned up by keep_bam_global="n")
			local out_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$SRR"
			local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
			local out_abund="$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo.tsv"

			if [[ -f "$out_gtf" && -f "$out_abund" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[HISAT2 DN] StringTie output exists for $SRR - skipping"
				continue
			fi

			local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/$SRR"
			mkdir -p "$HISAT2_DIR"

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ for $SRR not found - skipping"; continue; }

			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"

			if [[ -f "$bam" && ( -f "${bam}.bai" || -f "${bam}.csi" ) && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[HISAT2 ALIGN] BAM exists for $SRR - skipping alignment"
			else
				log_step "Aligning $SRR using HISAT2 De Novo"

				# Pipe hisat2 directly into samtools sort — eliminates 10-50GB SAM intermediate per sample
				# sort_threads, sort_mem, _sort_wi_flag pre-computed before loop
				local _align_log="$HISAT2_DIR/${SRR}_${fasta_tag}_alignment_summary.txt"
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "${THREADS}" --dta $hisat2_strand_opts -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" 2>"$_align_log" \
						| samtools sort -@ "$_seq_sort_threads" -m "$_seq_sort_mem" $_seq_sort_wi_flag -o "$bam"
				else
					hisat2 -p "${THREADS}" --dta $hisat2_strand_opts -x "$index_prefix" \
						-U "$trimmed1" 2>"$_align_log" \
						| samtools sort -@ "$_seq_sort_threads" -m "$_seq_sort_mem" $_seq_sort_wi_flag -o "$bam"
				fi
				local _ps=("${PIPESTATUS[@]}")
				# Display alignment summary
				[[ -s "$_align_log" ]] && sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g' "$_align_log"
				[[ ${_ps[0]} -ne 0 ]] && { log_error "[HISAT2] Alignment failed for $SRR (exit=${_ps[0]})"; rm -f "$bam"; ((_seq_failures++)) || true; continue; }
				[[ ${_ps[1]} -ne 0 ]] && { log_error "[SAMTOOLS] sort failed for $SRR"; rm -f "$bam"; ((_seq_failures++)) || true; continue; }

				# Only run separate index if --write-index was not used
				if [[ -z "$_seq_sort_wi_flag" ]]; then
					local _idx_t=$THREADS; (( _idx_t > 4 )) && _idx_t=4
					run_with_space_time_log samtools index -@ "$_idx_t" "$bam" \
						|| { log_error "[SAMTOOLS] index failed for $SRR"; rm -f "$bam"; ((_seq_failures++)) || true; continue; }
				fi
			fi

			# StringTie assembly (de novo - no reference GTF)
			# out_dir, out_gtf, out_abund already set above
			mkdir -p "$out_dir"

			if [[ -f "$out_gtf" && -f "$out_abund" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[STRINGTIE] De novo assembly exists for $SRR - skipping"
			else
				log_step "Assembling transcripts for $SRR (de novo)"
				run_with_space_time_log --input "$bam" --output "$out_dir" \
					stringtie -p "$THREADS" $stringtie_strand_opt "$bam" -o "$out_gtf" \
						-A "$out_abund" \
					|| { log_error "[STRINGTIE] Assembly failed for $SRR"; rm -f "$out_gtf" "$out_abund"; ((_seq_failures++)) || true; continue; }
			fi

			# Only delete BAM after verifying StringTie produced valid output
			if [[ "$keep_bam_global" != "y" ]]; then
				if [[ -f "$out_gtf" && -s "$out_gtf" ]]; then
					rm -f "$bam" "${bam}.bai" "${bam}.csi"
				else
					log_warn "[BAM] Retaining $SRR BAM — StringTie GTF missing or empty"
				fi
			fi

			log_info "[STRINGTIE] Done processing $SRR (de novo)"
		done

		if [[ $_seq_failures -gt 0 ]]; then
			log_warn "[SEQUENTIAL] $_seq_failures sample(s) failed during HISAT2 De Novo processing"
			return 1
		fi
	fi

	# POST-ALIGNMENT QC: batch-extract alignment rates from all summary files in single awk pass
	# (replaces N separate awk spawns with 1 process for N samples)
	local _qc_warn=0
	local _summary_files=()
	for SRR in "${rnaseq_list[@]}"; do
		local _sumf="$HISAT2_DE_NOVO_ROOT/$SRR/${SRR}_${fasta_tag}_alignment_summary.txt"
		[[ -f "$_sumf" ]] && _summary_files+=("$_sumf")
	done
	if [[ ${#_summary_files[@]} -gt 0 ]]; then
		# Single awk pass extracts sample name and overall rate from all summary files
		# Uses fasta_tag variable for exact suffix stripping (supports fasta_tags with underscores)
		local _qc_output
		_qc_output=$(awk -v ft="$fasta_tag" '
			FNR==1 { sample = FILENAME; sub(/.*\//, "", sample); sub("_" ft "_alignment_summary\\.txt$", "", sample) }
			/overall alignment rate/ { gsub(/%.*/, ""); print sample, $NF }
		' "${_summary_files[@]}" 2>/dev/null) || true
		while read -r _qc_srr _qc_rate; do
			[[ -z "$_qc_rate" ]] && continue
			local _ov_int=${_qc_rate%.*}
			if (( _ov_int < 50 )); then
				log_warn "[HISAT2 DN QC] $_qc_srr: Overall alignment ${_qc_rate}% — VERY LOW"
				((_qc_warn++)) || true
			elif (( _ov_int < 70 )); then
				log_warn "[HISAT2 DN QC] $_qc_srr: Overall alignment ${_qc_rate}% — below 70% threshold"
				((_qc_warn++)) || true
			fi
		done <<< "$_qc_output"
	fi
	if [[ $_qc_warn -gt 0 ]]; then
		log_warn "[HISAT2 DN QC] $_qc_warn warning(s) — review sample quality"
	else
		log_info "[HISAT2 DN QC] All samples passed alignment rate checks"
	fi

	log_step "HISAT2 de novo pipeline completed for $fasta_tag"
}

# ==============================================================================
# STRANDNESS INFERENCE HELPER
# ==============================================================================
# Infer library strandness from a BAM aligned to a transcriptome FASTA.
# All reference sequences in a transcriptome FASTA are on the + strand, so:
#   RF (dUTP/TruSeq):  read1 maps predominantly reverse  (flag 0x10)
#   FR (ligation):     read1 maps predominantly forward   (no 0x10)
#   Unstranded:        ~50/50 split
# Sets hisat2_strand_opts and stringtie_strand_opt in the caller's scope.
# Usage: _m2_infer_strand_from_bam <bam>
_m2_infer_strand_from_bam() {
	local bam="$1"
	local fwd_count rev_count total fwd_frac is_paired

	# Single samtools view + awk pass replaces 3 separate samtools calls.
	# Exclude unmapped (0x4), secondary (0x100), supplementary (0x800) = 0x904
	# Paired-end: check read1 (0x40) strand via 0x10 flag
	# Single-end: check overall strand via 0x10 flag
	# Uses POSIX-compatible int(flag/N)%2 instead of gawk-specific and()
	# Single samtools|awk pipeline computes counts AND strand decision
	# (1 pipe instead of 1 pipe + 3 extra awk invocations for float math)
	# O(min(N, 200000)) — awk exits early after 200K reads, eliminating the head subprocess
	eval "$(samtools view -F 0x904 "$bam" 2>/dev/null | awk '
		NR > 200000 { exit }
		BEGIN { paired=0; fwd=0; rev=0 }
		{
			flag = $2
			if (int(flag/2) % 2) {
				paired++
				if (int(flag/64) % 2) {   # read1 (0x40)
					if (int(flag/16) % 2) rev++
					else fwd++
				}
			} else {
				if (int(flag/16) % 2) rev++
				else fwd++
			}
		}
		END {
			total = fwd + rev
			frac = (total > 0) ? fwd / total : 0.5
			is_pe = (paired > 0) ? 1 : 0
			# strand: FR if frac>0.7, RF if frac<0.3, else unstranded
			strand = (frac > 0.7) ? "FR" : (frac < 0.3) ? "RF" : "unstranded"
			printf "is_paired=%d fwd_count=%d rev_count=%d total=%d fwd_frac=%.4f _detected_strand=%s", \
				paired, fwd, rev, total, frac, strand
		}
	')"

	if [[ "${total:-0}" -eq 0 ]]; then
		log_warn "[STRANDNESS] No mapped reads in BAM — running unstranded"
		_detected_strand="unstranded"
		return 0
	fi

	local _is_pe=$(( ${is_paired:-0} > 0 ? 1 : 0 ))
	log_info "[STRANDNESS] Read1 forward fraction: $fwd_frac ($total primary alignments, PE=$_is_pe)"

	if [[ "$_detected_strand" == "FR" ]]; then
		# HISAT2 requires FR for paired-end, F for single-end
		if [[ $_is_pe -eq 1 ]]; then
			hisat2_strand_opts="--rna-strandness FR"
		else
			hisat2_strand_opts="--rna-strandness F"
		fi
		stringtie_strand_opt="--fr"
		log_info "[STRANDNESS] Auto-detected: FR (ligation / forward-stranded)"
	elif [[ "$_detected_strand" == "RF" ]]; then
		# HISAT2 requires RF for paired-end, R for single-end
		if [[ $_is_pe -eq 1 ]]; then
			hisat2_strand_opts="--rna-strandness RF"
		else
			hisat2_strand_opts="--rna-strandness R"
		fi
		stringtie_strand_opt="--rf"
		log_info "[STRANDNESS] Auto-detected: RF (dUTP / TruSeq / reverse-stranded)"
	else
		_detected_strand="unstranded"
		log_info "[STRANDNESS] Auto-detected: Unstranded (fwd_frac=$fwd_frac)"
	fi
}
