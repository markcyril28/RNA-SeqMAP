#!/bin/bash
# ==============================================================================
# METHOD 1: HISAT2 REFERENCE GUIDED PIPELINE
# ==============================================================================
# HISAT2 Reference Guided alignment with StringTie assembly
# Uses reference GTF for splice site information
# Pure quantification mode: single StringTie pass with -e (no merge/re-estimation)
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M1_HISAT2_REF_SOURCED:-}" == "true" ]] && return 0
export M1_HISAT2_REF_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# POST-ALIGNMENT QC
# ==============================================================================
# Parse HISAT2 alignment summary files to flag low mapping rates, high
# multi-mapping (potential rRNA contamination), and cohort-level outliers.
# Mirrors _star_check_alignment_rates() in m3_star_alignment.sh.
# Usage: _hisat2_check_alignment_rates <align_dir> <fasta_tag> <srr1> [srr2 ...]

_hisat2_check_alignment_rates() {
	local align_dir="$1" ft="$2"; shift 2
	local srr_list=("$@")
	local warn_count=0

	# Collect per-sample stats for cohort-level outlier detection
	local -a sample_names=() overall_rates=() concordant_once=() concordant_multi=()

	for SRR in "${srr_list[@]}"; do
		local sumf="${align_dir}/${SRR}/${SRR}_${ft}_ref_guided_alignment_summary.txt"
		[[ ! -f "$sumf" ]] && continue

		# HISAT2 summary format (paired-end):
		#   <N> reads; of these:
		#     <N> (...) were paired; of these:
		#       <N> (...) aligned concordantly 0 times
		#       <N> (...) aligned concordantly exactly 1 time
		#       <N> (...) aligned concordantly >1 times
		#   <X>% overall alignment rate
		local overall conc1 concm
		overall=$(grep -oP '[0-9.]+(?=% overall alignment rate)' "$sumf" 2>/dev/null | head -1)
		conc1=$(grep 'aligned concordantly exactly 1 time' "$sumf" 2>/dev/null \
			| grep -oP '[0-9.]+(?=%)' | head -1)
		concm=$(grep 'aligned concordantly >1 time' "$sumf" 2>/dev/null \
			| grep -oP '[0-9.]+(?=%)' | head -1)

		# Single-end fallback: "aligned exactly 1 time" / "aligned >1 times"
		if [[ -z "$conc1" ]]; then
			conc1=$(grep 'aligned exactly 1 time' "$sumf" 2>/dev/null \
				| grep -oP '[0-9.]+(?=%)' | head -1)
			concm=$(grep 'aligned >1 time' "$sumf" 2>/dev/null \
				| grep -oP '[0-9.]+(?=%)' | head -1)
		fi

		[[ -z "$overall" ]] && continue

		sample_names+=("$SRR")
		overall_rates+=("$overall")
		concordant_once+=("${conc1:-0}")
		concordant_multi+=("${concm:-0}")

		# Per-sample checks
		if awk "BEGIN{exit !($overall < 50)}" 2>/dev/null; then
			log_warn "[HISAT2 QC] $SRR: Overall alignment ${overall}% — VERY LOW (check sample quality, adapter contamination, or genome mismatch)"
			((warn_count++)) || true
		elif awk "BEGIN{exit !($overall < 70)}" 2>/dev/null; then
			log_warn "[HISAT2 QC] $SRR: Overall alignment ${overall}% — below 70% threshold"
			((warn_count++)) || true
		fi

		# High multi-mapping (>20%) may indicate rRNA contamination
		if awk "BEGIN{exit !(${concm:-0} > 20)}" 2>/dev/null; then
			log_warn "[HISAT2 QC] $SRR: Multi-mapped ${concm}% — may indicate rRNA contamination or repetitive element enrichment"
			((warn_count++)) || true
		fi
	done

	# Cohort-level outlier detection: flag samples >2 SD below mean overall rate
	local n=${#overall_rates[@]}
	if [[ $n -ge 3 ]]; then
		# Single AWK pass for sum, mean, sd, threshold (replaces 2N+3 AWK spawns)
		local _stats
		_stats=$(printf '%s\n' "${overall_rates[@]}" | awk '{s+=$1; ss+=$1*$1} END{
			m=s/NR; v=ss/NR - m*m; sd=(v>0)?sqrt(v):0
			printf "%.2f %.2f %.2f", m, sd, m-2*sd
		}')
		local mean sd threshold
		read -r mean sd threshold <<< "$_stats"

		log_info "[HISAT2 QC] Cohort alignment stats: mean=${mean}%, SD=${sd}%, outlier threshold=${threshold}%"

		for ((i=0; i<n; i++)); do
			if awk "BEGIN{exit !(${overall_rates[$i]} < $threshold)}" 2>/dev/null; then
				log_warn "[HISAT2 QC] OUTLIER: ${sample_names[$i]} (${overall_rates[$i]}%) is >2 SD below cohort mean (${mean}%)"
				((warn_count++)) || true
			fi
		done
	fi

	# Summary table
	if [[ ${#sample_names[@]} -gt 0 ]]; then
		log_info "[HISAT2 QC] ┌───────────────────┬──────────┬──────────┬──────────┐"
		log_info "[HISAT2 QC] │ Sample            │ Overall  │ Conc.1x  │ Conc.>1x │"
		log_info "[HISAT2 QC] ├───────────────────┼──────────┼──────────┼──────────┤"
		for ((i=0; i<${#sample_names[@]}; i++)); do
			printf -v _row "[HISAT2 QC] │ %-17s │  %5s%%  │  %5s%%  │  %5s%%  │" \
				"${sample_names[$i]}" "${overall_rates[$i]}" "${concordant_once[$i]}" "${concordant_multi[$i]}"
			log_info "$_row"
		done
		log_info "[HISAT2 QC] └───────────────────┴──────────┴──────────┴──────────┘"
	fi

	if [[ $warn_count -gt 0 ]]; then
		log_warn "[HISAT2 QC] $warn_count warning(s) detected — review samples before proceeding"
	else
		log_info "[HISAT2 QC] All samples passed alignment rate checks"
	fi
}

# ==============================================================================
# FASTA/GTF CHROMOSOME VALIDATION
# ==============================================================================
# Verify that the GTF chromosome names match the FASTA sequence names.
# A mismatch means HISAT2 splice site guidance is silently ineffective.

_m1_validate_fasta_gtf_chromosomes() {
	local fasta="$1" gtf="$2"

	# Extract first 20 unique chromosome names from each file
	local fasta_chrs gtf_chrs
	fasta_chrs=$(grep '^>' "$fasta" | head -20 | sed 's/^>//; s/[[:space:]].*//' | sort)
	gtf_chrs=$(awk '$1 !~ /^#/ {print $1}' "$gtf" | sort -u | head -20)

	if [[ -z "$fasta_chrs" || -z "$gtf_chrs" ]]; then
		log_warn "[VALIDATE] Could not extract chromosome names from FASTA or GTF"
		return 0
	fi

	# Count overlapping chromosome names
	local overlap
	overlap=$(comm -12 <(echo "$fasta_chrs") <(echo "$gtf_chrs") | wc -l)
	local fasta_count gtf_count
	fasta_count=$(echo "$fasta_chrs" | wc -l)
	gtf_count=$(echo "$gtf_chrs" | wc -l)

	if [[ "$overlap" -eq 0 ]]; then
		log_warn "[VALIDATE] FASTA/GTF MISMATCH: No shared chromosome names between FASTA and GTF!"
		log_warn "[VALIDATE]   FASTA chromosomes (first 5): $(echo "$fasta_chrs" | head -5 | tr '\n' ' ')"
		log_warn "[VALIDATE]   GTF chromosomes (first 5):   $(echo "$gtf_chrs" | head -5 | tr '\n' ' ')"
		log_warn "[VALIDATE]   Splice site guidance will be INEFFECTIVE — alignment quality degraded"
	elif [[ "$overlap" -lt "$gtf_count" ]]; then
		log_info "[VALIDATE] FASTA/GTF partial overlap: $overlap of $gtf_count GTF chromosomes found in FASTA"
	else
		log_info "[VALIDATE] FASTA/GTF chromosome names match ($overlap shared)"
	fi
}

# ==============================================================================
# BAM METRICS: SOFT-CLIPPING & INSERT SIZE
# ==============================================================================
# Collect lightweight BAM metrics using samtools stats. Must run BEFORE BAM deletion.
# Saves a summary TSV per sample for post-hoc review.
# Usage: _m1_collect_bam_metrics <bam> <output_dir> [method_tag] [srr_tag]

_m1_collect_bam_metrics() {
	local bam="$1" out_dir="$2"
	local method="${3:-HISAT2_RG}" srr="${4:-SAMPLE}"
	local metrics_file="$out_dir/${srr}_bam_metrics.txt"

	if ! command -v samtools >/dev/null 2>&1 || [[ ! -f "$bam" ]]; then
		return 0
	fi

	# samtools stats is fast — runs in seconds even on large BAMs
	local stats
	stats=$(samtools stats "$bam" 2>/dev/null | grep '^SN\t') || return 0

	local total_bases bases_clipped insert_mean insert_sd
	total_bases=$(echo "$stats" | awk -F'\t' '/^SN\tbases mapped \(cigar\)/{print $3}')
	bases_clipped=$(echo "$stats" | awk -F'\t' '/^SN\tbases trimmed/{print $3}')
	insert_mean=$(echo "$stats" | awk -F'\t' '/^SN\tinsert size average/{print $3}')
	insert_sd=$(echo "$stats" | awk -F'\t' '/^SN\tinsert size standard deviation/{print $3}')

	# Save metrics to file for post-hoc review
	{
		echo "sample=$srr"
		echo "total_mapped_bases=$total_bases"
		echo "bases_soft_clipped=$bases_clipped"
		echo "insert_size_mean=$insert_mean"
		echo "insert_size_sd=$insert_sd"
	} > "$metrics_file"

	# Flag excessive soft-clipping (>10% of mapped bases)
	if [[ -n "$total_bases" && -n "$bases_clipped" && "$total_bases" -gt 0 ]]; then
		local clip_pct
		clip_pct=$(awk "BEGIN{printf \"%.1f\", 100*$bases_clipped/$total_bases}")
		if awk "BEGIN{exit !($clip_pct > 10)}" 2>/dev/null; then
			_parallel_log "$method" "$srr" WARN "Soft-clipped ${clip_pct}% of mapped bases — may indicate adapter contamination or index mismatch"
		fi
	fi

	# Flag abnormal insert size (outside 100-800 bp for typical RNA-seq)
	if [[ -n "$insert_mean" ]]; then
		local mean_int=${insert_mean%.*}
		if [[ "$mean_int" -gt 0 ]]; then
			if [[ "$mean_int" -lt 100 || "$mean_int" -gt 800 ]]; then
				_parallel_log "$method" "$srr" WARN "Unusual insert size: mean=${insert_mean}, SD=${insert_sd:-N/A} (expected 100-800 bp for RNA-seq)"
			fi
		fi
	fi
}

# ==============================================================================
# HISAT2 REFERENCE GUIDED PIPELINE
# ==============================================================================

hisat2_ref_guided_pipeline() {
	local fasta="" gtf="" strandness="" rnaseq_list=()

	# Parse arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--GTF) gtf="$2"; shift 2;;
			--STRANDNESS) strandness="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done

	# Validate inputs
	[[ -z "$fasta" ]] && { log_error "No FASTA file specified. Use --FASTA <genome_fasta>."; return 1; }
	[[ ! -f "$fasta" ]] && { log_error "FASTA file '$fasta' not found."; return 1; }
	[[ -z "$gtf" ]] && { log_error "No GTF file specified. Use --GTF <annotation_gtf>."; return 1; }
	[[ ! -f "$gtf" ]] && { log_error "GTF file '$gtf' not found."; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")

	# Derive strandness flags from --STRANDNESS (RF | FR | "" for unstranded)
	local hisat2_strand_opts="" stringtie_strand_opt=""
	if [[ -n "$strandness" ]]; then
		if [[ "$strandness" != "RF" && "$strandness" != "FR" ]]; then
			log_error "[STRANDNESS] Invalid value '$strandness'. Must be RF (dUTP/TruSeq) or FR (ligation protocol)."; return 1
		fi
		hisat2_strand_opts="--rna-strandness $strandness"
		[[ "$strandness" == "RF" ]] && stringtie_strand_opt="--rf"
		[[ "$strandness" == "FR" ]] && stringtie_strand_opt="--fr"
		log_info "[STRANDNESS] HISAT2: $hisat2_strand_opts | StringTie: $stringtie_strand_opt"
	else
		log_warn "[STRANDNESS] No strandness specified — running unstranded. Pass --STRANDNESS RF or FR for stranded libraries."
	fi

	local fasta_base fasta_tag index_prefix
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	set_fasta_output_dirs "$fasta_tag"
	index_prefix="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_ref_guided"

	# Validate FASTA/GTF chromosome name consistency
	_m1_validate_fasta_gtf_chromosomes "$fasta" "$gtf"

	# BUILD HISAT2 REFERENCE-GUIDED INDEX
	mkdir -p "$HISAT2_REF_GUIDED_INDEX_DIR"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1 && [[ "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[INDEX] Ref-Guided index exists - skipping build"
	else
		log_step "Building HISAT2 Ref-Guided index: $fasta_base"

		local splice_sites="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_splice_sites.txt"
		local exons="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_exons.txt"
		local build_opts=""

		# Extract splice sites and exons (may be empty for single-exon transcriptomes)
		if ! hisat2_extract_splice_sites.py "$gtf" > "$splice_sites" 2>&1; then
			log_warn "[INDEX] hisat2_extract_splice_sites.py failed for $gtf — continuing without splice sites"
			> "$splice_sites"
		fi
		if ! hisat2_extract_exons.py "$gtf" > "$exons" 2>&1; then
			log_warn "[INDEX] hisat2_extract_exons.py failed for $gtf — continuing without exons"
			> "$exons"
		fi

		# Validate annotation coordinates against FASTA sequence lengths.
		# Transcriptome FASTAs contain short spliced sequences, but GTFs may carry
		# genomic coordinates that exceed those lengths, causing hisat2-build to fail
		# with 'Nongraph exception'. Skip --ss/--exon when this mismatch is detected.
		if [[ -s "$splice_sites" ]]; then
			local _ss_chr _ss_pos _ss_seq_len
			read -r _ss_chr _ss_pos _ < "$splice_sites"
			if [[ -n "$_ss_chr" ]]; then
				_ss_seq_len=$(awk -v t="$_ss_chr" \
					'/^>/{if(f){print l; f=0; exit} n=$1; sub(/^>/,"",n); if(n==t){f=1; l=0}; next}
					 f{l+=length} END{if(f) print l}' "$fasta")
				if [[ -n "$_ss_seq_len" && "$_ss_pos" -gt "$_ss_seq_len" ]]; then
					log_warn "[INDEX] Splice site coords exceed FASTA sequence lengths (pos $_ss_pos > seq len $_ss_seq_len for $_ss_chr)"
					log_warn "[INDEX] Transcriptome FASTA detected - skipping --ss/--exon (not applicable)"
					> "$splice_sites"
					> "$exons"
				fi
			fi
		fi

		# Build options based on available annotation data
		if [[ -s "$splice_sites" ]]; then
			log_info "[INDEX] Extracted $(wc -l < "$splice_sites") splice sites"
			build_opts="$build_opts --ss $splice_sites"
		else
			log_warn "[INDEX] No splice sites found - GTF may contain only single-exon transcripts"
		fi

		if [[ -s "$exons" ]]; then
			log_info "[INDEX] Extracted $(wc -l < "$exons") exons"
			build_opts="$build_opts --exon $exons"
		else
			log_warn "[INDEX] No exons extracted from GTF"
		fi

		log_file_size "$fasta" "Input FASTA for HISAT2 index"
		run_with_space_time_log --input "$fasta" --output "$HISAT2_REF_GUIDED_INDEX_DIR" \
			hisat2-build -p "${THREADS}" $build_opts "$fasta" "$index_prefix"
		log_file_size "$HISAT2_REF_GUIDED_INDEX_DIR" "HISAT2 index output"
	fi

	# ALIGNMENT AND STRINGTIE ASSEMBLY (pure quantification: single pass with -e)
	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
		log_step "[PARALLEL] HISAT2 Ref-Guided Align+StringTie: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export fasta_tag index_prefix threads_per_job hisat2_strand_opts stringtie_strand_opt OVERWRITE_MODE
		local abs_hisat2_rg_root="$HISAT2_REF_GUIDED_ROOT"
		[[ "$abs_hisat2_rg_root" != /* ]] && abs_hisat2_rg_root="$(pwd)/$abs_hisat2_rg_root"
		local abs_stringtie_rg_root="$STRINGTIE_HISAT2_REF_GUIDED_ROOT"
		[[ "$abs_stringtie_rg_root" != /* ]] && abs_stringtie_rg_root="$(pwd)/$abs_stringtie_rg_root"
		local abs_gtf="$gtf"
		[[ "$abs_gtf" != /* ]] && abs_gtf="$(pwd)/$abs_gtf"
		export abs_hisat2_rg_root abs_stringtie_rg_root abs_gtf

		_m1_align_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log HISAT2_RG "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local HISAT2_DIR="$abs_hisat2_rg_root/$SRR"
			mkdir -p "$HISAT2_DIR"
			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
			local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped.sam"
			local out_gtf="$abs_stringtie_rg_root/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"

			if [[ -f "$bam" && -f "${bam}.bai" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_RG "$SRR" INFO "BAM exists - skipping alignment"
			elif [[ ! -f "$bam" && -f "$out_gtf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_RG "$SRR" INFO "GTF exists, BAM already cleaned for $SRR - skipping alignment"
			else
				_parallel_log HISAT2_RG "$SRR" INFO "Aligning with $threads_per_job threads"
				local align_exit=0
				local summary_file="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_alignment_summary.txt"
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "$threads_per_job" --dta $hisat2_strand_opts -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" -S "$sam" 2>"$summary_file"
					align_exit=$?
				else
					hisat2 -p "$threads_per_job" --dta $hisat2_strand_opts -x "$index_prefix" \
						-U "$trimmed1" -S "$sam" 2>"$summary_file"
					align_exit=$?
				fi
				# Display alignment summary (strip ANSI codes for clean log output)
				[[ -s "$summary_file" ]] && sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g' "$summary_file"
				[[ $align_exit -ne 0 ]] && { _parallel_log HISAT2_RG "$SRR" ERROR "HISAT2 failed (exit=$align_exit)"; rm -f "$sam"; return $align_exit; }

				samtools sort -@ "$threads_per_job" -o "$bam" "$sam" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				[[ ${PIPESTATUS[0]} -ne 0 ]] && { rm -f "$sam"; _parallel_log HISAT2_RG "$SRR" ERROR "samtools sort failed"; return 1; }
				samtools index -@ "$threads_per_job" "$bam" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g' || true
				rm -f "$sam"

				# Infer strandness once (lock-file ensures only first worker runs it)
				[[ -z "$hisat2_strand_opts" ]] && \
					_m1_infer_strandness "$bam" "$abs_gtf" "$(dirname "$index_prefix")" "HISAT2_RG" "$SRR"
			fi

			# StringTie quantification (ref-guided, single pass)
			local out_dir="$abs_stringtie_rg_root/$SRR"
			local ballgown_dir="$out_dir/ballgown"
			mkdir -p "$out_dir" "$ballgown_dir"

			if [[ -f "$out_gtf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_RG "$SRR" INFO "Assembly exists - skipping"
			else
				_parallel_log HISAT2_RG "$SRR" INFO "Quantifying transcripts (ref-guided)"
				stringtie -e $stringtie_strand_opt -p "$threads_per_job" "$bam" -G "$abs_gtf" -o "$out_gtf" \
					-A "$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv" \
					-B -C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf" 2>&1 || \
					{ _parallel_log HISAT2_RG "$SRR" ERROR "StringTie failed"; return 1; }
			fi

			# Collect BAM metrics before potential deletion
			[[ -f "$bam" ]] && _m1_collect_bam_metrics "$bam" "$HISAT2_DIR" "HISAT2_RG" "$SRR"

			if [[ "$keep_bam_global" != "y" && -f "$bam" ]]; then
				_parallel_log HISAT2_RG "$SRR" WARN "Deleting BAM to save disk (set keep_bam_global=y to retain): $(basename "$bam")"
				rm -f "$bam" "${bam}.bai"
			fi

			_parallel_log HISAT2_RG "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m1_align_parallel_worker _m1_infer_strandness _hisat2_check_alignment_rates _m1_collect_bam_metrics

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env index_prefix --env threads_per_job \
			--env abs_hisat2_rg_root --env abs_stringtie_rg_root --env abs_gtf \
			--env hisat2_strand_opts --env stringtie_strand_opt \
			--env OVERWRITE_MODE \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$HISAT2_REF_GUIDED_ROOT/parallel_hisat2_refguided_align.log" \
			_m1_align_parallel_worker {}

		local par_exit=$?
		log_info "[PARALLEL] HISAT2 Ref-Guided align+assembly complete (exit=$par_exit)"
		[[ $par_exit -ne 0 ]] && log_warn "[PARALLEL] Some jobs failed - check $HISAT2_REF_GUIDED_ROOT/parallel_hisat2_refguided_align.log"
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local HISAT2_DIR="$HISAT2_REF_GUIDED_ROOT/$SRR"
			mkdir -p "$HISAT2_DIR"

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ not found for $SRR - skipping"; continue; }

			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
			local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped.sam"
			local out_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"

			if [[ -f "$bam" && -f "${bam}.bai" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[ALIGN] BAM exists for $SRR/$fasta_tag - skipping"
			elif [[ ! -f "$bam" && -f "$out_gtf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[ALIGN] GTF exists, BAM already cleaned for $SRR - skipping alignment"
			else
				log_step "Aligning: $SRR -> $fasta_tag (HISAT2 Ref-Guided)"
				local summary_file="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_alignment_summary.txt"

				# Run hisat2 directly (not via run_with_space_time_log) to cleanly
				# capture alignment summary from stderr into a separate file for QC.
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "${THREADS}" --dta $hisat2_strand_opts -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" -S "$sam" 2>"$summary_file"
				else
					hisat2 -p "${THREADS}" --dta $hisat2_strand_opts -x "$index_prefix" \
						-U "$trimmed1" -S "$sam" 2>"$summary_file"
				fi
				# Display alignment summary
				[[ -s "$summary_file" ]] && cat "$summary_file"

				# Verify alignment produced a SAM file before proceeding
				if [[ ! -s "$sam" ]]; then
					log_error "[ALIGN] HISAT2 failed to produce SAM for $SRR — skipping sample"
					rm -f "$sam"
					continue
				fi

				log_info "[SAMTOOLS] Converting to sorted BAM..."
				run_with_space_time_log --input "$sam" --output "$bam" samtools sort -@ "${THREADS}" -o "$bam" "$sam"
				run_with_space_time_log samtools index -@ "${THREADS}" "$bam"
				rm -f "$sam"

				# Verify BAM was created before downstream steps
				if [[ ! -f "$bam" ]]; then
					log_error "[ALIGN] samtools sort/index failed for $SRR — skipping sample"
					continue
				fi

				# Infer strandness once on the first sample
				if [[ -z "$strandness" && -z "${_m1_strand_inferred:-}" ]]; then
					_m1_strand_inferred=1
					_m1_infer_strandness "$bam" "$gtf" "$HISAT2_REF_GUIDED_INDEX_DIR" "HISAT2_RG" "$SRR"
				fi
			fi

			# StringTie quantification (ref-guided, single pass)
			local out_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR"
			local ballgown_dir="$out_dir/ballgown"
			mkdir -p "$out_dir" "$ballgown_dir"

			if [[ -f "$out_gtf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[STRINGTIE] Quantification exists for $SRR/$fasta_tag - skipping"
			else
				log_step "Quantifying transcripts: $SRR -> $fasta_tag"
				run_with_space_time_log --input "$bam" --output "$out_dir" \
					stringtie -e $stringtie_strand_opt -p "$THREADS" "$bam" -G "$gtf" -o "$out_gtf" \
						-A "$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv" \
						-B -C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf"
			fi

			# Collect BAM metrics before potential deletion
			[[ -f "$bam" ]] && _m1_collect_bam_metrics "$bam" "$HISAT2_DIR" "HISAT2_RG" "$SRR"

			if [[ "$keep_bam_global" != "y" && -f "$bam" ]]; then
				log_warn "[BAM] Deleting $SRR BAM to save disk (set keep_bam_global=y to retain)"
				rm -f "$bam" "${bam}.bai"
			fi
		done
	fi

	# POST-ALIGNMENT QC: check alignment rates and flag outliers
	_hisat2_check_alignment_rates "$HISAT2_REF_GUIDED_ROOT" "$fasta_tag" "${rnaseq_list[@]}"

	# PREPARE COUNT MATRICES
	log_step "Preparing count matrices for DESeq2"
	local deseq2_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/deseq2_input"
	local prepde_sample_list="$deseq2_dir/sample_list.txt"
	local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
	local transcript_count_matrix="$deseq2_dir/transcript_count_matrix.csv"
	mkdir -p "$deseq2_dir"

	local prepde_list_content="" samples_found=0
	for SRR in "${rnaseq_list[@]}"; do
		local assembled_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		if [[ -f "$assembled_gtf" ]]; then
			prepde_list_content+="$SRR $assembled_gtf"$'\n'
			(( samples_found++ )) || true
		fi
	done

	[[ $samples_found -lt 2 ]] && { log_error "Insufficient samples: $samples_found (need ≥2)"; return 1; }
	printf '%s' "$prepde_list_content" > "$prepde_sample_list"

	if [[ ! -f "$gene_count_matrix" || "${OVERWRITE_MODE:-skip}" == "overwrite" ]]; then
		# Auto-detect read length from first available trimmed FASTQ
		local read_length=""
		for SRR in "${rnaseq_list[@]}"; do
			find_trimmed_fastq "$SRR"
			if [[ -n "$trimmed1" && -f "$trimmed1" ]]; then
				read_length=$(detect_read_length "$trimmed1" 0)
				if [[ -n "$read_length" && "$read_length" -gt 0 ]]; then
					break
				fi
				read_length=""
			fi
		done
		if [[ -z "$read_length" || "$read_length" -eq 0 ]]; then
			log_error "[PREPDE] Failed to auto-detect read length from any trimmed FASTQ. Cannot run prepDE.py without accurate read length."
			log_error "[PREPDE] Ensure trimmed FASTQs exist in $TRIM_DIR_ROOT for at least one sample."
			return 1
		fi
		log_info "[PREPDE] Auto-detected read length: ${read_length} bp (from first available trimmed FASTQ)"

		if command -v prepDE.py >/dev/null 2>&1; then
			run_with_space_time_log prepDE.py -i "$prepde_sample_list" \
				-g "$gene_count_matrix" -t "$transcript_count_matrix" -l "$read_length"
			if [[ ! -f "$gene_count_matrix" ]]; then
				log_error "[PREPDE] prepDE.py did not produce gene count matrix: $gene_count_matrix"
				return 1
			fi
		else
			log_error "prepDE.py not found"
			return 1
		fi
	fi

	# Create sample metadata
	local sample_metadata="$deseq2_dir/sample_metadata.csv"
	[[ ! -f "$sample_metadata" ]] && create_sample_metadata "$sample_metadata" "${rnaseq_list[@]}"

	# Validate outputs
	[[ -f "$gene_count_matrix" ]] && validate_count_matrix "$gene_count_matrix" "gene" 2

	log_step "HISAT2 reference-guided pipeline completed for $fasta_tag"
	log_info "Gene count matrix: $gene_count_matrix"
	log_info "Sample metadata: $sample_metadata"
}

# ==============================================================================
# STRANDNESS INFERENCE HELPER
# ==============================================================================
# Build BED12 from GTF (cached) and run infer_experiment.py on a BAM once.
# A lock file (infer_strandness.done) ensures only the first caller executes.
# val_1 = FR (1++,1--,2+-,2-+) / val_2 = RF (1+-,1-+,2++,2--)
# Usage: _m1_infer_strandness <bam> <gtf> <index_dir> [method_tag] [srr_tag]
_m1_infer_strandness() {
	local bam="$1" gtf="$2" index_dir="$3"
	local method="${4:-STRANDNESS}" srr="${5:-PIPELINE}"
	local sentinel="$index_dir/infer_strandness.done"
	local bed12="$index_dir/annotation_infer_exp.bed"

	# Atomic lock — only the first concurrent caller proceeds
	( set -o noclobber; : > "$sentinel" ) 2>/dev/null || return 0

	if ! command -v infer_experiment.py >/dev/null 2>&1; then
		_parallel_log "$method" "$srr" WARN "infer_experiment.py not found — skipping strandness check"
		return 0
	fi

	# Build BED12 from GTF (cached in index dir)
	if [[ ! -s "$bed12" ]]; then
		if command -v gtfToGenePred >/dev/null 2>&1 && command -v genePredToBed >/dev/null 2>&1; then
			gtfToGenePred "$gtf" /dev/stdout 2>/dev/null | genePredToBed /dev/stdin "$bed12" 2>/dev/null
		else
			_parallel_log "$method" "$srr" WARN "gtfToGenePred not found — cannot build BED12 for strandness check"
			return 0
		fi
	fi

	[[ ! -s "$bed12" ]] && { _parallel_log "$method" "$srr" WARN "BED12 conversion from GTF failed"; return 0; }

	_parallel_log "$method" "$srr" INFO "Running infer_experiment.py on: $(basename "$bam")"
	local result
	result=$(infer_experiment.py -i "$bam" -r "$bed12" 2>/dev/null)
	printf '%s\n' "$result" >> "$sentinel"

	# Parse: val_1 = FR (forward-stranded), val_2 = RF (reverse-stranded / dUTP)
	# Paired-end: "1++,1--,2+-,2-+" (FR) / "1+-,1-+,2++,2--" (RF)
	# Single-end:       "++,--"      (FR) /       "+-,-+"       (RF)
	local val1 val2
	val1=$(printf '%s' "$result" | grep -i '"1++,1--,2+-,2-+"' | grep -oP '[0-9]+\.[0-9]+' | tail -1)
	val2=$(printf '%s' "$result" | grep -i '"1+-,1-+,2++,2--"' | grep -oP '[0-9]+\.[0-9]+' | tail -1)
	# Fallback to single-end patterns if paired-end yielded nothing
	if [[ -z "$val1" && -z "$val2" ]]; then
		val1=$(printf '%s' "$result" | grep -i '"++,--"' | grep -oP '[0-9]+\.[0-9]+' | tail -1)
		val2=$(printf '%s' "$result" | grep -i '"+-,-+"' | grep -oP '[0-9]+\.[0-9]+' | tail -1)
	fi

	local rec
	if   [[ -n "$val1" ]] && awk "BEGIN{exit !($val1 > 0.6)}"; then
		rec="FR  → add --STRANDNESS FR to your config (ligation / forward-stranded)"
	elif [[ -n "$val2" ]] && awk "BEGIN{exit !($val2 > 0.6)}"; then
		rec="RF  → add --STRANDNESS RF to your config (dUTP / TruSeq / reverse-stranded)"
	else
		rec="Unstranded  → omit --STRANDNESS (already the default)"
	fi

	_parallel_log "$method" "$srr" INFO "[STRANDNESS]   _val_1 (FR  1++,1--,2+-,2-+): ${val1:-N/A}"
	_parallel_log "$method" "$srr" INFO "[STRANDNESS]   _val_2 (RF  1+-,1-+,2++,2--): ${val2:-N/A}"
	_parallel_log "$method" "$srr" WARN "[STRANDNESS]   Recommendation: $rec"
}
