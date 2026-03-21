#!/bin/bash
# ==============================================================================
# METHOD 3: STAR SPLICE-AWARE ALIGNMENT PIPELINE
# ==============================================================================
# STAR splice-aware alignment with Salmon quantification and tximport for DESeq2
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M3_STAR_SOURCED:-}" == "true" ]] && return 0
export M3_STAR_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_M3_SCRIPT_DIR="$SCRIPT_DIR"
source "$SCRIPT_DIR/shared_utils_method.sh"

# Binary availability cached in shared_utils_method.sh: _SHARED_HAS_PARALLEL, _SHARED_HAS_SAMTOOLS

# ==============================================================================
# STAR CONFIGURATION - IMPORTANT PARAMETERS (tweak here)
# ==============================================================================

# Genome loading mode: NoSharedMemory (default/safe), LoadAndKeep (multi-run HPC)
# Already set by global_config_method.sh; re-state default here for clarity.
STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-NoSharedMemory}"

# sjdbOverhang = read_length - 1; set to actual read length for best splice detection
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"

# Delete intermediate STAR files after alignment to save disk space (true/false)
# Set to "false" to keep 2-pass genome, unsorted BAM, etc. for debugging
STAR_DELETE_TRANSIENT="${STAR_DELETE_TRANSIENT:-true}"

# Use shared pigz detection from shared_utils_method.sh; set STAR-specific aliases
_PIGZ_DC="${_SHARED_GZIP_DC:-gzip -dc}"
_STAR_READ_CMD="$_PIGZ_DC"

# _samtools_sort_mem() is now defined in shared_utils_method.sh (shared across M1-M5)

# Calculate total RAM budget (bytes) for STAR's internal BAM sorting.
# Used with --limitBAMsortRAM when --outSAMtype BAM SortedByCoordinate.
# Reserves 30% of available RAM for STAR alignment + OS, gives the rest to sorting.
# Usage: _star_sort_ram [parallel_jobs]
_star_sort_ram() {
	local parallel_jobs="${1:-1}"
	local avail_mb
	avail_mb=$(_get_available_ram_mb)

	# 70% of available for sorting, split across concurrent jobs
	local per_job_mb=$(( avail_mb * 70 / 100 / parallel_jobs ))
	# Floor at 2GB, cap at 32GB per job
	[[ $per_job_mb -lt 2048 ]] && per_job_mb=2048
	[[ $per_job_mb -gt 32768 ]] && per_job_mb=32768

	# STAR wants bytes
	echo $(( per_job_mb * 1048576 ))
}

# REPRODUCIBILITY NOTE: Salmon's EM algorithm convergence is thread-schedule dependent.
# For exact reproducibility, always use the same --threads value across runs.
# Salmon does not expose a --seed flag for its internal EM algorithm.

# GTF annotation file for splice junction detection (resolved at runtime in pipeline)
# NOTE: STAR_GTF_FILE is resolved inside star_alignment_pipeline() from $gtf_file

# Transcriptome FASTA for Salmon quantification step (auto-detected if not set)
# When using a genome FASTA for STAR, point this to the transcript-level FASTA.
# e.g. STAR_TRANSCRIPTOME_FASTA="inputs/fasta/reference_genomes/GPE001970_transcripts.fa"
# Auto-detect: looks for <genome_basename>_transcripts.fa; falls back to --FASTA if absent.

# STAR temp directory: always uses the output directory per sample
# (${star_genome_dir}/_STARtmp_${SRR}) to keep all I/O on the same filesystem.

# ==============================================================================
# HELPER FUNCTIONS (called by the pipeline)
# ==============================================================================

# Detect actual read length from the second line of the first trimmed FASTQ.
# Returns the integer length; prints nothing on failure so callers can test -z.
# Usage: len=$(_star_detect_read_length "SRR123456")
_star_detect_read_length() {
	local srr="$1"
	find_trimmed_fastq "$srr" 2>/dev/null
	[[ -z "$trimmed1" || ! -f "$trimmed1" ]] && return
	local _seq
	if [[ "$trimmed1" == *.gz ]]; then
		_seq=$($_PIGZ_DC "$trimmed1" 2>/dev/null | sed -n '2p')
	else
		_seq=$(sed -n '2p' "$trimmed1" 2>/dev/null)
	fi
	[[ "${#_seq}" -gt 0 ]] && echo "${#_seq}"
}

# Map STAR_STRAND_SPECIFIC (None/Forward/Reverse) to a Salmon library-type string.
# Salmon's -l A (auto-detect) is used for None so Salmon empirically confirms the
# strandedness rather than accepting a mis-set config value.
# Usage: lib=$(_get_salmon_lib_type "$STAR_STRAND_SPECIFIC" "true|false")
#   paired = "true"  → paired-end prefixes (I*)
#   paired = "false" → single-end prefixes
_get_salmon_lib_type() {
	local strand="${1:-None}"
	local paired="${2:-false}"
	case "$strand" in
		Forward) [[ "$paired" == "true" ]] && echo "ISF" || echo "SF" ;;
		Reverse) [[ "$paired" == "true" ]] && echo "ISR" || echo "SR" ;;
		*)       echo "A" ;;   # auto-detect for unstranded or unknown
	esac
}

# Post-alignment QC: parse Log.final.out files to flag alignment rate outliers,
# abnormal unmapped read categories, and potential rRNA contamination.
# Usage: _star_check_alignment_rates <alignment_dir> <srr1> [srr2 ...]
_star_check_alignment_rates() {
	local align_dir="$1"; shift
	local srr_list=("$@")
	local warn_count=0

	# Collect per-sample stats for cohort-level outlier detection
	local -a sample_names=() input_reads_arr=() unique_rates=() multi_rates=() unmapped_short=() unmapped_mismatch=() unmapped_other=()

	for SRR in "${srr_list[@]}"; do
		local logf="${align_dir}/${SRR}_Log.final.out"
		[[ ! -f "$logf" ]] && continue

		# Extract all key metrics in single AWK pass (replaces 6 separate AWK invocations)
		local uniq_pct multi_pct short_pct mismatch_pct other_pct input_reads
		local _metrics
		_metrics=$(awk -F'|' '
			BEGIN {ir=""; uq=""; ml="0"; sh="0"; mm="0"; ot="0"}
			/Number of input reads/           {gsub(/[[:space:]]/, "", $2); ir=$2}
			/Uniquely mapped reads %/         {gsub(/[[:space:]%]/, "", $2); uq=$2}
			/% of reads mapped to multiple/   {gsub(/[[:space:]%]/, "", $2); ml=$2}
			/% of reads unmapped: too short/  {gsub(/[[:space:]%]/, "", $2); sh=$2}
			/% of reads unmapped: too many/   {gsub(/[[:space:]%]/, "", $2); mm=$2}
			/% of reads unmapped: other/      {gsub(/[[:space:]%]/, "", $2); ot=$2}
			END {printf "%s %s %s %s %s %s", ir, uq, ml, sh, mm, ot}
		' "$logf" 2>/dev/null)
		read -r input_reads uniq_pct multi_pct short_pct mismatch_pct other_pct <<< "$_metrics"

		[[ -z "$uniq_pct" ]] && continue

		sample_names+=("$SRR")
		input_reads_arr+=("${input_reads:-0}")
		unique_rates+=("$uniq_pct")
		multi_rates+=("${multi_pct:-0}")
		unmapped_short+=("${short_pct:-0}")
		unmapped_mismatch+=("${mismatch_pct:-0}")
		unmapped_other+=("${other_pct:-0}")

		# Per-sample checks (truncate to integer for fast bash arithmetic — no awk spawns)
		# 1. Low unique mapping rate (<50% is very concerning, <70% is a warning)
		local _uq_int=${uniq_pct%.*}
		if (( _uq_int < 50 )); then
			log_warn "[STAR QC] $SRR: Uniquely mapped only ${uniq_pct}% — VERY LOW (check sample quality, adapter contamination, or genome mismatch)"
			((warn_count++)) || true
		elif (( _uq_int < 70 )); then
			log_warn "[STAR QC] $SRR: Uniquely mapped ${uniq_pct}% — below 70% threshold"
			((warn_count++)) || true
		fi

		# 2. High multi-mapping (>20%) may indicate rRNA contamination or repetitive sequences
		local _mp_int=${multi_pct:-0}; _mp_int=${_mp_int%.*}
		if (( _mp_int > 20 )); then
			log_warn "[STAR QC] $SRR: Multi-mapped ${multi_pct}% — high rate may indicate rRNA contamination or repetitive element enrichment"
			((warn_count++)) || true
		fi

		# 3. High unmapped-too-short (>15%) suggests adapter contamination or degraded RNA
		local _sp_int=${short_pct:-0}; _sp_int=${_sp_int%.*}
		if (( _sp_int > 15 )); then
			log_warn "[STAR QC] $SRR: Unmapped (too short) ${short_pct}% — check adapter trimming or RNA degradation"
			((warn_count++)) || true
		fi

		# 4. High unmapped-too-many-mismatches (>5%) suggests genome version mismatch
		local _mm_int=${mismatch_pct:-0}; _mm_int=${_mm_int%.*}
		if (( _mm_int > 5 )); then
			log_warn "[STAR QC] $SRR: Unmapped (mismatches) ${mismatch_pct}% — may indicate genome/species mismatch"
			((warn_count++)) || true
		fi
	done

	# Cohort-level outlier detection: flag samples >2 SD below mean unique mapping rate
	local n=${#unique_rates[@]}
	if [[ $n -ge 3 ]]; then
		# O(S) — single AWK pass over S sample rates computes mean, SD, and outlier indices simultaneously
		local _stats_and_outliers
		_stats_and_outliers=$(printf '%s\n' "${unique_rates[@]}" | awk '{
			vals[NR-1]=$1; s+=$1; ss+=$1*$1
		} END {
			m=s/NR; v=ss/NR - m*m; sd=(v>0)?sqrt(v):0; thr=m-2*sd
			printf "%.2f %.2f %.2f", m, sd, thr
			for (i=0; i<NR; i++) if (vals[i]+0 < thr+0) printf " %d", i
		}')
		local mean sd threshold _outlier_rest
		read -r mean sd threshold _outlier_rest <<< "$_stats_and_outliers"

		log_info "[STAR QC] Cohort alignment stats: mean=${mean}%, SD=${sd}%, outlier threshold=${threshold}%"

		# Extract outlier indices (captured by _outlier_rest from read above)
		local _oi
		for _oi in $_outlier_rest; do
			log_warn "[STAR QC] OUTLIER: ${sample_names[$_oi]} (${unique_rates[$_oi]}%) is >2 SD below cohort mean (${mean}%)"
			((warn_count++)) || true
		done
	fi

	# Summary table
	if [[ ${#sample_names[@]} -gt 0 ]]; then
		log_info "[STAR QC] ┌───────────────────┬────────────┬────────┬────────┬──────────┬──────────┬────────┐"
		log_info "[STAR QC] │ Sample            │ Input Rds  │ Unique │ Multi  │ Unmap:Sh │ Unmap:MM │ Unmap:O│"
		log_info "[STAR QC] ├───────────────────┼────────────┼────────┼────────┼──────────┼──────────┼────────┤"
		for ((i=0; i<${#sample_names[@]}; i++)); do
			printf -v _row "[STAR QC] │ %-17s │ %10s │ %5s%% │ %5s%% │   %5s%% │   %5s%% │ %5s%%│" \
				"${sample_names[$i]}" "${input_reads_arr[$i]}" "${unique_rates[$i]}" "${multi_rates[$i]}" \
				"${unmapped_short[$i]}" "${unmapped_mismatch[$i]}" "${unmapped_other[$i]}"
			log_info "$_row"
		done
		log_info "[STAR QC] └───────────────────┴────────────┴────────┴────────┴──────────┴──────────┴────────┘"
	fi

	if [[ $warn_count -gt 0 ]]; then
		log_warn "[STAR QC] $warn_count warning(s) detected — review samples before proceeding"
	else
		log_info "[STAR QC] All samples passed alignment rate checks"
	fi
}

# ==============================================================================
# MAIN STAR ALIGNMENT PIPELINE
# ==============================================================================

star_alignment_pipeline() {
	# Resolve GTF at runtime so the config's gtf_file is available
	STAR_GTF_FILE="${STAR_GTF_FILE:-${gtf_file:-}}"

	local fasta="" transcriptome_fasta="" rnaseq_list=() tissue_tag=""

	# Check for GNU parallel (uses shared cache from shared_utils_method.sh)
	if ! $_SHARED_HAS_PARALLEL; then
		log_warn "[PERFORMANCE] GNU parallel not found - Salmon quantification will run sequentially"
	fi

	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--TRANSCRIPTOME) transcriptome_fasta="$2"; shift 2;;
			--TISSUE_TAG) tissue_tag="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done

	[[ -z "$fasta" ]] && { log_error "No FASTA file specified. Use --FASTA <fasta_file>."; return 1; }
	[[ ! -f "$fasta" ]] && { log_error "FASTA file '$fasta' not found."; return 1; }

	# Resolve transcriptome FASTA for Salmon: explicit arg > config var > auto-detect > fallback to genome
	if [[ -z "$transcriptome_fasta" ]]; then
		transcriptome_fasta="${STAR_TRANSCRIPTOME_FASTA:-}"
	fi
	if [[ -z "$transcriptome_fasta" ]]; then
		# Auto-detect: look for transcriptome FASTA alongside the genome FASTA.
		# Try multiple naming conventions:
		#   1. <basename>_transcripts.fa       (e.g. GPE001970_genome_transcripts.fa)
		#   2. <prefix>_transcripts.fa          (strip _genome suffix: GPE001970_transcripts.fa)
		#   3. <prefix>_transcripts.function.fa (e.g. Eggplant_V4.1_transcripts.function.fa)
		local fasta_dir fasta_stem auto_tx=""
		# Pure bash parameter expansion — avoids 2 subshell forks
		# Edge cases: bare filename → ".", root path "/file.fa" → "/", normal path → dirname
		fasta_dir="${fasta%/*}"
		if [[ "$fasta_dir" == "$fasta" ]]; then
			fasta_dir="."
		elif [[ -z "$fasta_dir" ]]; then
			fasta_dir="/"
		fi
		local _fasta_name="${fasta##*/}"
		fasta_stem="${_fasta_name%.*}"
		local prefix="${fasta_stem%_genome}"  # strip _genome suffix if present
		for candidate in \
			"${fasta_dir}/${fasta_stem}_transcripts.fa" \
			"${fasta_dir}/${prefix}_transcripts.fa" \
			"${fasta_dir}/${prefix}_transcripts.function.fa"; do
			if [[ -f "$candidate" ]]; then
				auto_tx="$candidate"
				break
			fi
		done

		if [[ -n "$auto_tx" ]]; then
			transcriptome_fasta="$auto_tx"
			log_info "[STAR] Auto-detected transcriptome FASTA: $transcriptome_fasta"
		else
			transcriptome_fasta="$fasta"
			log_warn "[STAR] No transcriptome FASTA found. Salmon will index the genome FASTA."
			log_warn "[STAR] Set STAR_TRANSCRIPTOME_FASTA or use --TRANSCRIPTOME to fix this."
		fi
	else
		if [[ ! -f "$transcriptome_fasta" ]]; then
			log_error "Transcriptome FASTA '$transcriptome_fasta' not found."
			return 1
		fi
		log_info "[STAR] Using transcriptome FASTA for Salmon: $transcriptome_fasta"
	fi
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	[[ ${#rnaseq_list[@]} -eq 0 ]] && { log_error "No RNA-seq samples provided."; return 1; }

	local fasta_base fasta_tag star_index_dir star_genome_dir
	# Pure bash parameter expansion — avoids 1 subshell fork vs $(basename ...)
	fasta_base="${fasta##*/}"
	fasta_tag="${fasta_base%.*}"
	set_fasta_output_dirs "$fasta_tag"

	# Resolve absolute paths (STAR_INDEX_ROOT and STAR_ALIGN_ROOT are set by
	# set_fasta_output_dirs → global_config_method.sh which already converts to
	# absolute paths via PROJECT_ROOT/pwd; just guard the edge case).
	local abs_star_index_root="$STAR_INDEX_ROOT"
	local abs_star_align_root="$STAR_ALIGN_ROOT"
	[[ "$abs_star_index_root" != /* ]] && abs_star_index_root="${PROJECT_ROOT:-$(pwd)}/$abs_star_index_root"
	[[ "$abs_star_align_root" != /* ]] && abs_star_align_root="${PROJECT_ROOT:-$(pwd)}/$abs_star_align_root"
	abs_star_index_root="${abs_star_index_root//\/\//\/}"
	abs_star_align_root="${abs_star_align_root//\/\//\/}"

	mkdir -p "$abs_star_index_root" "$abs_star_align_root" 2>/dev/null || true

	# Set directories based on fasta_tag and tissue_tag (using absolute paths).
	# STAR_INDEX_ROOT already contains {fasta_tag} (set by set_fasta_output_dirs).
	# STAR_ALIGN_ROOT does NOT contain {fasta_tag}; we embed it explicitly below to
	# isolate BAM outputs and Salmon quant per reference, preventing multi-reference collisions.
	#
	# The STAR genome index depends only on the genome FASTA + GTF, NOT on tissue.
	# Always use the shared base path so star_tissue_specific_pipeline() never
	# rebuilds the same 30+ GB index for every tissue (analogous to salmon_idx).
	star_index_dir="${abs_star_index_root}/5_star/index"
	if [[ -n "$tissue_tag" ]]; then
		star_genome_dir="${abs_star_align_root}/${fasta_tag}/5_star/alignments/${tissue_tag}"
		log_info "[STAR] Tissue-specific alignment for: $tissue_tag (${#rnaseq_list[@]} samples)"
	else
		star_genome_dir="${abs_star_align_root}/${fasta_tag}/5_star/alignments"
		log_info "[STAR] Pooled alignment: ${#rnaseq_list[@]} samples"
	fi

	# Clean up any double slashes in final paths
	star_index_dir="${star_index_dir//\/\//\/}"
	star_genome_dir="${star_genome_dir//\/\//\/}"

	# Log resolved paths for debugging
	log_info "[STAR] Index directory (absolute): $star_index_dir"
	log_info "[STAR] Output directory (absolute): $star_genome_dir"

	# Auto-detect actual read length from the first sample's trimmed FASTQ.
	# This ensures sjdbOverhang = readLength-1 is always correct regardless of
	# whether STAR_READ_LENGTH was explicitly set (e.g. 150 bp NovaSeq vs 100 bp HiSeq).
	local _det_len
	_det_len=$(_star_detect_read_length "${rnaseq_list[0]}" 2>/dev/null || true)
	if [[ -n "$_det_len" && "$_det_len" -gt 20 ]]; then
		if [[ "$_det_len" -ne "$STAR_READ_LENGTH" ]]; then
			log_warn "[STAR] Detected read length (${_det_len} bp) differs from STAR_READ_LENGTH=${STAR_READ_LENGTH}"
			log_info "[STAR] Using detected length for optimal sjdbOverhang. Set STAR_READ_LENGTH after auto-detect to override."
			STAR_READ_LENGTH="$_det_len"
		else
			log_info "[STAR] Read length confirmed by FASTQ: ${STAR_READ_LENGTH} bp"
		fi
	else
		log_info "[STAR] Using configured read length: ${STAR_READ_LENGTH} bp (FASTQ not yet available for detection)"
	fi

	local star_overhang=$((STAR_READ_LENGTH - 1))
	log_info "[STAR] sjdbOverhang = $star_overhang (read length ${STAR_READ_LENGTH} bp)"
	log_info "[STAR] CPU allocation: Total=$THREADS threads"

	# Pre-compute Salmon library type strings from strandedness config.
	# Computed once here so all per-sample calls (sequential + parallel) are consistent.
	local _sal_lib_pe _sal_lib_se
	_sal_lib_pe=$(_get_salmon_lib_type "${STAR_STRAND_SPECIFIC:-None}" "true")
	_sal_lib_se=$(_get_salmon_lib_type "${STAR_STRAND_SPECIFIC:-None}" "false")
	log_info "[STAR] Salmon library type: PE=${_sal_lib_pe}  SE=${_sal_lib_se}  (STAR_STRAND_SPECIFIC=${STAR_STRAND_SPECIFIC:-None})"

	# Pre-compute STAR strandedness and splice-junction args (used in both parallel + sequential).
	# --outSAMstrandField intronMotif: infers XS strand tag from splice site motifs.
	#   Appropriate for UNSTRANDED libraries; for stranded protocols STAR derives XS from
	#   the read orientation, so adding intronMotif would override with a less-reliable signal.
	local star_strand_args=()
	if [[ "${STAR_STRAND_SPECIFIC:-None}" == "None" ]]; then
		star_strand_args+=(--outSAMstrandField intronMotif)
		log_info "[STAR] Unstranded: using --outSAMstrandField intronMotif for XS tags"
	else
		log_info "[STAR] Stranded (${STAR_STRAND_SPECIFIC}): XS tag derived from read orientation"
	fi

	# Explicit splice junction filtering thresholds (STAR defaults documented here
	# so the pipeline behaviour is transparent across STAR versions).
	# Format: 4 values for canonical GT/AG, semi-canonical CT/AC or GT/AT, non-canonical, other.
	# --outSJfilterCountUniqueMin   3 1 1 1   (min unique reads per junction type)
	# --outSJfilterCountTotalMin    3 1 1 1   (min total reads per junction type)
	# --outSJfilterOverhangMin      30 12 12 12 (min overhang for reported junctions)
	# --outSJfilterIntronMaxVsReadN 50000 100000 200000 (max intron vs read length)
	local star_sj_filter_args=(
		--outSJfilterCountUniqueMin 3 1 1 1
		--outSJfilterCountTotalMin  3 1 1 1
		--outSJfilterOverhangMin    30 12 12 12
		--outSJfilterIntronMaxVsReadN 50000 100000 200000
	)

	# Resolve effective genome load. --twopassMode Basic is incompatible with LoadAndKeep;
	# override and warn rather than silently degrading to 1-pass alignment.
	# effective_genome_load is not declared local so it can be exported for parallel workers.
	# OPTIMIZATION NOTE: For multi-sample sequential runs against the same genome,
	# LoadAndRemove on the last sample (or LoadAndKeep + STAR --genomeRemove at cleanup)
	# can save 15-25 min of genome loading per reference. Requires restructuring the
	# alignment loop to track first/last sample, so deferred for now.
	effective_genome_load="${STAR_GENOME_LOAD:-NoSharedMemory}"
	if [[ "$effective_genome_load" == "LoadAndKeep" ]]; then
		log_warn "[STAR] STAR_GENOME_LOAD=LoadAndKeep is incompatible with --twopassMode Basic. Overriding to NoSharedMemory."
		effective_genome_load="NoSharedMemory"
	fi
	log_info "[STAR] Genome load mode: $effective_genome_load"

	# Validate GTF early - used by both index build and tx2gene creation.
	# Do this outside the index-skip branch so the check runs even when the index already exists.
	if [[ ! -f "${STAR_GTF_FILE:-}" ]]; then
		log_error "[STAR] GTF annotation file not found: '${STAR_GTF_FILE:-<unset>}'"
		log_error "[STAR] Set STAR_GTF_FILE (or gtf_file) in your configuration before running M3."
		return 1
	fi
	log_info "[STAR] GTF annotation: $STAR_GTF_FILE"

	# STEP 1: BUILD STAR GENOME INDEX
	log_step "STAR genome index generation for $fasta_tag"

	if [[ -f "$star_index_dir/SAindex" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[STAR INDEX] Index for $fasta_tag already exists. Skipping."
	else
		log_info "[STAR INDEX] Building STAR genome index from genome FASTA..."
		mkdir -p "$star_index_dir"

		log_file_size "$fasta" "Input genome FASTA for STAR indexing"

		# Compute recommended SAindex size: min(14, floor(log2(GenomeLength)/2 - 1))
		local genome_sa_index=14
		local genome_chr_bin=18
		local genome_bp num_seqs
		# Single awk pass: extract genome metrics AND compute STAR index parameters
		# (1 process instead of 3 — eliminates two extra awk BEGIN invocations)
		eval "$(awk '
			/^>/{n++} !/^>/{bp+=length($0)}
			END{
				if(n<1) n=1; if(bp<1) bp=0
				sa=14; cb=18
				if(bp>0){
					v=int(log(bp)/log(2)/2-1); sa=(v<14)?v:14; if(sa<1) sa=1
					v=int(log(bp/n)/log(2)-1);  cb=(v<18)?v:18; if(cb<1) cb=1
				}
				printf "genome_bp=%d num_seqs=%d genome_sa_index=%d genome_chr_bin=%d", bp, n, sa, cb
			}' "$fasta" 2>/dev/null)"
		log_info "[STAR INDEX] Genome size: ${genome_bp:-unknown} bp, $num_seqs sequences -> genomeSAindexNbases=$genome_sa_index, genomeChrBinNbits=$genome_chr_bin"

		# O(G_bp × log(G_bp)) — STAR suffix-array construction scales super-linearly with genome size;
		# suffix array sort dominates at ~15 min for a 1 Gbp genome; run once per reference
		run_with_space_time_log --input "$fasta" --output "$star_index_dir" \
			STAR --runMode genomeGenerate \
				--genomeDir "$star_index_dir" \
				--genomeFastaFiles "$fasta" \
				--sjdbGTFfile "$STAR_GTF_FILE" \
				--sjdbOverhang "$star_overhang" \
				--genomeSAindexNbases "$genome_sa_index" \
				--genomeChrBinNbits "$genome_chr_bin" \
				--runThreadN "$THREADS"

		[[ ! -f "$star_index_dir/SAindex" ]] && { log_error "STAR index generation failed"; return 1; }
		log_info "[STAR INDEX] Index built successfully: $star_index_dir"
	fi

	# STEP 2: ALIGN READS WITH STAR
	mkdir -p "$star_genome_dir"

	# Verify directory exists and is writable
	if [[ ! -d "$star_genome_dir" ]]; then
		log_error "Failed to create STAR output directory: $star_genome_dir"
		return 1
	fi
	if [[ ! -w "$star_genome_dir" ]]; then
		log_error "STAR output directory is not writable: $star_genome_dir"
		return 1
	fi
	log_info "[STAR] Output directory verified: $star_genome_dir"

	log_step "STAR splice-aware alignment for $fasta_tag samples"

	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	# Adaptive: cap parallel jobs at sample count to maximize per-job thread allocation
	local _n_samples=${#rnaseq_list[@]}
	(( parallel_jobs > _n_samples )) && parallel_jobs=$_n_samples
	(( parallel_jobs < 1 )) && parallel_jobs=1
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	# Calculate --limitBAMsortRAM for STAR's internal coordinate sorting.
	# STAR sorts in-process when --outSAMtype BAM SortedByCoordinate, avoiding
	# a separate samtools sort step (saves one full BAM read+write pass).
	local _star_sort_ram_bytes
	_star_sort_ram_bytes=$(_star_sort_ram "$parallel_jobs")
	export _star_sort_ram_bytes

	# Use cached parallel availability (set at module load) to avoid per-invocation command -v
	if $_SHARED_HAS_PARALLEL && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
		log_step "[PARALLEL] STAR alignment: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		log_warn "[PARALLEL] STAR is memory-intensive (~30GB/instance). Ensure sufficient RAM for $parallel_jobs concurrent jobs."
		_prepare_parallel_env

		export fasta_tag threads_per_job
		export star_index_dir star_genome_dir
		export STAR_DELETE_TRANSIENT PROJECT_ROOT
		export effective_genome_load
		export _STAR_READ_CMD _star_sort_ram_bytes
		# Serialize array args for parallel workers (bash can't export arrays)
		export _star_strand_args_str="${star_strand_args[*]}"
		export _star_sj_filter_args_str="${star_sj_filter_args[*]}"

		_m3_star_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log STAR "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"

			# Check if BAM exists and has content (shared helper: _bam_is_valid)
			if [[ -f "$bam_output" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				local bam_size
				bam_size=$(_bam_is_valid "$bam_output" 1000) && {
					_parallel_log STAR "$SRR" INFO "BAM exists (${bam_size} bytes) - skipping"
					return 0
				}
				rm -f "$bam_output"
			fi

			# Clean up stale files from previous runs
			rm -f "${star_genome_dir}/${SRR}_Log.out" "${star_genome_dir}/${SRR}_Log.progress.out" \
				"${star_genome_dir}/${SRR}_Log.final.out" "${star_genome_dir}/${SRR}_SJ.out.tab" 2>/dev/null || true
			rm -rf "${star_genome_dir}/${SRR}__STARgenome" "${star_genome_dir}/${SRR}__STARpass1" \
				"${star_genome_dir}/${SRR}_STARtmp" "${star_genome_dir}/_STARtmp_${SRR}" 2>/dev/null || true

			local star_reads_args=("$trimmed1")
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				star_reads_args+=("$trimmed2")
			fi

			# Only use --readFilesCommand for compressed input; STAR reads uncompressed natively
			local _par_read_cmd_args=()
			if [[ "$trimmed1" == *.gz ]]; then
				_par_read_cmd_args=(--readFilesCommand "$_STAR_READ_CMD")
			fi

			local star_tmp_dir="${star_genome_dir}/_STARtmp_${SRR}"
			rm -rf "$star_tmp_dir" 2>/dev/null || true

			local out_prefix="${star_genome_dir}/${SRR}_"
			out_prefix="${out_prefix//\/\//\/}"

			_parallel_log STAR "$SRR" INFO "Aligning with $threads_per_job threads (STAR internal sort, limitBAMsortRAM=${_star_sort_ram_bytes})"
			_parallel_log STAR "$SRR" INFO "--- BEGIN STAR OUTPUT ---"

			# Deserialize array args from exported strings
			local _par_strand_args=($_star_strand_args_str)
			local _par_sj_args=($_star_sj_filter_args_str)

			STAR --runMode alignReads \
				--genomeDir "$star_index_dir" \
				--readFilesIn "${star_reads_args[@]}" \
				${_par_read_cmd_args[@]:+"${_par_read_cmd_args[@]}"} \
				--outFileNamePrefix "$out_prefix" \
				--outTmpDir "$star_tmp_dir" \
				--outSAMtype BAM SortedByCoordinate \
				--limitBAMsortRAM "$_star_sort_ram_bytes" \
				${_par_strand_args[@]:+"${_par_strand_args[@]}"} \
				--outSAMattributes NH HI AS NM MD \
				--outSAMunmapped Within \
				--twopassMode Basic \
				${_par_sj_args[@]:+"${_par_sj_args[@]}"} \
				--genomeLoad "$effective_genome_load" \
				--runThreadN "$threads_per_job" 2>&1 || \
				{ _parallel_log STAR "$SRR" ERROR "--- END STAR OUTPUT (FAILED) ---"; _parallel_log STAR "$SRR" ERROR "STAR alignment failed"; return 1; }

			_parallel_log STAR "$SRR" INFO "--- END STAR OUTPUT ---"

			if [[ ! -f "$bam_output" ]]; then
				_parallel_log STAR "$SRR" ERROR "Sorted BAM not created by STAR"
				return 1
			fi

			local final_bam_size
			final_bam_size=$(_bam_is_valid "$bam_output" 1000) || {
				_parallel_log STAR "$SRR" ERROR "Sorted BAM is empty/corrupt (${final_bam_size} bytes)"
				return 1
			}

			# Index threads: cap at 4 — samtools index is I/O-bound, extra threads add overhead
			local _idx_threads=$threads_per_job
			(( _idx_threads > 4 )) && _idx_threads=4
			samtools index -@ "$_idx_threads" "$bam_output" 2>&1 || true

			# Clean up transient files
			rm -rf "$star_tmp_dir" 2>/dev/null || true
			if [[ "${STAR_DELETE_TRANSIENT:-true}" == "true" ]]; then
				rm -rf "${star_genome_dir}/${SRR}__STARgenome" "${star_genome_dir}/${SRR}__STARpass1" \
					"${star_genome_dir}/${SRR}_STARtmp" "${star_genome_dir}/_STARtmp_${SRR}" 2>/dev/null || true
				rm -f "${star_genome_dir}/${SRR}_Log.progress.out" 2>/dev/null || true
			fi

			_parallel_log STAR "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m3_star_parallel_worker _bam_is_valid

		# O(S/parallel_jobs × (N log N + G_bp)) — S samples dispatched across parallel_jobs slots;
		# each worker runs STAR 2-pass O(N log N) alignment + in-process BAM sort O(N log N)
		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env threads_per_job --env parallel_jobs \
			--env star_index_dir --env star_genome_dir \
			--env STAR_DELETE_TRANSIENT --env PROJECT_ROOT \
			--env OVERWRITE_MODE --env effective_genome_load \
			--env _star_strand_args_str --env _star_sj_filter_args_str \
			--env _STAR_READ_CMD --env _star_sort_ram_bytes \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$star_genome_dir/parallel_star_align.log" \
			_m3_star_parallel_worker {}

		local par_exit=$?
		log_info "[PARALLEL] STAR alignment complete (exit=$par_exit)"
		[[ $par_exit -ne 0 ]] && log_warn "[PARALLEL] Some jobs failed - check $star_genome_dir/parallel_star_align.log"
	else
		# Sequential fallback — recompute RAM budget for a single concurrent job
		# (the initial _star_sort_ram was computed for parallel_jobs instances)
		_star_sort_ram_bytes=$(_star_sort_ram 1)

		# Pre-flight checks: disk space and write permissions (once, not per sample)
		mkdir -p "$star_genome_dir"
		if [[ "${STAR_NFS_SYNC:-false}" == "true" ]]; then
			sync
			sleep 1
		fi
		if [[ ! -d "$star_genome_dir" ]]; then
			log_error "[STAR] FATAL: Cannot create output directory: $star_genome_dir"
			return 1
		fi
		local _preflight_test="${star_genome_dir}/_star_write_test_$$"
		if ! touch "$_preflight_test" 2>/dev/null; then
			log_error "[STAR] FATAL: Output directory not writable: $star_genome_dir"
			return 1
		fi
		rm -f "$_preflight_test"
		local _preflight_space
		_preflight_space=$(df -P "$star_genome_dir" 2>/dev/null | awk 'NR==2 {print $4}')
		if [[ -n "$_preflight_space" ]]; then
			local _preflight_gb=$((_preflight_space / 1024 / 1024))
			if [[ $_preflight_gb -lt 30 ]]; then
				log_error "[STAR] FATAL: Less than 30GB available in $star_genome_dir - STAR needs ~30GB per run"
				return 1
			fi
		fi
		# Clean up any stray project-root temp dirs once
		rm -rf "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true

		for SRR in "${rnaseq_list[@]}"; do
			local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"

			# Check if BAM exists AND has content (shared helper: _bam_is_valid)
			if [[ -f "$bam_output" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				local bam_size
				bam_size=$(_bam_is_valid "$bam_output" 1000) && {
					log_info "[STAR] Alignment for $SRR already exists (${bam_size} bytes). Skipping."
					continue
				}
				log_warn "[STAR] Found empty/corrupt BAM for $SRR (${bam_size} bytes) - removing and re-running"
				rm -f "$bam_output"
			fi

			# Clean up stale files from previous failed STAR runs — single find pass
			# replaces 2 separate rm invocations (O(1) find vs O(N) rm calls)
			log_info "[STAR] Cleaning up stale files for $SRR..."
			find "$star_genome_dir" -maxdepth 1 \( \
				-name "${SRR}_Log.out" -o -name "${SRR}_Log.progress.out" \
				-o -name "${SRR}_Log.final.out" -o -name "${SRR}_SJ.out.tab" \
				-o -name "${SRR}__STARgenome" -o -name "${SRR}__STARpass1" \
				-o -name "${SRR}_STARtmp" -o -name "${SRR}_*.tmp" \
				-o -name "_STARtmp_${SRR}" \
			\) -exec rm -rf {} + 2>/dev/null || true

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ for $SRR not found; skipping."; continue; }

			log_info "[STAR] Aligning: $SRR"

			local star_reads_args=("$trimmed1")
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				star_reads_args+=("$trimmed2")
				log_info "[STAR] Processing paired-end reads for $SRR"
			else
				log_info "[STAR] Processing single-end reads for $SRR"
			fi

			# Only use --readFilesCommand for compressed input; STAR reads uncompressed natively
			local _seq_read_cmd_args=()
			if [[ "$trimmed1" == *.gz ]]; then
				_seq_read_cmd_args=(--readFilesCommand "$_STAR_READ_CMD")
			fi

			# Configure STAR temp directory - ALWAYS use output directory for temp
			# This ensures all STAR operations happen on the same filesystem
			local star_tmp_dir="${star_genome_dir}/_STARtmp_${SRR}"
			rm -rf "$star_tmp_dir" 2>/dev/null || true
			rm -rf "${PROJECT_ROOT}/_STARtmp_${SRR}" 2>/dev/null || true
			log_info "[STAR] Using temp directory: $star_tmp_dir"

			# Construct output prefix - ensure no double slashes and path is clean
			local out_prefix="${star_genome_dir}/${SRR}_"
			# Remove any double slashes
			out_prefix="${out_prefix//\/\//\/}"

			# Run STAR alignment - use internal coordinate sorting (SortedByCoordinate)
			# to avoid a separate samtools sort pass (saves one full BAM read+write cycle).
			# This matches the parallel mode and saves ~15-25 min per sample.
			log_info "[STAR] Aligning with STAR internal sort (limitBAMsortRAM=${_star_sort_ram_bytes})"

			run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$star_genome_dir" \
				STAR --runMode alignReads \
					--genomeDir "$star_index_dir" \
					--readFilesIn "${star_reads_args[@]}" \
					${_seq_read_cmd_args[@]:+"${_seq_read_cmd_args[@]}"} \
					--outFileNamePrefix "$out_prefix" \
					--outTmpDir "$star_tmp_dir" \
					--outSAMtype BAM SortedByCoordinate \
					--limitBAMsortRAM "$_star_sort_ram_bytes" \
					${star_strand_args[@]:+"${star_strand_args[@]}"} \
					--outSAMattributes NH HI AS NM MD \
					--outSAMunmapped Within \
					--twopassMode Basic \
					${star_sj_filter_args[@]:+"${star_sj_filter_args[@]}"} \
					--genomeLoad "$effective_genome_load" \
					--runThreadN "$THREADS"

			# Check if sorted BAM was created — cache log tail to avoid duplicate reads
			local _star_log_tail=""
			if [[ ! -f "$bam_output" ]]; then
				log_error "[STAR] FATAL: Sorted BAM file not created for $SRR"
				log_error "[STAR] Check STAR log: ${out_prefix}Log.out"
				_star_log_tail=$(tail -30 "${out_prefix}Log.out" 2>/dev/null)
				log_error "[STAR] Last 30 lines of STAR log:"
				log_error "$_star_log_tail"
				return 1
			fi

			local final_bam_size
			final_bam_size=$(_bam_is_valid "$bam_output" 1000) || {
				log_error "[STAR] FATAL: Sorted BAM is empty/corrupt for $SRR (${final_bam_size} bytes)"
				[[ -z "$_star_log_tail" ]] && _star_log_tail=$(tail -30 "${out_prefix}Log.out" 2>/dev/null)
				log_error "[STAR] Last 30 lines of STAR log:"
				log_error "$_star_log_tail"
				return 1
			}

			log_info "[STAR] BAM sorted successfully: $final_bam_size bytes"

			# Index the BAM (cap threads at 4 — samtools index is I/O-bound)
			log_info "[STAR] Indexing BAM..."
			local _idx_threads=$THREADS
			(( _idx_threads > 4 )) && _idx_threads=4
			samtools index -@ "$_idx_threads" "$bam_output" 2>&1 || log_warn "[STAR] BAM indexing failed (non-fatal)"

			# Clean up temp directories after successful alignment
			rm -rf "$star_tmp_dir" "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true

			# Delete transient big files if enabled (saves significant disk space)
			if [[ "${STAR_DELETE_TRANSIENT:-true}" == "true" ]]; then
				log_info "[STAR] Cleaning up transient files for $SRR..."

				# Remove 2-pass intermediate dirs, leftover temp dirs, and progress log in one call
				rm -rf "${star_genome_dir}/${SRR}__STARgenome" \
					"${star_genome_dir}/${SRR}__STARpass1" \
					"${star_genome_dir}/${SRR}_STARtmp" \
					"${star_genome_dir}/_STARtmp_${SRR}" \
					"${star_genome_dir}/${SRR}_Log.progress.out" 2>/dev/null || true

				log_info "[STAR] Transient files cleaned up for $SRR"
			fi

			log_info "[STAR] Successfully aligned: $SRR"
		done
	fi

	log_info "[STAR] All samples aligned successfully"

	# POST-ALIGNMENT QC: Parse Log.final.out for alignment rate outliers and anomalies
	log_step "STAR post-alignment quality check"
	_star_check_alignment_rates "$star_genome_dir" "${rnaseq_list[@]}"

	# STEP 3: SALMON QUANTIFICATION
	# Include fasta_tag in paths to prevent multi-reference collisions.
	# tissue_tag further subdivides quant outputs within a given reference run.
	local ref_suffix=""
	[[ -n "${tissue_tag:-}" ]] && ref_suffix="${tissue_tag}"
	# The Salmon index depends only on the transcriptome FASTA, not on tissue.
	# Do NOT include ref_suffix here; otherwise star_tissue_specific_pipeline()
	# rebuilds the same expensive index once per tissue.
	local salmon_idx="${abs_star_align_root}/${fasta_tag}/6_salmon/index"
	local quant_root="${abs_star_align_root}/${fasta_tag}/6_salmon/quant${ref_suffix:+/${ref_suffix}}"
	# Clean up double slashes
	salmon_idx="${salmon_idx//\/\//\/}"
	quant_root="${quant_root//\/\//\/}"
	mkdir -p "$quant_root"

	# Build Salmon index (using transcriptome FASTA, not genome)
	if [[ -f "$salmon_idx/versionInfo.json" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[SALMON INDEX] STAR Salmon index exists - skipping build"
	else
		log_step "Building Salmon index for transcriptome"
		log_info "[SALMON INDEX] Indexing: $transcriptome_fasta"
		run_with_space_time_log --input "$transcriptome_fasta" --output "$salmon_idx" \
			salmon index -t "$transcriptome_fasta" -i "$salmon_idx" -k 31 --threads "$THREADS"
		[[ ! -f "$salmon_idx/versionInfo.json" ]] && { log_error "[SALMON INDEX] Index build failed - versionInfo.json not found in $salmon_idx"; return 1; }
		log_info "[SALMON INDEX] Index built successfully: $salmon_idx"
	fi

	# Validate transcriptome FASTA IDs match GTF transcript IDs (prevents tximport failures)
	# Cache result per index directory to skip expensive GTF/FASTA parsing on re-runs
	local _tx_validation_cache="$salmon_idx/.tx_id_validated"
	if [[ -f "$_tx_validation_cache" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[SALMON] Transcript ID validation: using cached result"
	elif [[ -f "${STAR_GTF_FILE:-}" && -f "$transcriptome_fasta" ]]; then
		# Single awk pass over both files: extract first 20 FASTA IDs, build GTF transcript
		# ID set, compute overlap (replaces grep|head|sed + awk|sort -u + grep -cFxf chain)
		local _tx_result
		_tx_result=$(awk '
			FILENAME == ARGV[1] && /^>/ && fn < 20 {
				id = $0; sub(/^>/, "", id); sub(/ .*/, "", id)
				fasta_ids[fn] = id; fn++
			}
			FILENAME == ARGV[2] && $3 == "transcript" {
				for (i=9; i<=NF; i++) if ($i == "transcript_id") {
					id = $(i+1); gsub(/[";]/, "", id); gtf[id] = 1
				}
			}
			END {
				overlap = 0
				for (i=0; i<fn; i++) if (fasta_ids[i] in gtf) overlap++
				# Format: overlap fasta_count fasta_first3 | gtf_first3
				printf "%d %d ", overlap, fn
				for (i=0; i<fn && i<3; i++) printf "%s,", fasta_ids[i]
				printf "| "
				n=0; for (id in gtf) { if (n<3) printf "%s,", id; n++; if (n>=3) break }
			}
		' "$transcriptome_fasta" "$STAR_GTF_FILE" 2>/dev/null)
		local _overlap _fasta_count _fasta_ex _gtf_ex
		read -r _overlap _fasta_count _fasta_ex _ _gtf_ex <<< "${_tx_result//|/ | }"
		if [[ "${_overlap:-0}" -eq 0 && "${_fasta_count:-0}" -gt 0 ]]; then
			log_warn "[SALMON] Transcript ID mismatch: transcriptome FASTA IDs do not match GTF transcript_id attributes"
			log_warn "[SALMON]   FASTA example: ${_fasta_ex:-N/A}"
			log_warn "[SALMON]   GTF example:   ${_gtf_ex:-N/A}"
			log_warn "[SALMON]   Gene-level tximport will FAIL. Isoform-level will still work."
			log_warn "[SALMON]   Fix: use a transcriptome FASTA derived from the same annotation as the GTF,"
			log_warn "[SALMON]   or set STAR_TRANSCRIPTOME_FASTA to a FASTA whose IDs match: ${STAR_GTF_FILE}"
		fi
		touch "$_tx_validation_cache"
	fi

	# Quantify samples
	log_info "[SALMON] Starting quantification for ${#rnaseq_list[@]} samples"

	# Use cached parallel availability (set at module load) to avoid per-invocation command -v
	if $_SHARED_HAS_PARALLEL && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
		log_step "[PARALLEL] Salmon quant (STAR): ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export fasta_tag threads_per_job salmon_idx quant_root
		export _sal_lib_pe _sal_lib_se

		_m3_salmon_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log SALMON_STAR "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local quant_dir="$quant_root/$SRR"
			[[ -f "$quant_dir/quant.sf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]] && { _parallel_log SALMON_STAR "$SRR" INFO "quant.sf exists - skipping"; return 0; }

			mkdir -p "$quant_dir"
			_parallel_log SALMON_STAR "$SRR" INFO "Quantifying with $threads_per_job threads"

			local quant_exit=0
			# Redirect to log file — avoids sed pipe (saves 1 process per sample) and
			# simplifies exit code capture (direct $? vs PIPESTATUS).
			# Parallel already captures worker stdout, so piping through sed doubled I/O.
			local _salmon_log="$quant_dir/salmon_quant.log"
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				salmon quant -p "$threads_per_job" -i "$salmon_idx" -o "$quant_dir" \
					--gcBias --seqBias -l "${_sal_lib_pe:-A}" -1 "$trimmed1" -2 "$trimmed2" \
					> "$_salmon_log" 2>&1
				quant_exit=$?
			else
				salmon quant -p "$threads_per_job" -i "$salmon_idx" -o "$quant_dir" \
					--gcBias --seqBias -l "${_sal_lib_se:-A}" -r "$trimmed1" \
					> "$_salmon_log" 2>&1
				quant_exit=$?
			fi
			[[ $quant_exit -ne 0 ]] && { _parallel_log SALMON_STAR "$SRR" ERROR "Salmon quant failed (exit=$quant_exit)"; return $quant_exit; }

			_parallel_log SALMON_STAR "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m3_salmon_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env threads_per_job --env salmon_idx --env quant_root \
			--env _sal_lib_pe --env _sal_lib_se \
			--env OVERWRITE_MODE \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$quant_root/parallel_salmon_star_quant.log" \
			_m3_salmon_parallel_worker {}

		local par_exit_salmon=$?
		log_info "[PARALLEL] Salmon quant complete (exit=$par_exit_salmon)"
		[[ $par_exit_salmon -ne 0 ]] && log_warn "[PARALLEL] Some Salmon jobs failed - check $quant_root/parallel_salmon_star_quant.log"
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local quant_dir="$quant_root/$SRR"

			[[ -f "$quant_dir/quant.sf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]] && { log_info "[SALMON] Quantification for $SRR exists - skipping"; continue; }

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "No trimmed reads for $SRR - skipping"; continue; }

			log_info "[SALMON] Quantifying: $SRR"
			mkdir -p "$quant_dir"

			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				run_with_space_time_log salmon quant -p "$THREADS" -i "$salmon_idx" -o "$quant_dir" \
					--gcBias --seqBias -l "${_sal_lib_pe:-A}" -1 "$trimmed1" -2 "$trimmed2"
			else
				run_with_space_time_log salmon quant -p "$THREADS" -i "$salmon_idx" -o "$quant_dir" \
					--gcBias --seqBias -l "${_sal_lib_se:-A}" -r "$trimmed1"
			fi

			if [[ -f "$quant_dir/quant.sf" ]]; then
				log_info "[SALMON] Successfully quantified: $SRR"
			else
				log_warn "[SALMON] quant.sf not produced for $SRR - check Salmon output above"
			fi
		done
	fi

	# STEP 4: PREPARE TXIMPORT FILES FOR DESEQ2
	log_step "Preparing tximport input for DESeq2 (STAR + Salmon)"

	# Include tissue/ref suffix so tissue-specific runs get isolated matrix dirs.
	# Without this, star_tissue_specific_pipeline() tissue runs overwrite each
	# other's sample_info.tsv, tximport script, and count matrix TSVs.
	local matrix_dir="${STAR_MATRIX_ROOT}${ref_suffix:+/${ref_suffix}}"
	# master_ref drives output filenames in tximport_star_helper.R; include the
	# tissue suffix so per-tissue TSV files have distinct, non-colliding names.
	local master_ref="${fasta_tag}${ref_suffix:+_${ref_suffix}}"
	mkdir -p "$matrix_dir"

	local sample_metadata="$matrix_dir/sample_info.tsv"
	# tx2gene is derived from the GTF (same for all tissues of this reference);
	# keep it in the base STAR_MATRIX_ROOT so it is shared and not re-generated
	# for every tissue.  Name it only by fasta_tag.
	local tx2gene_file="${STAR_MATRIX_ROOT}/tx2gene_${fasta_tag}.tsv"
	mkdir -p "$STAR_MATRIX_ROOT"

	# Create tx2gene mapping from GTF (transcript_id -> gene_id attributes)
	# This correctly handles multi-transcript genes; version-stripping FASTA headers is not reliable.
	if [[ ! -f "$tx2gene_file" || "${OVERWRITE_MODE:-skip}" == "overwrite" ]]; then
		log_info "[TXIMPORT] Creating transcript-to-gene mapping from GTF: $STAR_GTF_FILE"
		awk '$3=="transcript" {
			tid=""; gid=""
			for (i=9; i<=NF; i++) {
				if ($i == "transcript_id") { gsub(/[";]/, "", $(i+1)); tid=$(i+1) }
				if ($i == "gene_id")       { gsub(/[";]/, "", $(i+1)); gid=$(i+1) }
			}
			if (tid != "" && gid != "") print tid "\t" gid
		}' "$STAR_GTF_FILE" | awk '!seen[$0]++' > "$tx2gene_file"
		# O(n) awk dedup vs O(n log n) sort -u; awk END{print NR} vs wc -l subprocess
		local tx2gene_count
		tx2gene_count=$(awk 'END{print NR}' "$tx2gene_file")
		if [[ "$tx2gene_count" -eq 0 ]]; then
			log_error "[TXIMPORT] tx2gene mapping is empty - check GTF has 'transcript' features with transcript_id/gene_id attributes"
			return 1
		fi
		log_info "[TXIMPORT] Created tx2gene mapping: $tx2gene_count transcripts"
	fi

	# Collect samples that actually have quant.sf.
	# Only passing successful samples to create_sample_metadata prevents tximport_star_helper.R
	# from stop()-ing on the first missing quant.sf when partial failures occurred.
	local quant_srrs=()
	for SRR in "${rnaseq_list[@]}"; do
		[[ -f "$quant_root/$SRR/quant.sf" ]] && quant_srrs+=("$SRR")
	done

	[[ ${#quant_srrs[@]} -lt 2 ]] && { log_error "Insufficient quantifications: ${#quant_srrs[@]} (need ≥2)"; return 1; }

	# Create sample metadata from successfully quantified samples only.
	create_sample_metadata "$sample_metadata" "${quant_srrs[@]}"

	# Generate tximport R script (always refresh so helper updates propagate)
	local tximport_script="$matrix_dir/run_tximport_star_salmon.R"
	generate_tximport_star_script "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" "$tximport_script" || return 1

	# Run tximport if R is available
	if command -v Rscript >/dev/null 2>&1; then
		log_step "Running tximport to import Salmon quantifications"
		# Note: stdout already goes through tee via exec redirect; do NOT pipe to tee -a "$LOG_FILE" (causes double-logging)
		if Rscript "$tximport_script" "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" "$master_ref" 2>&1; then
			log_info "[TXIMPORT] Successfully imported counts for DESeq2"
			local count_matrix="$matrix_dir/gene_level/${master_ref}_NumReads_Gene_ID_from_${master_ref}_gene_level.tsv"
			[[ -f "$count_matrix" ]] && validate_count_matrix "$count_matrix" "gene" 2
		else
			log_warn "[TXIMPORT] tximport failed - see R error output above for details"
			log_warn "[TXIMPORT] Common causes: transcript ID mismatch between quant.sf and tx2gene, or missing R packages (BiocManager::install('tximport'))"
		fi
	else
		log_warn "[TXIMPORT] Rscript not found - run manually: Rscript $tximport_script"
	fi

	log_step "STAR alignment pipeline completed for $fasta_tag"
	log_info "STAR alignments: $star_genome_dir"
	log_info "Salmon quantifications: $quant_root"
	log_info "Tximport outputs: $matrix_dir"
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Run tximport for STAR+Salmon using external R helper
# Usage: run_tximport_star <quant_dir> <metadata_file> <tx2gene_file> [output_dir] [master_ref]
run_tximport_star() {
	local quant_dir="$1"
	local metadata_file="$2"
	local tx2gene_file="$3"
	local output_dir="${4:-$(dirname "$metadata_file")}"
	local master_ref="${5:-$(basename "$output_dir")}"
	local helper_script="$_M3_SCRIPT_DIR/../c_post_processing/preprocessing/STAR/tximport_star_helper.R"

	if [[ ! -f "$helper_script" ]]; then
		log_error "tximport_star_helper.R not found: $helper_script"
		return 1
	fi

	log_info "[TXIMPORT] Running STAR+Salmon import..."
	Rscript "$helper_script" "$quant_dir" "$metadata_file" "$tx2gene_file" "$output_dir" "$master_ref"
}

# Generate tximport script (copies helper to output location)
# Usage: generate_tximport_star_script <output_script>
# Legacy callers pass extra args (quant_root, metadata, tx2gene, matrix_dir)
# before output_script — accept and ignore them for backward compatibility.
generate_tximport_star_script() {
	local output_script="${!#}"  # last argument
	local helper_script="$_M3_SCRIPT_DIR/../c_post_processing/preprocessing/STAR/tximport_star_helper.R"

	if [[ -f "$helper_script" ]]; then
		cp "$helper_script" "$output_script"
		chmod +x "$output_script"
		log_info "[TXIMPORT] Copied helper to: $output_script"
	else
		log_error "tximport_star_helper.R not found: $helper_script"
		return 1
	fi
}

# ==============================================================================
# ALTERNATIVE FLOW: TISSUE-SPECIFIC ALIGNMENT
# ==============================================================================
# Splits samples by condition/tissue and runs star_alignment_pipeline() per group.
# Falls back to pooled alignment when only one tissue type is detected.

star_tissue_specific_pipeline() {
	local fasta="" rnaseq_list=() metadata_file=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--METADATA) metadata_file="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done

	[[ -z "$fasta" || ! -f "$fasta" ]] && { log_error "Valid FASTA required"; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")

	declare -A sample_metadata
	if [[ -n "$metadata_file" && -f "$metadata_file" ]]; then
		load_sample_metadata "$metadata_file" sample_metadata || {
			log_warn "Metadata load failed - running pooled alignment"
			star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
			return $?
		}
	else
		log_warn "No metadata - running pooled alignment"
		star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi

	# Group samples by tissue
	declare -A tissue_samples
	for SRR in "${rnaseq_list[@]}"; do
		local tissue="${sample_metadata[${SRR}_condition]:-unknown}"
		tissue_samples["$tissue"]+="$SRR "
	done

	local tissue_count=${#tissue_samples[@]}
	log_info "[STAR TISSUE] Found $tissue_count tissue types"

	if [[ $tissue_count -lt 2 ]]; then
		log_warn "Single tissue detected - using pooled alignment"
		star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi

	log_info "[STAR TISSUE] Running tissue-specific alignments for better isoform detection"

	# Run STAR per tissue
	for tissue in "${!tissue_samples[@]}"; do
		IFS=' ' read -ra tissue_srrs <<< "${tissue_samples[$tissue]}"
		log_step "STAR alignment for tissue: $tissue (${#tissue_srrs[@]} samples)"

		star_alignment_pipeline \
			--FASTA "$fasta" \
			--RNASEQ_LIST "${tissue_srrs[@]}" \
			--TISSUE_TAG "$tissue"
	done
}
