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
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# STAR CONFIGURATION - IMPORTANT PARAMETERS (tweak here)
# ==============================================================================

# Genome loading mode: NoSharedMemory (default/safe), LoadAndKeep (multi-run HPC)
STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-$(get_star_genome_load 2>/dev/null || echo NoSharedMemory)}"

# sjdbOverhang = read_length - 1; set to actual read length for best splice detection
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"

# Delete intermediate STAR files after alignment to save disk space (true/false)
# Set to "false" to keep 2-pass genome, unsorted BAM, etc. for debugging
STAR_DELETE_TRANSIENT="${STAR_DELETE_TRANSIENT:-true}"

# REPRODUCIBILITY NOTE: Salmon's EM algorithm convergence is thread-schedule dependent.
# For exact reproducibility, always use the same --threads value across runs.
# Salmon does not expose a --seed flag for its internal EM algorithm.

# GTF annotation file for splice junction detection (resolved at runtime in pipeline)
# NOTE: STAR_GTF_FILE is resolved inside star_alignment_pipeline() from $gtf_file

# Transcriptome FASTA for Salmon quantification step (auto-detected if not set)
# When using a genome FASTA for STAR, point this to the transcript-level FASTA.
# e.g. STAR_TRANSCRIPTOME_FASTA="inputs/fasta/reference_genomes/GPE001970_transcripts.fa"
# Auto-detect: looks for <genome_basename>_transcripts.fa; falls back to --FASTA if absent.

# STAR temp directory: "system"=/tmp, "local"=output dir, "cwd"=current dir, "none"=STAR default
STAR_TMP_MODE="${STAR_TMP_MODE:-cwd}"

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
		_seq=$(zcat -f "$trimmed1" 2>/dev/null | sed -n '2p')
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
	local -a sample_names=() unique_rates=() multi_rates=() unmapped_short=() unmapped_mismatch=() unmapped_other=()

	for SRR in "${srr_list[@]}"; do
		local logf="${align_dir}/${SRR}_Log.final.out"
		[[ ! -f "$logf" ]] && continue

		# Extract key metrics from STAR Log.final.out
		local uniq_pct multi_pct short_pct mismatch_pct other_pct input_reads
		input_reads=$(awk -F'|' '/Number of input reads/{gsub(/[[:space:]]/, "", $2); print $2}' "$logf" 2>/dev/null)
		uniq_pct=$(awk -F'|' '/Uniquely mapped reads %/{gsub(/[[:space:]%]/, "", $2); print $2}' "$logf" 2>/dev/null)
		multi_pct=$(awk -F'|' '/% of reads mapped to multiple loci/{gsub(/[[:space:]%]/, "", $2); print $2}' "$logf" 2>/dev/null)
		short_pct=$(awk -F'|' '/% of reads unmapped: too short/{gsub(/[[:space:]%]/, "", $2); print $2}' "$logf" 2>/dev/null)
		mismatch_pct=$(awk -F'|' '/% of reads unmapped: too many mismatches/{gsub(/[[:space:]%]/, "", $2); print $2}' "$logf" 2>/dev/null)
		other_pct=$(awk -F'|' '/% of reads unmapped: other/{gsub(/[[:space:]%]/, "", $2); print $2}' "$logf" 2>/dev/null)

		[[ -z "$uniq_pct" ]] && continue

		sample_names+=("$SRR")
		unique_rates+=("$uniq_pct")
		multi_rates+=("${multi_pct:-0}")
		unmapped_short+=("${short_pct:-0}")
		unmapped_mismatch+=("${mismatch_pct:-0}")
		unmapped_other+=("${other_pct:-0}")

		# Per-sample checks
		# 1. Low unique mapping rate (<50% is very concerning, <70% is a warning)
		if awk "BEGIN{exit !($uniq_pct < 50)}" 2>/dev/null; then
			log_warn "[STAR QC] $SRR: Uniquely mapped only ${uniq_pct}% — VERY LOW (check sample quality, adapter contamination, or genome mismatch)"
			((warn_count++))
		elif awk "BEGIN{exit !($uniq_pct < 70)}" 2>/dev/null; then
			log_warn "[STAR QC] $SRR: Uniquely mapped ${uniq_pct}% — below 70% threshold"
			((warn_count++))
		fi

		# 2. High multi-mapping (>20%) may indicate rRNA contamination or repetitive sequences
		if awk "BEGIN{exit !(${multi_pct:-0} > 20)}" 2>/dev/null; then
			log_warn "[STAR QC] $SRR: Multi-mapped ${multi_pct}% — high rate may indicate rRNA contamination or repetitive element enrichment"
			((warn_count++))
		fi

		# 3. High unmapped-too-short (>15%) suggests adapter contamination or degraded RNA
		if awk "BEGIN{exit !(${short_pct:-0} > 15)}" 2>/dev/null; then
			log_warn "[STAR QC] $SRR: Unmapped (too short) ${short_pct}% — check adapter trimming or RNA degradation"
			((warn_count++))
		fi

		# 4. High unmapped-too-many-mismatches (>5%) suggests genome version mismatch
		if awk "BEGIN{exit !(${mismatch_pct:-0} > 5)}" 2>/dev/null; then
			log_warn "[STAR QC] $SRR: Unmapped (mismatches) ${mismatch_pct}% — may indicate genome/species mismatch"
			((warn_count++))
		fi
	done

	# Cohort-level outlier detection: flag samples >2 SD below mean unique mapping rate
	local n=${#unique_rates[@]}
	if [[ $n -ge 3 ]]; then
		local sum=0 sum_sq=0
		for rate in "${unique_rates[@]}"; do
			sum=$(awk "BEGIN{printf \"%.4f\", $sum + $rate}")
			sum_sq=$(awk "BEGIN{printf \"%.4f\", $sum_sq + ($rate * $rate)}")
		done
		local mean=$(awk "BEGIN{printf \"%.2f\", $sum / $n}")
		local sd=$(awk "BEGIN{v=($sum_sq/$n) - ($sum/$n)^2; printf \"%.2f\", (v>0)?sqrt(v):0}")
		local threshold=$(awk "BEGIN{printf \"%.2f\", $mean - 2 * $sd}")

		log_info "[STAR QC] Cohort alignment stats: mean=${mean}%, SD=${sd}%, outlier threshold=${threshold}%"

		for ((i=0; i<n; i++)); do
			if awk "BEGIN{exit !(${unique_rates[$i]} < $threshold)}" 2>/dev/null; then
				log_warn "[STAR QC] OUTLIER: ${sample_names[$i]} (${unique_rates[$i]}%) is >2 SD below cohort mean (${mean}%)"
				((warn_count++))
			fi
		done
	fi

	# Summary table
	if [[ ${#sample_names[@]} -gt 0 ]]; then
		log_info "[STAR QC] ┌───────────────────┬────────┬────────┬──────────┬──────────┬────────┐"
		log_info "[STAR QC] │ Sample            │ Unique │ Multi  │ Unmap:Sh │ Unmap:MM │ Unmap:O│"
		log_info "[STAR QC] ├───────────────────┼────────┼────────┼──────────┼──────────┼────────┤"
		for ((i=0; i<${#sample_names[@]}; i++)); do
			printf -v _row "[STAR QC] │ %-17s │ %5s%% │ %5s%% │   %5s%% │   %5s%% │ %5s%%│" \
				"${sample_names[$i]}" "${unique_rates[$i]}" "${multi_rates[$i]}" \
				"${unmapped_short[$i]}" "${unmapped_mismatch[$i]}" "${unmapped_other[$i]}"
			log_info "$_row"
		done
		log_info "[STAR QC] └───────────────────┴────────┴────────┴──────────┴──────────┴────────┘"
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
	STAR_GTF_FILE="${STAR_GTF_FILE:-$gtf_file}"

	local fasta="" transcriptome_fasta="" rnaseq_list=() tissue_tag=""

	# Check for GNU parallel
	if ! command -v parallel >/dev/null 2>&1; then
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
		fasta_dir="$(dirname "$fasta")"
		fasta_stem="$(basename "${fasta%.*}")"
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
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	set_fasta_output_dirs "$fasta_tag"

	# Get absolute paths - using realpath for robustness, fallback to manual resolution
	local abs_star_index_root abs_star_align_root

	# First ensure the directories exist
	mkdir -p "$STAR_INDEX_ROOT" "$STAR_ALIGN_ROOT" 2>/dev/null || true

	# Get absolute paths using realpath if available, otherwise use cd/pwd
	if command -v realpath >/dev/null 2>&1; then
		abs_star_index_root="$(realpath -m "$STAR_INDEX_ROOT" 2>/dev/null)" || abs_star_index_root=""
		abs_star_align_root="$(realpath -m "$STAR_ALIGN_ROOT" 2>/dev/null)" || abs_star_align_root=""
	fi

	# Fallback: use cd/pwd method
	if [[ -z "$abs_star_index_root" ]]; then
		abs_star_index_root="$(cd "$STAR_INDEX_ROOT" 2>/dev/null && pwd)" || abs_star_index_root="$STAR_INDEX_ROOT"
	fi
	if [[ -z "$abs_star_align_root" ]]; then
		abs_star_align_root="$(cd "$STAR_ALIGN_ROOT" 2>/dev/null && pwd)" || abs_star_align_root="$STAR_ALIGN_ROOT"
	fi

	# If still relative, prepend PROJECT_ROOT
	[[ "$abs_star_index_root" != /* ]] && abs_star_index_root="${PROJECT_ROOT}/${abs_star_index_root}"
	[[ "$abs_star_align_root" != /* ]] && abs_star_align_root="${PROJECT_ROOT}/${abs_star_align_root}"

	# Clean up any double slashes
	abs_star_index_root="${abs_star_index_root//\/\//\/}"
	abs_star_align_root="${abs_star_align_root//\/\//\/}"

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
			log_info "[STAR] Overriding to detected length. Export STAR_READ_LENGTH explicitly to keep configured value."
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
		genome_bp=$(awk '!/^>/{total+=length($0)} END{printf "%d", total+0}' "$fasta" 2>/dev/null)
		num_seqs=$(grep -c '^>' "$fasta" 2>/dev/null || echo "1")
		if [[ -n "$genome_bp" && "$genome_bp" -gt 0 ]]; then
			genome_sa_index=$(awk "BEGIN{v=int(log($genome_bp)/log(2)/2-1); print (v<14)?v:14}")
			[[ "$genome_sa_index" -lt 1 ]] && genome_sa_index=1
			# genomeChrBinNbits: min(18, floor(log2(GenomeLength/NumberOfSequences) - 1))
			# Required for genomes with many scaffolds (e.g., draft assemblies)
			genome_chr_bin=$(awk "BEGIN{v=int(log($genome_bp/$num_seqs)/log(2)-1); print (v<18)?v:18}")
			[[ "$genome_chr_bin" -lt 1 ]] && genome_chr_bin=1
		fi
		log_info "[STAR INDEX] Genome size: ${genome_bp:-unknown} bp, $num_seqs sequences -> genomeSAindexNbases=$genome_sa_index, genomeChrBinNbits=$genome_chr_bin"

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
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
		log_step "[PARALLEL] STAR alignment: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		log_warn "[PARALLEL] STAR is memory-intensive (~30GB/instance). Ensure sufficient RAM for $parallel_jobs concurrent jobs."
		_prepare_parallel_env

		export fasta_tag threads_per_job
		export star_index_dir star_genome_dir
		export STAR_DELETE_TRANSIENT PROJECT_ROOT
		export effective_genome_load
		# Serialize array args for parallel workers (bash can't export arrays)
		export _star_strand_args_str="${star_strand_args[*]}"
		export _star_sj_filter_args_str="${star_sj_filter_args[*]}"

		_m3_star_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log STAR "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"

			# Check if BAM exists and has content
			if [[ -f "$bam_output" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				local bam_size
				bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
				if [[ "$bam_size" -gt 1000 ]]; then
					_parallel_log STAR "$SRR" INFO "BAM exists (${bam_size} bytes) - skipping"
					return 0
				fi
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

			local star_tmp_dir="${star_genome_dir}/_STARtmp_${SRR}"
			rm -rf "$star_tmp_dir" 2>/dev/null || true

			local out_prefix="${star_genome_dir}/${SRR}_"
			out_prefix="${out_prefix//\/\//\/}"
			local unsorted_bam="${out_prefix}Aligned.out.bam"

			_parallel_log STAR "$SRR" INFO "Aligning with $threads_per_job threads"
			_parallel_log STAR "$SRR" INFO "--- BEGIN STAR OUTPUT ---"

			# Deserialize array args from exported strings
			local _par_strand_args=($_star_strand_args_str)
			local _par_sj_args=($_star_sj_filter_args_str)

			STAR --runMode alignReads \
				--genomeDir "$star_index_dir" \
				--readFilesIn "${star_reads_args[@]}" \
				--readFilesCommand "zcat -f" \
				--outFileNamePrefix "$out_prefix" \
				--outTmpDir "$star_tmp_dir" \
				--outSAMtype BAM Unsorted \
				${_par_strand_args[@]:+"${_par_strand_args[@]}"} \
				--outSAMattributes NH HI AS NM MD \
				--outSAMunmapped Within \
				--twopassMode Basic \
				${_par_sj_args[@]:+"${_par_sj_args[@]}"} \
				--genomeLoad "$effective_genome_load" \
				--runThreadN "$threads_per_job" 2>&1 || \
				{ _parallel_log STAR "$SRR" ERROR "--- END STAR OUTPUT (FAILED) ---"; _parallel_log STAR "$SRR" ERROR "STAR alignment failed"; return 1; }

			_parallel_log STAR "$SRR" INFO "--- END STAR OUTPUT ---"

			if [[ ! -f "$unsorted_bam" ]]; then
				_parallel_log STAR "$SRR" ERROR "Unsorted BAM not created"
				return 1
			fi

			local unsorted_size
			unsorted_size=$(stat -c%s "$unsorted_bam" 2>/dev/null || stat -f%z "$unsorted_bam" 2>/dev/null || echo "0")
			if [[ "$unsorted_size" -lt 1000 ]]; then
				_parallel_log STAR "$SRR" ERROR "Unsorted BAM is empty/corrupt (${unsorted_size} bytes)"
				return 1
			fi

			_parallel_log STAR "$SRR" INFO "Sorting BAM with samtools"
			samtools sort -@ "$threads_per_job" -m 2G -o "$bam_output" "$unsorted_bam" 2>&1 | sed 's/^/\t/'
			local sort_exit=${PIPESTATUS[0]}
			[[ $sort_exit -ne 0 ]] && { _parallel_log STAR "$SRR" ERROR "samtools sort failed"; return 1; }

			local final_bam_size
			final_bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
			if [[ "$final_bam_size" -lt 1000 ]]; then
				_parallel_log STAR "$SRR" ERROR "Sorted BAM is empty/corrupt (${final_bam_size} bytes)"
				return 1
			fi

			rm -f "$unsorted_bam"
			samtools index -@ "$threads_per_job" "$bam_output" 2>&1 || true

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
		export -f _m3_star_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env threads_per_job \
			--env star_index_dir --env star_genome_dir \
			--env STAR_DELETE_TRANSIENT --env PROJECT_ROOT \
			--env OVERWRITE_MODE --env effective_genome_load \
			--env _star_strand_args_str --env _star_sj_filter_args_str \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$star_genome_dir/parallel_star_align.log" \
			_m3_star_parallel_worker {}

		local par_exit=$?
		log_info "[PARALLEL] STAR alignment complete (exit=$par_exit)"
		[[ $par_exit -ne 0 ]] && log_warn "[PARALLEL] Some jobs failed - check $star_genome_dir/parallel_star_align.log"
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"

			# Check if BAM exists AND has content (not 0 bytes from failed run)
			if [[ -f "$bam_output" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				local bam_size
				bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
				if [[ "$bam_size" -gt 1000 ]]; then
					log_info "[STAR] Alignment for $SRR already exists (${bam_size} bytes). Skipping."
					continue
				else
					log_warn "[STAR] Found empty/corrupt BAM for $SRR (${bam_size} bytes) - removing and re-running"
					rm -f "$bam_output"
				fi
			fi

			# Clean up stale files from previous failed STAR runs
			log_info "[STAR] Cleaning up stale files for $SRR..."
			rm -f "${star_genome_dir}/${SRR}_Log.out" "${star_genome_dir}/${SRR}_Log.progress.out" \
				"${star_genome_dir}/${SRR}_Log.final.out" "${star_genome_dir}/${SRR}_SJ.out.tab" 2>/dev/null || true
			rm -rf "${star_genome_dir}/${SRR}__STARgenome" "${star_genome_dir}/${SRR}__STARpass1" \
				"${star_genome_dir}/${SRR}_STARtmp" "${star_genome_dir}/${SRR}_"*.tmp \
				"${star_genome_dir}/_STARtmp_${SRR}" "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true

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

			# CRITICAL: Ensure output directory exists and is writable RIGHT BEFORE STAR runs
			# This fixes "could not create output file" errors on HPC/server environments
			mkdir -p "$star_genome_dir"
			# On NFS/HPC systems, sync and a brief pause help avoid race conditions.
			# Set STAR_NFS_SYNC=true in your config if running on NFS storage.
			if [[ "${STAR_NFS_SYNC:-false}" == "true" ]]; then
				sync
				sleep 1
			fi

			if [[ ! -d "$star_genome_dir" ]]; then
				log_error "[STAR] FATAL: Cannot create output directory: $star_genome_dir"
				return 1
			fi

			# Test write permissions by creating a test file with the exact output name pattern
			local test_file="${star_genome_dir}/${SRR}_test_write_$$"
			if ! touch "$test_file" 2>/dev/null; then
				log_error "[STAR] FATAL: Cannot create output file: $test_file"
				log_error "[STAR] Check disk space and permissions for: $star_genome_dir"
				return 1
			fi
			rm -f "$test_file"

			# Configure STAR temp directory - ALWAYS use output directory for temp
			# This ensures all STAR operations happen on the same filesystem
			local star_tmp_dir="${star_genome_dir}/_STARtmp_${SRR}"
			rm -rf "$star_tmp_dir" 2>/dev/null || true
			# Also clean up any stray temp dirs
			rm -rf "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true
			rm -rf "${PROJECT_ROOT}/_STARtmp_${SRR}" 2>/dev/null || true
			log_info "[STAR] Using temp directory: $star_tmp_dir"

			# Construct output prefix - ensure no double slashes and path is clean
			local out_prefix="${star_genome_dir}/${SRR}_"
			# Remove any double slashes
			out_prefix="${out_prefix//\/\//\/}"

			# Check disk space before running STAR (needs ~30GB per run)
			local available_space
			available_space=$(df -P "$star_genome_dir" 2>/dev/null | awk 'NR==2 {print $4}')
			if [[ -n "$available_space" ]]; then
				local available_gb=$((available_space / 1024 / 1024))
				if [[ $available_gb -lt 30 ]]; then
					log_error "[STAR] FATAL: Less than 30GB available in $star_genome_dir - STAR needs ~30GB per run"
					return 1
				fi
			fi

			# Run STAR alignment - output UNSORTED BAM first, then sort with samtools
			local unsorted_bam="${out_prefix}Aligned.out.bam"

			run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$star_genome_dir" \
				STAR --runMode alignReads \
					--genomeDir "$star_index_dir" \
					--readFilesIn "${star_reads_args[@]}" \
					--readFilesCommand "zcat -f" \
					--outFileNamePrefix "$out_prefix" \
					--outTmpDir "$star_tmp_dir" \
					--outSAMtype BAM Unsorted \
					${star_strand_args[@]:+"${star_strand_args[@]}"} \
					--outSAMattributes NH HI AS NM MD \
					--outSAMunmapped Within \
					--twopassMode Basic \
					${star_sj_filter_args[@]:+"${star_sj_filter_args[@]}"} \
					--genomeLoad "$effective_genome_load" \
					--runThreadN "$THREADS"

			# Check if unsorted BAM was created
			if [[ ! -f "$unsorted_bam" ]]; then
				log_error "[STAR] FATAL: Unsorted BAM file not created for $SRR"
				log_error "[STAR] Check STAR log: ${out_prefix}Log.out"
				# Show last 30 lines of STAR log
				log_error "[STAR] Last 30 lines of STAR log:"
				tail -30 "${out_prefix}Log.out" 2>/dev/null | while IFS= read -r line; do log_error "  $line"; done
				return 1
			fi

			local unsorted_size
			unsorted_size=$(stat -c%s "$unsorted_bam" 2>/dev/null || stat -f%z "$unsorted_bam" 2>/dev/null || echo "0")
			if [[ "$unsorted_size" -lt 1000 ]]; then
				log_error "[STAR] FATAL: Unsorted BAM is empty/corrupt for $SRR (${unsorted_size} bytes)"
				log_error "[STAR] Last 30 lines of STAR log:"
				tail -30 "${out_prefix}Log.out" 2>/dev/null | while IFS= read -r line; do log_error "  $line"; done
				return 1
			fi

			# Sort BAM with samtools
			log_info "[STAR] Sorting BAM with samtools..."

			if ! samtools sort -@ "$THREADS" -m 2G -o "$bam_output" "$unsorted_bam" 2>&1; then
				log_error "[STAR] FATAL: samtools sort failed for $SRR"
				return 1
			fi

			# Verify sorted BAM
			local final_bam_size
			final_bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
			if [[ "$final_bam_size" -lt 1000 ]]; then
				log_error "[STAR] FATAL: Sorted BAM file is empty/corrupt for $SRR (${final_bam_size} bytes)"
				return 1
			fi

			log_info "[STAR] BAM sorted successfully: $final_bam_size bytes"

			# Remove unsorted BAM to save space
			rm -f "$unsorted_bam"
			log_info "[STAR] Removed unsorted BAM to save space"

			# Index the BAM
			log_info "[STAR] Indexing BAM..."
			samtools index -@ "$THREADS" "$bam_output" 2>&1 || log_warn "[STAR] BAM indexing failed (non-fatal)"

			# Clean up temp directories after successful alignment
			rm -rf "$star_tmp_dir" "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true

			# Delete transient big files if enabled (saves significant disk space)
			if [[ "${STAR_DELETE_TRANSIENT:-true}" == "true" ]]; then
				log_info "[STAR] Cleaning up transient files for $SRR..."

				# Remove 2-pass intermediate directories (can be several GB each)
				rm -rf "${star_genome_dir}/${SRR}__STARgenome" 2>/dev/null || true
				rm -rf "${star_genome_dir}/${SRR}__STARpass1" 2>/dev/null || true

				# Remove any leftover temp directories
				rm -rf "${star_genome_dir}/${SRR}_STARtmp" 2>/dev/null || true
				rm -rf "${star_genome_dir}/_STARtmp_${SRR}" 2>/dev/null || true

				# Remove progress log (Log.out and Log.final.out are kept for diagnostics)
				rm -f "${star_genome_dir}/${SRR}_Log.progress.out" 2>/dev/null || true

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

	# Quantify samples
	log_info "[SALMON] Starting quantification for ${#rnaseq_list[@]} samples"

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
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
			# Strip ANSI escape codes and carriage returns (Salmon uses colored progress bars)
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				salmon quant -p "$threads_per_job" -i "$salmon_idx" -o "$quant_dir" \
					--validateMappings --gcBias -l "${_sal_lib_pe:-A}" -1 "$trimmed1" -2 "$trimmed2" 2>&1 | \
					sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				quant_exit=${PIPESTATUS[0]}
			else
				salmon quant -p "$threads_per_job" -i "$salmon_idx" -o "$quant_dir" \
					--validateMappings --gcBias -l "${_sal_lib_se:-A}" -r "$trimmed1" 2>&1 | \
					sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				quant_exit=${PIPESTATUS[0]}
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
					--validateMappings --gcBias -l "${_sal_lib_pe:-A}" -1 "$trimmed1" -2 "$trimmed2"
			else
				run_with_space_time_log salmon quant -p "$THREADS" -i "$salmon_idx" -o "$quant_dir" \
					--validateMappings --gcBias -l "${_sal_lib_se:-A}" -r "$trimmed1"
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
	if [[ ! -f "$tx2gene_file" ]]; then
		log_info "[TXIMPORT] Creating transcript-to-gene mapping from GTF: $STAR_GTF_FILE"
		awk '$3=="transcript" {
			tid=""; gid=""
			for (i=9; i<=NF; i++) {
				if ($i == "transcript_id") { gsub(/[";]/, "", $(i+1)); tid=$(i+1) }
				if ($i == "gene_id")       { gsub(/[";]/, "", $(i+1)); gid=$(i+1) }
			}
			if (tid != "" && gid != "") print tid "\t" gid
		}' "$STAR_GTF_FILE" | sort -u > "$tx2gene_file"
		local tx2gene_count
		tx2gene_count=$(wc -l < "$tx2gene_file")
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
	local helper_script="$SCRIPT_DIR/../c_post_processing/preprocessing/STAR/tximport_star_helper.R"

	if [[ ! -f "$helper_script" ]]; then
		log_error "tximport_star_helper.R not found: $helper_script"
		return 1
	fi

	log_info "[TXIMPORT] Running STAR+Salmon import..."
	Rscript "$helper_script" "$quant_dir" "$metadata_file" "$tx2gene_file" "$output_dir" "$master_ref"
}

# Generate tximport script (legacy - copies helper)
generate_tximport_star_script() {
	local quant_root="$1"
	local sample_metadata="$2"
	local tx2gene_file="$3"
	local matrix_dir="$4"
	local output_script="$5"
	local helper_script="$SCRIPT_DIR/../c_post_processing/preprocessing/STAR/tximport_star_helper.R"

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
		local tissue_srrs=(${tissue_samples[$tissue]})
		log_step "STAR alignment for tissue: $tissue (${#tissue_srrs[@]} samples)"

		star_alignment_pipeline \
			--FASTA "$fasta" \
			--RNASEQ_LIST "${tissue_srrs[@]}" \
			--TISSUE_TAG "$tissue"
	done
}
