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
M1_HISAT2_REF_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/b_main_methods}"
if [[ -z "$SCRIPT_DIR" ]]; then
	SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
	SCRIPT_DIR="$(cd "$SCRIPT_DIR" 2>/dev/null && pwd)"
fi
source "$SCRIPT_DIR/shared_utils_method.sh"

# Binary availability cached in shared_utils_method.sh: _SHARED_HAS_SAMTOOLS, _SHARED_HAS_PARALLEL

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
		# Single awk pass extracts all 3 metrics (replaces 5-7 grep|grep pipelines)
		local overall conc1 concm
		# Reset before eval — bash 'local' doesn't re-initialize already-local vars
		overall="" conc1="" concm=""
		# Uses POSIX-portable match()+substr() instead of gawk-only capture groups
		eval "$(awk '
			/overall alignment rate/ { gsub(/%.*/, ""); printf "overall=%s ", $NF }
			/aligned concordantly exactly 1 time/ { if (match($0, /[0-9.]+%/)) { v=substr($0, RSTART, RLENGTH-1); printf "conc1=%s ", v; conc1_set=1 } }
			/aligned concordantly >1 time/        { if (match($0, /[0-9.]+%/)) { v=substr($0, RSTART, RLENGTH-1); printf "concm=%s ", v; concm_set=1 } }
			# Single-end fallback patterns (guard prevents overwrite by non-concordant mate lines in PE data)
			/aligned exactly 1 time/ && !/concordantly/ { if (!conc1_set && match($0, /[0-9.]+%/)) { v=substr($0, RSTART, RLENGTH-1); printf "conc1=%s ", v; conc1_set=1 } }
			/aligned >1 time/ && !/concordantly/        { if (!concm_set && match($0, /[0-9.]+%/)) { v=substr($0, RSTART, RLENGTH-1); printf "concm=%s ", v; concm_set=1 } }
		' "$sumf" 2>/dev/null)"

		[[ -z "$overall" ]] && continue

		sample_names+=("$SRR")
		overall_rates+=("$overall")
		concordant_once+=("${conc1:-0}")
		concordant_multi+=("${concm:-0}")

		# Per-sample checks (truncate to integer for fast bash arithmetic — no awk spawns)
		local _ov_int=${overall%.*}
		if (( _ov_int < 50 )); then
			log_warn "[HISAT2 QC] $SRR: Overall alignment ${overall}% — VERY LOW (check sample quality, adapter contamination, or genome mismatch)"
			((warn_count++)) || true
		elif (( _ov_int < 70 )); then
			log_warn "[HISAT2 QC] $SRR: Overall alignment ${overall}% — below 70% threshold"
			((warn_count++)) || true
		fi

		# High multi-mapping (>20%) may indicate rRNA contamination
		local _cm_int=${concm:-0}; _cm_int=${_cm_int%.*}
		if (( _cm_int > 20 )); then
			log_warn "[HISAT2 QC] $SRR: Multi-mapped ${concm}% — may indicate rRNA contamination or repetitive element enrichment"
			((warn_count++)) || true
		fi
	done

	# Cohort-level outlier detection: flag samples >2 SD below mean overall rate
	local n=${#overall_rates[@]}
	if [[ $n -ge 3 ]]; then
		# Single AWK pass: compute stats AND find outliers (1 process instead of 2)
		local _stats_and_outliers
		# Herestring avoids printf subprocess fork — awk splits on RS=' '
		_stats_and_outliers=$(awk -v RS=' ' 'NF{
			vals[n]=$1+0; s+=$1; ss+=$1*$1; n++
		} END {
			m=s/n; v=ss/n - m*m; sd=(v>0)?sqrt(v):0; thr=m-2*sd
			printf "%.2f %.2f %.2f", m, sd, thr
			for (i=0; i<n; i++) if (vals[i] < thr) printf " %d", i
		}' <<< "${overall_rates[*]}")
		local mean sd threshold _outlier_tail
		read -r mean sd threshold _outlier_tail <<< "$_stats_and_outliers"

		log_info "[HISAT2 QC] Cohort alignment stats: mean=${mean}%, SD=${sd}%, outlier threshold=${threshold}%"

		# Extract outlier indices (everything after the 3 stats fields)
		local _outlier_indices="${_stats_and_outliers#* * * }"
		[[ "$_outlier_indices" == "$_stats_and_outliers" ]] && _outlier_indices=""
		local _oi
		for _oi in $_outlier_indices; do
			log_warn "[HISAT2 QC] OUTLIER: ${sample_names[$_oi]} (${overall_rates[$_oi]}%) is >2 SD below cohort mean (${mean}%)"
			((warn_count++)) || true
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

	# Cache validation result: skip re-scanning multi-GB files on resume runs.
	# Sentinel is invalidated when either input file is newer.
	local _sentinel="${HISAT2_REF_GUIDED_INDEX_DIR:-.}/.fasta_gtf_validated"
	if [[ -f "$_sentinel" && "$fasta" -ot "$_sentinel" && "$gtf" -ot "$_sentinel" ]]; then
		return 0
	fi

	# Use .fai index if available (a few KB vs scanning multi-GB FASTA for >headers)
	local fasta_source="$fasta"
	local use_fai=false
	if [[ -f "${fasta}.fai" ]]; then
		fasta_source="${fasta}.fai"
		use_fai=true
	fi

	# Single AWK pass: extract chromosome names, compute overlap
	local result
	if [[ "$use_fai" == "true" ]]; then
		result=$(awk '
			# Pass 1: .fai file (col1 = sequence name, first 20)
			FILENAME == ARGV[1] && fasta_n < 20 {
				if (!($1 in fasta_seen)) { fasta_seen[$1]=1; fasta_n++ }
			}
			# Pass 2: GTF chromosome names (first 20 unique, skip comments)
			FILENAME == ARGV[2] && !/^#/ && gtf_n < 20 {
				if (!($1 in gtf_seen)) { gtf_seen[$1]=1; gtf_n++ }
			}
			END {
				overlap = 0
				for (c in gtf_seen) if (c in fasta_seen) overlap++
				printf "%d %d %d ", overlap, fasta_n, gtf_n
				n=0; for (c in fasta_seen) { if (n<5) printf "%s ", c; n++ }
				printf "| "
				n=0; for (c in gtf_seen) { if (n<5) printf "%s ", c; n++ }
			}
		' "$fasta_source" "$gtf" 2>/dev/null)
	else
		# Two-step approach: (1) awk on FASTA with early exit after 20 unique headers
		# (avoids full scan of multi-GB files), (2) single awk on GTF that also computes
		# overlap using FASTA chroms passed via variable. 2 processes instead of 3.
		# O(header_positions) for FASTA; O(GTF_lines_until_20_unique) for GTF.
		local _fasta_chroms
		_fasta_chroms=$(awk '/^>/ { sub(/^>/, ""); sub(/[[:space:]].*/, ""); if (!seen[$0]++) { print; if (++n >= 20) exit } }' "$fasta" 2>/dev/null)

		# Single awk: load FASTA chroms from stdin (NR==FNR), then scan GTF and compute overlap
		result=$(awk '
			NR == FNR { if ($0 != "") { fasta_seen[$0]=1; fasta_n++ }; next }
			!/^#/ && gtf_n < 20 { if (!($1 in gtf_seen)) { gtf_seen[$1]=1; gtf_n++ } }
			END {
				overlap = 0
				for (c in gtf_seen) if (c in fasta_seen) overlap++
				printf "%d %d %d ", overlap, fasta_n, gtf_n
				n=0; for (c in fasta_seen) { if (n<5) printf "%s ", c; n++ }
				printf "| "
				n=0; for (c in gtf_seen) { if (n<5) printf "%s ", c; n++ }
			}
		' <(printf '%s\n' "$_fasta_chroms") "$gtf" 2>/dev/null)
	fi

	if [[ -z "$result" ]]; then
		log_warn "[VALIDATE] Could not extract chromosome names from FASTA or GTF"
		return 0
	fi

	local overlap fasta_count gtf_count _fasta_chrs
	read -r overlap fasta_count gtf_count _fasta_chrs <<< "${result%% |*}"
	local gtf_first5="${result#* | }"
	local fasta_first5="$_fasta_chrs"

	if [[ "$overlap" -eq 0 ]]; then
		log_warn "[VALIDATE] FASTA/GTF MISMATCH: No shared chromosome names between FASTA and GTF!"
		log_warn "[VALIDATE]   FASTA chromosomes (first 5): $fasta_first5"
		log_warn "[VALIDATE]   GTF chromosomes (first 5):   $gtf_first5"
		log_warn "[VALIDATE]   Splice site guidance will be INEFFECTIVE — alignment quality degraded"
	elif [[ "$overlap" -lt "$gtf_count" ]]; then
		log_info "[VALIDATE] FASTA/GTF partial overlap: $overlap of $gtf_count GTF chromosomes found in FASTA"
	else
		log_info "[VALIDATE] FASTA/GTF chromosome names match ($overlap shared)"
	fi

	# Write sentinel so subsequent runs skip this validation
	touch "$_sentinel" 2>/dev/null || true
}

# ==============================================================================
# BAM METRICS: BASE TRIMMING & INSERT SIZE
# ==============================================================================
# Collect lightweight BAM metrics using samtools stats. Must run BEFORE BAM deletion.
# Saves a summary TSV per sample for post-hoc review.
# Usage: _m1_collect_bam_metrics <bam> <output_dir> [method_tag] [srr_tag]

_m1_collect_bam_metrics() {
	local bam="$1" out_dir="$2"
	local method="${3:-HISAT2_RG}" srr="${4:-SAMPLE}"
	local metrics_file="$out_dir/${srr}_bam_metrics.txt"

	# Use cached samtools availability (set in shared_utils_method.sh) to avoid per-BAM command -v
	if ! $_SHARED_HAS_SAMTOOLS || [[ ! -f "$bam" ]]; then
		return 0
	fi

	# Skip if metrics already collected and BAM hasn't changed (resume optimization)
	# O(1) mtime comparison avoids O(N) samtools stats scan on unchanged BAMs
	if [[ -f "$metrics_file" && "$metrics_file" -nt "$bam" ]]; then
		return 0
	fi

	# samtools stats with multi-threading (cap at 4 — I/O bound beyond that)
	local _st_threads=${threads_per_job:-${THREADS:-4}}
	(( _st_threads > 4 )) && _st_threads=4
	local stats
	# O(N) — samtools stats performs a linear scan of all N alignments in the BAM.
	# Single awk pass directly on samtools output — eliminates redundant grep subprocess
	# since awk already filters on /^SN\t/ patterns.
	local total_bases bases_clipped insert_mean insert_sd
	eval "$(samtools stats -@ "$_st_threads" "$bam" 2>/dev/null | awk -F'\t' '
		/^SN\tbases mapped \(cigar\)/               { printf "total_bases=%s ", $3 }
		/^SN\tbases trimmed/                         { printf "bases_clipped=%s ", $3 }
		/^SN\tinsert size average/                   { printf "insert_mean=%s ", $3 }
		/^SN\tinsert size standard deviation/        { printf "insert_sd=%s ", $3 }
	')" || return 0

	# Save metrics to file for post-hoc review
	{
		echo "sample=$srr"
		echo "total_mapped_bases=$total_bases"
		echo "bases_trimmed=$bases_clipped"
		echo "insert_size_mean=$insert_mean"
		echo "insert_size_sd=$insert_sd"
	} > "$metrics_file"

	# Flag excessive base trimming (>10% of mapped bases)
	# O(1) bash integer math: 1000*clipped/total vs threshold 100 (=10.0%)
	# Eliminates awk subprocess for float comparison
	if [[ -n "$total_bases" && -n "$bases_clipped" && "$total_bases" -gt 0 ]]; then
		local clip_pct_x10=$(( (1000 * bases_clipped) / total_bases ))
		local clip_pct="$(( clip_pct_x10 / 10 )).$(( clip_pct_x10 % 10 ))"
		if (( clip_pct_x10 > 100 )); then
			if [[ -n "${abs_error_warn_file:-}" ]]; then
				_parallel_log "$method" "$srr" WARN "Bases trimmed ${clip_pct}% of mapped bases — may indicate quality issues or adapter contamination"
			else
				log_warn "[$method] $srr: Bases trimmed ${clip_pct}% of mapped bases — may indicate quality issues or adapter contamination"
			fi
		fi
	fi

	# Flag abnormal insert size (outside 100-800 bp for typical RNA-seq)
	if [[ -n "$insert_mean" ]]; then
		local mean_int=${insert_mean%.*}
		if [[ "$mean_int" -gt 0 ]]; then
			if [[ "$mean_int" -lt 100 || "$mean_int" -gt 800 ]]; then
				if [[ -n "${abs_error_warn_file:-}" ]]; then
					_parallel_log "$method" "$srr" WARN "Unusual insert size: mean=${insert_mean}, SD=${insert_sd:-N/A} (expected 100-800 bp for RNA-seq)"
				else
					log_warn "[$method] $srr: Unusual insert size: mean=${insert_mean}, SD=${insert_sd:-N/A} (expected 100-800 bp for RNA-seq)"
				fi
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
	[[ ${#rnaseq_list[@]} -eq 0 ]] && { log_error "No RNA-seq samples provided."; return 1; }

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
	# Pure bash parameter expansion — avoids basename subshell fork
	fasta_base="${fasta##*/}"
	fasta_tag="${fasta_base%.*}"
	set_fasta_output_dirs "$fasta_tag"
	index_prefix="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_ref_guided"

	# Validate FASTA/GTF chromosome name consistency
	_m1_validate_fasta_gtf_chromosomes "$fasta" "$gtf"

	# BUILD HISAT2 REFERENCE-GUIDED INDEX
	mkdir -p "$HISAT2_REF_GUIDED_INDEX_DIR"
	# Check for existing index via specific file (avoids ls glob subprocess + edge cases)
	if [[ -f "${index_prefix}.1.ht2" ]] && [[ "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[INDEX] Ref-Guided index exists - skipping build"
	else
		log_step "Building HISAT2 Ref-Guided index: $fasta_base"

		local splice_sites="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_splice_sites.txt"
		local exons="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_exons.txt"
		local build_opts=""

		# Extract splice sites and exons (may be empty for single-exon transcriptomes)
		if ! hisat2_extract_splice_sites.py "$gtf" > "$splice_sites" 2>/dev/null; then
			log_warn "[INDEX] hisat2_extract_splice_sites.py failed for $gtf — continuing without splice sites"
			> "$splice_sites"
		fi
		if ! hisat2_extract_exons.py "$gtf" > "$exons" 2>/dev/null; then
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
				# Use samtools faidx index for O(1) lookup when available (avoids scanning entire FASTA)
				if [[ -f "${fasta}.fai" ]] || ($_SHARED_HAS_SAMTOOLS && samtools faidx "$fasta" 2>/dev/null); then
					_ss_seq_len=$(awk -v t="$_ss_chr" '$1==t{print $2;exit}' "${fasta}.fai" 2>/dev/null)
				fi
				# Fallback: scan FASTA directly (for small FASTAs or if samtools unavailable)
				if [[ -z "${_ss_seq_len:-}" ]]; then
					_ss_seq_len=$(awk -v t="$_ss_chr" \
						'/^>/{if(f){print l; f=0; exit} n=$1; sub(/^>/,"",n); if(n==t){f=1; l=0}; next}
						 f{l+=length} END{if(f) print l}' "$fasta")
				fi
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
			# O(1) fork vs O(N) bash loop — wc -l is a single C-level scan
			local _ss_n; _ss_n=$(wc -l < "$splice_sites")
			log_info "[INDEX] Extracted $_ss_n splice sites"
			build_opts="$build_opts --ss $splice_sites"
		else
			log_warn "[INDEX] No splice sites found - GTF may contain only single-exon transcripts"
		fi

		if [[ -s "$exons" ]]; then
			# O(1) fork vs O(N) bash loop — wc -l is a single C-level scan
			local _ex_n; _ex_n=$(wc -l < "$exons")
			log_info "[INDEX] Extracted $_ex_n exons"
			build_opts="$build_opts --exon $exons"
		else
			log_warn "[INDEX] No exons extracted from GTF"
		fi

		log_file_size "$fasta" "Input FASTA for HISAT2 index"
		run_with_space_time_log --input "$fasta" --output "$HISAT2_REF_GUIDED_INDEX_DIR" \
			hisat2-build -p "${THREADS}" $build_opts "$fasta" "$index_prefix" \
			|| { log_error "[INDEX] HISAT2 ref-guided index build failed for $fasta_tag"; rm -f "${index_prefix}".*.ht2; return 1; }
		log_file_size "$HISAT2_REF_GUIDED_INDEX_DIR" "HISAT2 index output"
	fi

	# ALIGNMENT AND STRINGTIE ASSEMBLY (pure quantification: single pass with -e)
	# HISAT2+StringTie optimal: 16 threads per job (scales well up to ~16)
	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	if [[ "${_JOBS_MODE:-}" == "auto" || "${_JOBS_MODE:-}" == "AUTO" ]]; then
		parallel_jobs=$(( THREADS / 16 ))
		(( parallel_jobs < 1 )) && parallel_jobs=1
	fi
	# Adaptive: don't spawn more parallel jobs than samples — wastes thread allocation
	# O(1) min check avoids reserving e.g. 32 threads/job when only 3 samples exist on 64-thread node
	local _n_samples=${#rnaseq_list[@]}
	(( parallel_jobs > _n_samples )) && parallel_jobs=$_n_samples
	(( parallel_jobs < 1 )) && parallel_jobs=1
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if $_SHARED_HAS_PARALLEL && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
		log_step "[PARALLEL] HISAT2 Ref-Guided Align+StringTie: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env "M1"

		# Pre-compute index directory outside worker (avoids dirname subshell per sample)
		local abs_index_dir="${index_prefix%/*}"
		export fasta_tag index_prefix abs_index_dir threads_per_job hisat2_strand_opts stringtie_strand_opt OVERWRITE_MODE
		local abs_hisat2_rg_root="$HISAT2_REF_GUIDED_ROOT"
		[[ "$abs_hisat2_rg_root" != /* ]] && abs_hisat2_rg_root="${PROJECT_ROOT:-$PWD}/$abs_hisat2_rg_root"
		local abs_stringtie_rg_root="$STRINGTIE_HISAT2_REF_GUIDED_ROOT"
		[[ "$abs_stringtie_rg_root" != /* ]] && abs_stringtie_rg_root="${PROJECT_ROOT:-$PWD}/$abs_stringtie_rg_root"
		local abs_gtf="$gtf"
		[[ "$abs_gtf" != /* ]] && abs_gtf="${PROJECT_ROOT:-$PWD}/$abs_gtf"
		export abs_hisat2_rg_root abs_stringtie_rg_root abs_gtf

		_m1_align_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log HISAT2_RG "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local HISAT2_DIR="$abs_hisat2_rg_root/$SRR"
			mkdir -p "$HISAT2_DIR"
			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
			local out_gtf="$abs_stringtie_rg_root/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"

			if [[ -f "$bam" && ( -f "${bam}.bai" || -f "${bam}.csi" ) && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_RG "$SRR" INFO "BAM exists - skipping alignment"
			elif [[ -f "$out_gtf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_RG "$SRR" INFO "GTF exists for $SRR - skipping alignment"
			else
				_parallel_log HISAT2_RG "$SRR" INFO "Aligning with $threads_per_job threads"
				local summary_file="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_alignment_summary.txt"
				# Pipe hisat2 directly into samtools sort — eliminates 10-50GB SAM intermediate per sample
				# Cap sort threads at 4: samtools sort is I/O-bound beyond ~4 threads,
				# and sort memory is per-thread, so capping saves RAM without losing speed.
				# Cap sort threads: min(threads_per_job, 4); at least 1
				local sort_threads=$(( threads_per_job < 4 ? threads_per_job : 4 ))
				(( sort_threads < 1 )) && sort_threads=1
				local sort_mem _sort_wi_flag=""
				_samtools_sort_mem "$sort_threads" "$parallel_jobs" sort_mem
				_samtools_has_write_index && _sort_wi_flag="--write-index"
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "$threads_per_job" --dta $hisat2_strand_opts -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" 2>"$summary_file" \
						| samtools sort -@ "$sort_threads" -m "$sort_mem" $_sort_wi_flag -o "$bam"
				else
					hisat2 -p "$threads_per_job" --dta $hisat2_strand_opts -x "$index_prefix" \
						-U "$trimmed1" 2>"$summary_file" \
						| samtools sort -@ "$sort_threads" -m "$sort_mem" $_sort_wi_flag -o "$bam"
				fi
				local _ps=("${PIPESTATUS[@]}")
				# Display alignment summary (strip CR for clean log output)
				if [[ -s "$summary_file" ]]; then
					while IFS= read -r _line || [[ -n "$_line" ]]; do
						_line="${_line//$'\r'/}"
						printf '%s\n' "$_line"
					done < "$summary_file"
				fi
				[[ ${_ps[0]} -ne 0 ]] && { _parallel_log HISAT2_RG "$SRR" ERROR "HISAT2 failed (exit=${_ps[0]})"; rm -f "$bam"; return ${_ps[0]}; }
				[[ ${_ps[1]} -ne 0 ]] && { _parallel_log HISAT2_RG "$SRR" ERROR "samtools sort failed"; rm -f "$bam"; return 1; }
				# Only run separate index if --write-index was not used
				# Cap index threads at 4 — samtools index is I/O-bound
				if [[ -z "$_sort_wi_flag" ]]; then
					local _idx_t=$threads_per_job; (( _idx_t > 4 )) && _idx_t=4
					if ! samtools index -@ "$_idx_t" "$bam"; then
						_parallel_log HISAT2_RG "$SRR" WARN "samtools index failed — BAM may need re-indexing"
					fi
				fi

				# Infer strandness once (lock-file ensures only first worker runs it)
				[[ -z "$hisat2_strand_opts" ]] && \
					_m1_infer_strandness "$bam" "$abs_gtf" "$abs_index_dir" "HISAT2_RG" "$SRR"
			fi

			# StringTie quantification (ref-guided, single pass)
			local out_dir="$abs_stringtie_rg_root/$SRR"
			mkdir -p "$out_dir"

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

			# Only delete BAM after verifying StringTie produced valid output
			# Prevents unrecoverable data loss if StringTie exits 0 but produces empty/missing GTF
			if [[ "$keep_bam_global" != "y" && -f "$bam" ]]; then
				if [[ -f "$out_gtf" && -s "$out_gtf" ]]; then
					_parallel_log HISAT2_RG "$SRR" WARN "Deleting BAM to save disk (set keep_bam_global=y to retain): ${bam##*/}"
					rm -f "$bam" "${bam}.bai" "${bam}.csi"
				else
					_parallel_log HISAT2_RG "$SRR" WARN "Retaining BAM — StringTie GTF missing or empty: ${out_gtf##*/}"
				fi
			fi

			_parallel_log HISAT2_RG "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m1_align_parallel_worker _m1_infer_strandness _m1_collect_bam_metrics
		# Pre-detect capabilities so parallel workers don't each test independently
		_samtools_has_write_index || true
		_get_available_ram_mb > /dev/null

		# O(S/parallel_jobs × (N×log N + T×G)) — S samples dispatched across parallel_jobs slots;
		# each worker runs HISAT2 O(N log N) + samtools sort O(N log N) + StringTie O(T×G)
		parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env _CONDA_PROFILE_SCRIPT --env WF_MANAGED_ENV \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env index_prefix --env threads_per_job \
			--env abs_hisat2_rg_root --env abs_stringtie_rg_root --env abs_gtf \
			--env hisat2_strand_opts --env stringtie_strand_opt \
			--env OVERWRITE_MODE --env _SAMTOOLS_HAS_WRITE_INDEX --env _CACHED_AVAIL_MB \
			-j "$parallel_jobs" \
			--halt soon,fail,1 \
			--joblog "$HISAT2_REF_GUIDED_ROOT/parallel_hisat2_refguided_align_${BASHPID:-$$}.log" \
			_m1_align_parallel_worker {} \
			< <(printf '%s\n' "${rnaseq_list[@]}")

		local par_exit=$?
		log_info "[PARALLEL] HISAT2 Ref-Guided align+assembly complete (exit=$par_exit)"
		if [[ $par_exit -ne 0 ]]; then
			log_error "[PARALLEL] Some jobs failed - check $HISAT2_REF_GUIDED_ROOT/parallel_hisat2_refguided_align_${BASHPID:-$$}.log"
			return $par_exit
		fi
	else
		# Sequential fallback
		# Pre-compute sort params once (invariant across samples in sequential mode)
		local _seq_sort_threads=$(( THREADS < 4 ? THREADS : 4 ))
		(( _seq_sort_threads < 1 )) && _seq_sort_threads=1
		local _seq_sort_mem _seq_sort_wi_flag=""
		_samtools_sort_mem "$_seq_sort_threads" 1 _seq_sort_mem
		_samtools_has_write_index && _seq_sort_wi_flag="--write-index"
		for SRR in "${rnaseq_list[@]}"; do
			local HISAT2_DIR="$HISAT2_REF_GUIDED_ROOT/$SRR"
			mkdir -p "$HISAT2_DIR"

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ not found for $SRR - skipping"; continue; }

			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
			local out_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"

			if [[ -f "$bam" && ( -f "${bam}.bai" || -f "${bam}.csi" ) && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[ALIGN] BAM exists for $SRR/$fasta_tag - skipping"
			elif [[ -f "$out_gtf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[ALIGN] GTF exists for $SRR - skipping alignment"
			else
				log_step "Aligning: $SRR -> $fasta_tag (HISAT2 Ref-Guided)"
				local summary_file="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_alignment_summary.txt"

				# Pipe hisat2 directly into samtools sort — eliminates 10-50GB SAM intermediate per sample
				# sort_threads, sort_mem, _sort_wi_flag pre-computed before loop
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "${THREADS}" --dta $hisat2_strand_opts -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" 2>"$summary_file" \
						| samtools sort -@ "$_seq_sort_threads" -m "$_seq_sort_mem" $_seq_sort_wi_flag -o "$bam"
				else
					hisat2 -p "${THREADS}" --dta $hisat2_strand_opts -x "$index_prefix" \
						-U "$trimmed1" 2>"$summary_file" \
						| samtools sort -@ "$_seq_sort_threads" -m "$_seq_sort_mem" $_seq_sort_wi_flag -o "$bam"
				fi
				local _ps=("${PIPESTATUS[@]}")
				# Display alignment summary (strip CR for clean log output)
				if [[ -s "$summary_file" ]]; then
					while IFS= read -r _line || [[ -n "$_line" ]]; do
						_line="${_line//$'\r'/}"
						printf '%s\n' "$_line"
					done < "$summary_file"
				fi

				if [[ ${_ps[0]} -ne 0 ]]; then
					log_error "[ALIGN] HISAT2 failed for $SRR (exit=${_ps[0]}) — skipping sample"
					rm -f "$bam"
					continue
				fi
				if [[ ${_ps[1]} -ne 0 || ! -f "$bam" ]]; then
					log_error "[ALIGN] samtools sort failed for $SRR — skipping sample"
					rm -f "$bam"
					continue
				fi

				# Only run separate index if --write-index was not used
				# Cap index threads at 4 — samtools index is I/O-bound
				if [[ -z "$_seq_sort_wi_flag" ]]; then
					local _idx_t=$THREADS; (( _idx_t > 4 )) && _idx_t=4
					run_with_space_time_log samtools index -@ "$_idx_t" "$bam"
				fi

				# Infer strandness once on the first sample
				if [[ -z "$strandness" && -z "${_m1_strand_inferred:-}" ]]; then
					_m1_strand_inferred=1
					_m1_infer_strandness "$bam" "$gtf" "$HISAT2_REF_GUIDED_INDEX_DIR" "HISAT2_RG" "$SRR"
				fi
			fi

			# StringTie quantification (ref-guided, single pass)
			local out_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR"
			mkdir -p "$out_dir"

			if [[ -f "$out_gtf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[STRINGTIE] Quantification exists for $SRR/$fasta_tag - skipping"
			else
				log_step "Quantifying transcripts: $SRR -> $fasta_tag"
				run_with_space_time_log --input "$bam" --output "$out_dir" \
					stringtie -e $stringtie_strand_opt -p "$THREADS" "$bam" -G "$gtf" -o "$out_gtf" \
						-A "$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv" \
						-B -C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf" \
					|| { log_error "[STRINGTIE] Ref-guided quantification failed for $SRR"; rm -f "$out_gtf"; continue; }
			fi

			# Collect BAM metrics before potential deletion
			[[ -f "$bam" ]] && _m1_collect_bam_metrics "$bam" "$HISAT2_DIR" "HISAT2_RG" "$SRR"

			# Only delete BAM after verifying StringTie produced valid output
			if [[ "$keep_bam_global" != "y" && -f "$bam" ]]; then
				if [[ -f "$out_gtf" && -s "$out_gtf" ]]; then
					log_warn "[BAM] Deleting $SRR BAM to save disk (set keep_bam_global=y to retain)"
					rm -f "$bam" "${bam}.bai" "${bam}.csi"
				else
					log_warn "[BAM] Retaining $SRR BAM — StringTie GTF missing or empty"
				fi
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

	# Use array + printf to avoid O(n²) string reallocation from repeated += on large sample sets.
	# O(S) array appends, then single O(total_chars) printf write.
	local prepde_lines=() samples_found=0
	for SRR in "${rnaseq_list[@]}"; do
		local assembled_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		if [[ -f "$assembled_gtf" ]]; then
			prepde_lines+=("$SRR $assembled_gtf")
			(( samples_found++ )) || true
		fi
	done

	[[ $samples_found -lt 2 ]] && { log_error "Insufficient samples: $samples_found (need ≥2)"; return 1; }
	printf '%s\n' "${prepde_lines[@]}" > "$prepde_sample_list"

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

		# Cache prepDE.py availability (avoid per-call PATH scan)
		if [[ -z "${_HAS_PREPDE:-}" ]]; then
			_HAS_PREPDE=false; command -v prepDE.py >/dev/null 2>&1 && _HAS_PREPDE=true
		fi
		if [[ "$_HAS_PREPDE" == "true" ]]; then
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
	if [[ -f "$gene_count_matrix" ]]; then
		validate_count_matrix "$gene_count_matrix" "gene" 2 || \
			log_warn "Gene count matrix validation failed: $gene_count_matrix"
	fi

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

	# Cache availability probes (avoid per-worker PATH scan in parallel mode)
	if [[ -z "${_M1_HAS_INFER_EXP:-}" ]]; then
		_M1_HAS_INFER_EXP=false; command -v infer_experiment.py >/dev/null 2>&1 && _M1_HAS_INFER_EXP=true
	fi
	if [[ "$_M1_HAS_INFER_EXP" != "true" ]]; then
		_parallel_log "$method" "$srr" WARN "infer_experiment.py not found — skipping strandness check"
		return 0
	fi

	# Build BED12 from GTF (cached in index dir, rebuilt if GTF is newer)
	if [[ ! -s "$bed12" || "$gtf" -nt "$bed12" ]]; then
		if [[ -z "${_M1_HAS_GTF2BED:-}" ]]; then
			_M1_HAS_GTF2BED=false
			command -v gtfToGenePred >/dev/null 2>&1 && command -v genePredToBed >/dev/null 2>&1 && _M1_HAS_GTF2BED=true
		fi
		if [[ "$_M1_HAS_GTF2BED" == "true" ]]; then
			gtfToGenePred "$gtf" /dev/stdout 2>/dev/null | genePredToBed /dev/stdin "$bed12" 2>/dev/null
		else
			_parallel_log "$method" "$srr" WARN "gtfToGenePred not found — cannot build BED12 for strandness check"
			return 0
		fi
	fi

	[[ ! -s "$bed12" ]] && { _parallel_log "$method" "$srr" WARN "BED12 conversion from GTF failed"; return 0; }

	_parallel_log "$method" "$srr" INFO "Running infer_experiment.py on: ${bam##*/}"
	local result
	result=$(infer_experiment.py -i "$bam" -r "$bed12" 2>/dev/null)
	printf '%s\n' "$result" >> "$sentinel"

	# Parse: val_1 = FR (forward-stranded), val_2 = RF (reverse-stranded / dUTP)
	# Paired-end: "1++,1--,2+-,2-+" (FR) / "1+-,1-+,2++,2--" (RF)
	# Single-end:       "++,--"      (FR) /       "+-,-+"       (RF)
	# Single AWK pass replaces 4-8 grep pipelines (saves 8-16 process spawns)
	local val1 val2
	# Herestring avoids printf subprocess — $result is already in memory
	read -r val1 val2 < <(awk '
		/1\+\+,1--,2\+-,2-\+/ || /"\+\+,--"/ {
			match($0, /[0-9]+\.[0-9]+/); if (RSTART) v1 = substr($0, RSTART, RLENGTH)
		}
		/1\+-,1-\+,2\+\+,2--/ || /"\+-,-\+"/ {
			match($0, /[0-9]+\.[0-9]+/); if (RSTART) v2 = substr($0, RSTART, RLENGTH)
		}
		END { print (v1 ? v1 : ""), (v2 ? v2 : "") }
	' <<< "$result")

	# Bash integer arithmetic for float comparison (×1000): avoids 2 awk forks.
	# Right-pad decimal to 3 digits so 0.6→600, 0.75→750, 0.123→123.
	local rec _v1_int=0 _v2_int=0
	if [[ -n "$val1" && "$val1" == *.* ]]; then
		local _d1="${val1#*.}000"; _v1_int="${val1%%.*}${_d1:0:3}"
	fi
	if [[ -n "$val2" && "$val2" == *.* ]]; then
		local _d2="${val2#*.}000"; _v2_int="${val2%%.*}${_d2:0:3}"
	fi
	if   [[ -n "$val1" ]] && (( 10#${_v1_int} > 600 )); then
		rec="FR  → add --STRANDNESS FR to your config (ligation / forward-stranded)"
	elif [[ -n "$val2" ]] && (( 10#${_v2_int} > 600 )); then
		rec="RF  → add --STRANDNESS RF to your config (dUTP / TruSeq / reverse-stranded)"
	else
		rec="Unstranded  → omit --STRANDNESS (already the default)"
	fi

	_parallel_log "$method" "$srr" INFO "[STRANDNESS]   _val_1 (FR  1++,1--,2+-,2-+): ${val1:-N/A}"
	_parallel_log "$method" "$srr" INFO "[STRANDNESS]   _val_2 (RF  1+-,1-+,2++,2--): ${val2:-N/A}"
	_parallel_log "$method" "$srr" WARN "[STRANDNESS]   Recommendation: $rec"
}
