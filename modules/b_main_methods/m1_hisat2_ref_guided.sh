#!/bin/bash
# ==============================================================================
# METHOD 1: HISAT2 REFERENCE GUIDED PIPELINE
# ==============================================================================
# HISAT2 Reference Guided alignment with StringTie assembly
# Uses reference GTF for splice site information
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M1_HISAT2_REF_SOURCED:-}" == "true" ]] && return 0
export M1_HISAT2_REF_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# HISAT2 REFERENCE GUIDED PIPELINE
# ==============================================================================

hisat2_ref_guided_pipeline() {
	local fasta="" gtf="" rnaseq_list=()
	
	# Parse arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--GTF) gtf="$2"; shift 2;;
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

	local fasta_base fasta_tag index_prefix
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	index_prefix="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_ref_guided"

	# BUILD HISAT2 REFERENCE-GUIDED INDEX
	mkdir -p "$HISAT2_REF_GUIDED_INDEX_DIR"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1; then
		log_info "[INDEX] Ref-Guided index exists - skipping build"
	else
		log_step "Building HISAT2 Ref-Guided index: $fasta_base"
		
		local splice_sites="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_splice_sites.txt"
		local exons="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_exons.txt"
		local build_opts=""
		
		# Extract splice sites and exons (may be empty for single-exon transcriptomes)
		hisat2_extract_splice_sites.py "$gtf" > "$splice_sites" 2>/dev/null || true
		hisat2_extract_exons.py "$gtf" > "$exons" 2>/dev/null || true

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

	# ALIGNMENT AND STRINGTIE ASSEMBLY
	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]]; then
		log_step "[PARALLEL] HISAT2 Ref-Guided Align+StringTie: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export fasta_tag index_prefix threads_per_job
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

			if [[ -f "$bam" && -f "${bam}.bai" ]]; then
				_parallel_log HISAT2_RG "$SRR" INFO "BAM exists - skipping alignment"
			else
				_parallel_log HISAT2_RG "$SRR" INFO "Aligning with $threads_per_job threads"
				local align_exit=0
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "$threads_per_job" -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" -S "$sam" 2>&1 || align_exit=$?
				else
					hisat2 -p "$threads_per_job" -x "$index_prefix" \
						-U "$trimmed1" -S "$sam" 2>&1 || align_exit=$?
				fi
				[[ $align_exit -ne 0 ]] && { _parallel_log HISAT2_RG "$SRR" ERROR "HISAT2 failed (exit=$align_exit)"; return $align_exit; }

				samtools sort -@ "$threads_per_job" -o "$bam" "$sam" 2>&1 || { _parallel_log HISAT2_RG "$SRR" ERROR "samtools sort failed"; return 1; }
				samtools index -@ "$threads_per_job" "$bam" 2>&1 || true
				rm -f "$sam"
			fi

			# StringTie assembly (ref-guided)
			local out_dir="$abs_stringtie_rg_root/$SRR"
			local ballgown_dir="$out_dir/ballgown"
			local out_gtf="$out_dir/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
			mkdir -p "$out_dir" "$ballgown_dir"

			if [[ -f "$out_gtf" ]]; then
				_parallel_log HISAT2_RG "$SRR" INFO "Assembly exists - skipping"
			else
				_parallel_log HISAT2_RG "$SRR" INFO "Assembling transcripts (ref-guided)"
				stringtie -e -p "$threads_per_job" "$bam" -G "$abs_gtf" -o "$out_gtf" \
					-A "$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv" \
					-B -C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf" 2>&1 || \
					{ _parallel_log HISAT2_RG "$SRR" ERROR "StringTie failed"; return 1; }
			fi

			_parallel_log HISAT2_RG "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m1_align_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env index_prefix --env threads_per_job \
			--env abs_hisat2_rg_root --env abs_stringtie_rg_root --env abs_gtf \
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

			if [[ -f "$bam" && -f "${bam}.bai" ]]; then
				log_info "[ALIGN] BAM exists for $SRR/$fasta_tag - skipping"
			else
				log_step "Aligning: $SRR -> $fasta_tag (HISAT2 Ref-Guided)"

				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$HISAT2_DIR" \
						hisat2 -p "${THREADS}" -x "$index_prefix" -1 "$trimmed1" -2 "$trimmed2" -S "$sam"
				else
					run_with_space_time_log hisat2 -p "${THREADS}" -x "$index_prefix" -U "$trimmed1" -S "$sam"
				fi

				log_info "[SAMTOOLS] Converting to sorted BAM..."
				run_with_space_time_log --input "$sam" --output "$bam" samtools sort -@ "${THREADS}" -o "$bam" "$sam"
				run_with_space_time_log samtools index -@ "${THREADS}" "$bam"
				rm -f "$sam"
			fi

			# StringTie assembly
			local out_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR"
			local out_gtf="$out_dir/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
			local ballgown_dir="$out_dir/ballgown"
			mkdir -p "$out_dir" "$ballgown_dir"

			if [[ -f "$out_gtf" ]]; then
				log_info "[STRINGTIE] Assembly exists for $SRR/$fasta_tag - skipping"
			else
				log_step "Assembling transcripts: $SRR -> $fasta_tag"
				run_with_space_time_log --input "$bam" --output "$out_dir" \
					stringtie -e -p "$THREADS" "$bam" -G "$gtf" -o "$out_gtf" \
						-A "$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv" \
						-B -C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf"
			fi
		done
	fi
	
	# MERGE GTF FILES
	log_step "Creating merged GTF file"
	local merge_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/merged"
	local merged_gtf="$merge_dir/${fasta_tag}_ref_guided_merged.gtf"
	local gtf_list="$merge_dir/gtf_list.txt"
	mkdir -p "$merge_dir"
	
	true > "$gtf_list"
	for SRR in "${rnaseq_list[@]}"; do
		local out_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		[[ -f "$out_gtf" ]] && echo "$out_gtf" >> "$gtf_list"
	done
	
	if [[ ! -f "$merged_gtf" ]]; then
		log_info "[STRINGTIE MERGE] Merging GTF files..."
		run_with_space_time_log stringtie --merge -p "$THREADS" -G "$gtf" -o "$merged_gtf" "$gtf_list"
	fi
	
	# RE-ESTIMATE ABUNDANCES WITH MERGED GTF
	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]]; then
		log_step "[PARALLEL] Re-estimating abundances: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export fasta_tag threads_per_job
		local abs_hisat2_rg_root="$HISAT2_REF_GUIDED_ROOT"
		[[ "$abs_hisat2_rg_root" != /* ]] && abs_hisat2_rg_root="$(pwd)/$abs_hisat2_rg_root"
		local abs_stringtie_rg_root="$STRINGTIE_HISAT2_REF_GUIDED_ROOT"
		[[ "$abs_stringtie_rg_root" != /* ]] && abs_stringtie_rg_root="$(pwd)/$abs_stringtie_rg_root"
		local abs_merged_gtf="$merged_gtf"
		[[ "$abs_merged_gtf" != /* ]] && abs_merged_gtf="$(pwd)/$abs_merged_gtf"
		export abs_hisat2_rg_root abs_stringtie_rg_root abs_merged_gtf

		_m1_reestimate_parallel_worker() {
			local SRR="$1"
			# Conda reactivation
			if [[ -n "${CONDA_PREFIX:-}" ]]; then
				eval "$(conda shell.bash hook 2>/dev/null)" && conda activate "${CONDA_DEFAULT_ENV:-base}" 2>/dev/null || true
			fi

			local bam="$abs_hisat2_rg_root/$SRR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
			local final_dir="$abs_stringtie_rg_root/$SRR/final"
			local final_gtf="$final_dir/${SRR}_${fasta_tag}_ref_guided_final.gtf"
			mkdir -p "$final_dir"

			if [[ ! -f "$bam" && -f "$final_gtf" ]]; then
				_parallel_log HISAT2_RG_RE "$SRR" INFO "Already complete (no BAM, final exists) - skipping"
				return 0
			fi

			if [[ -f "$bam" ]]; then
				_parallel_log HISAT2_RG_RE "$SRR" INFO "Re-estimating abundances with $threads_per_job threads"
				stringtie -p "$threads_per_job" -e -B -G "$abs_merged_gtf" \
					-A "$final_dir/${SRR}_${fasta_tag}_ref_guided_final_abundances.tsv" \
					-o "$final_gtf" "$bam" 2>&1 || \
					{ _parallel_log HISAT2_RG_RE "$SRR" ERROR "StringTie re-estimation failed"; return 1; }

				[[ "$keep_bam_global" != "y" ]] && rm -f "$bam" "${bam}.bai"
			fi

			_parallel_log HISAT2_RG_RE "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m1_reestimate_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env threads_per_job \
			--env abs_hisat2_rg_root --env abs_stringtie_rg_root --env abs_merged_gtf \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$HISAT2_REF_GUIDED_ROOT/parallel_hisat2_refguided_reestimate.log" \
			_m1_reestimate_parallel_worker {}

		local par_exit2=$?
		log_info "[PARALLEL] Re-estimation complete (exit=$par_exit2)"
		[[ $par_exit2 -ne 0 ]] && log_warn "[PARALLEL] Some re-estimation jobs failed - check $HISAT2_REF_GUIDED_ROOT/parallel_hisat2_refguided_reestimate.log"
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local bam="$HISAT2_REF_GUIDED_ROOT/$SRR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
			local final_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR/final"
			local final_gtf="$final_dir/${SRR}_${fasta_tag}_ref_guided_final.gtf"
			mkdir -p "$final_dir"

			if [[ ! -f "$bam" && -f "$final_gtf" ]]; then continue; fi

			if [[ -f "$bam" ]]; then
				log_step "Re-estimating abundances for $SRR"
				run_with_space_time_log stringtie -p "$THREADS" -e -B -G "$merged_gtf" \
					-A "$final_dir/${SRR}_${fasta_tag}_ref_guided_final_abundances.tsv" \
					-o "$final_gtf" "$bam"

				[[ "$keep_bam_global" != "y" ]] && rm -f "$bam" "${bam}.bai"
			fi
		done
	fi
	
	# PREPARE COUNT MATRICES
	log_step "Preparing count matrices for DESeq2"
	local deseq2_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/deseq2_input"
	local prepde_sample_list="$deseq2_dir/sample_list.txt"
	local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
	local transcript_count_matrix="$deseq2_dir/transcript_count_matrix.csv"
	mkdir -p "$deseq2_dir"
	
	true > "$prepde_sample_list"
	local samples_found=0
	for SRR in "${rnaseq_list[@]}"; do
		local final_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$SRR/final/${SRR}_${fasta_tag}_ref_guided_final.gtf"
		if [[ -f "$final_gtf" ]]; then
			echo "$SRR $final_gtf" >> "$prepde_sample_list"
			((samples_found++))
		fi
	done
	
	[[ $samples_found -lt 2 ]] && { log_error "Insufficient samples: $samples_found (need ≥2)"; return 1; }
	
	if [[ ! -f "$gene_count_matrix" ]]; then
		# Detect read length
		local read_length=150
		for SRR in "${rnaseq_list[@]}"; do
			find_trimmed_fastq "$SRR"
			if [[ -n "$trimmed1" ]]; then
				read_length=$(detect_read_length "$trimmed1" 150)
				break
			fi
		done
		
		if command -v prepDE.py >/dev/null 2>&1; then
			run_with_space_time_log prepDE.py -i "$prepde_sample_list" \
				-g "$gene_count_matrix" -t "$transcript_count_matrix" -l "$read_length"
		else
			log_error "prepDE.py not found"
			return 1
		fi
	fi
	
	# Create sample metadata
	local sample_metadata="$deseq2_dir/sample_metadata.csv"
	[[ ! -f "$sample_metadata" ]] && create_sample_metadata "$sample_metadata" rnaseq_list[@]
	
	# Validate outputs
	[[ -f "$gene_count_matrix" ]] && validate_count_matrix "$gene_count_matrix" "gene" 2
	
	log_step "HISAT2 reference-guided pipeline completed for $fasta_tag"
	log_info "Merged GTF: $merged_gtf"
	log_info "Gene count matrix: $gene_count_matrix"
	log_info "Sample metadata: $sample_metadata"
}
