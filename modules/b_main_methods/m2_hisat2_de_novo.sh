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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	set_fasta_output_dirs "$fasta_tag"
	index_prefix="$HISAT2_DE_NOVO_INDEX_DIR/${fasta_tag}_index"

	# BUILD HISAT2 INDEX
	mkdir -p "$HISAT2_DE_NOVO_INDEX_DIR"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1 && [[ "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[HISAT2 INDEX] De novo index exists - skipping build"
	else
		log_step "Building HISAT2 de novo index from $fasta"
		log_file_size "$fasta" "Input FASTA for HISAT2 de novo index"
		run_with_space_time_log --input "$fasta" --output "$HISAT2_DE_NOVO_INDEX_DIR" \
			hisat2-build -p "${THREADS}" "$fasta" "$index_prefix" \
			|| { log_error "HISAT2 de novo index build failed for $fasta_tag"; rm -f "${index_prefix}".*.ht2; return 1; }
		log_file_size "$HISAT2_DE_NOVO_INDEX_DIR" "HISAT2 de novo index output"
	fi

	# ALIGNMENT AND STRINGTIE ASSEMBLY
	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]] && [[ "${USE_GNU_PARALLEL:-TRUE}" != "FALSE" ]]; then
		log_step "[PARALLEL] HISAT2 De Novo: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export fasta_tag index_prefix threads_per_job
		local abs_hisat2_dn_root="$HISAT2_DE_NOVO_ROOT"
		[[ "$abs_hisat2_dn_root" != /* ]] && abs_hisat2_dn_root="$(pwd)/$abs_hisat2_dn_root"
		local abs_stringtie_dn_root="$STRINGTIE_HISAT2_DE_NOVO_ROOT"
		[[ "$abs_stringtie_dn_root" != /* ]] && abs_stringtie_dn_root="$(pwd)/$abs_stringtie_dn_root"
		export abs_hisat2_dn_root abs_stringtie_dn_root

		_m2_align_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log HISAT2_DN "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local HISAT2_DIR="$abs_hisat2_dn_root/$SRR"
			mkdir -p "$HISAT2_DIR"
			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
			local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped.sam"

			if [[ -f "$bam" && -f "${bam}.bai" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_DN "$SRR" INFO "BAM exists - skipping alignment"
			else
				_parallel_log HISAT2_DN "$SRR" INFO "Aligning with $threads_per_job threads"
				local align_exit=0
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					hisat2 -p "$threads_per_job" -x "$index_prefix" \
						-1 "$trimmed1" -2 "$trimmed2" -S "$sam" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
					align_exit=${PIPESTATUS[0]}
				else
					hisat2 -p "$threads_per_job" -x "$index_prefix" \
						-U "$trimmed1" -S "$sam" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
					align_exit=${PIPESTATUS[0]}
				fi
				[[ $align_exit -ne 0 ]] && { _parallel_log HISAT2_DN "$SRR" ERROR "HISAT2 failed (exit=$align_exit)"; rm -f "$sam"; return $align_exit; }

				samtools sort -@ "$threads_per_job" -o "$bam" "$sam" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				[[ ${PIPESTATUS[0]} -ne 0 ]] && { _parallel_log HISAT2_DN "$SRR" ERROR "samtools sort failed"; rm -f "$sam" "$bam"; return 1; }
				samtools index -@ "$threads_per_job" "$bam" 2>&1 | sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				[[ ${PIPESTATUS[0]} -ne 0 ]] && { _parallel_log HISAT2_DN "$SRR" ERROR "samtools index failed"; rm -f "$sam" "$bam"; return 1; }
				rm -f "$sam"
			fi

			# StringTie assembly (de novo)
			local out_dir="$abs_stringtie_dn_root/$SRR"
			local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
			local out_abund="$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo.tsv"
			mkdir -p "$out_dir"

			if [[ -f "$out_gtf" && -f "$out_abund" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				_parallel_log HISAT2_DN "$SRR" INFO "De novo assembly exists - skipping"
			else
				_parallel_log HISAT2_DN "$SRR" INFO "Assembling transcripts (de novo)"
				stringtie -p "$threads_per_job" "$bam" -o "$out_gtf" \
					-A "$out_abund" 2>&1 || \
					{ _parallel_log HISAT2_DN "$SRR" ERROR "StringTie failed"; rm -f "$out_gtf" "$out_abund"; return 1; }
			fi

			if [[ "$keep_bam_global" != "y" ]]; then
				rm -f "$bam" "${bam}.bai"
			fi
			_parallel_log HISAT2_DN "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m2_align_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env index_prefix --env threads_per_job \
			--env abs_hisat2_dn_root --env abs_stringtie_dn_root \
			--env OVERWRITE_MODE \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$HISAT2_DE_NOVO_ROOT/parallel_hisat2_denovo.log" \
			_m2_align_parallel_worker {}

		local par_exit=$?
		log_info "[PARALLEL] HISAT2 De Novo complete (exit=$par_exit)"
		[[ $par_exit -ne 0 ]] && log_warn "[PARALLEL] Some jobs failed - check $HISAT2_DE_NOVO_ROOT/parallel_hisat2_denovo.log" && return $par_exit
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/$SRR"
			mkdir -p "$HISAT2_DIR"

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ for $SRR not found - skipping"; continue; }

			local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
			local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped.sam"

			if [[ -f "$bam" && -f "${bam}.bai" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[HISAT2 ALIGN] BAM exists for $SRR - skipping alignment"
			else
				log_step "Aligning $SRR using HISAT2 De Novo"

				local align_exit=0
				if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
					run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$HISAT2_DIR" \
						hisat2 -p "${THREADS}" -x "$index_prefix" -1 "$trimmed1" -2 "$trimmed2" -S "$sam"
				else
					run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$HISAT2_DIR" \
						hisat2 -p "${THREADS}" -x "$index_prefix" -U "$trimmed1" -S "$sam"
				fi
				align_exit=$?
				[[ $align_exit -ne 0 ]] && { log_error "[HISAT2] Alignment failed for $SRR (exit=$align_exit)"; rm -f "$sam"; continue; }

				log_info "[SAMTOOLS] Converting SAM to sorted BAM..."
				run_with_space_time_log --input "$sam" --output "$bam" samtools sort -@ "${THREADS}" -o "$bam" "$sam" \
					|| { log_error "[SAMTOOLS] sort failed for $SRR"; rm -f "$sam" "$bam"; continue; }
				run_with_space_time_log samtools index -@ "${THREADS}" "$bam" \
					|| { log_error "[SAMTOOLS] index failed for $SRR"; rm -f "$sam" "$bam"; continue; }
				rm -f "$sam"
			fi

			# StringTie assembly (de novo - no reference GTF)
			local out_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$SRR"
			local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
			local out_abund="$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo.tsv"
			mkdir -p "$out_dir"

			if [[ -f "$out_gtf" && -f "$out_abund" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				log_info "[STRINGTIE] De novo assembly exists for $SRR - skipping"
			else
				log_step "Assembling transcripts for $SRR (de novo)"
				run_with_space_time_log --input "$bam" --output "$out_dir" \
					stringtie -p "$THREADS" "$bam" -o "$out_gtf" \
						-A "$out_abund" \
					|| { log_error "[STRINGTIE] Assembly failed for $SRR"; rm -f "$out_gtf" "$out_abund"; continue; }
			fi

			# Cleanup BAM files if configured
			if [[ "$keep_bam_global" != "y" ]]; then
				rm -f "$bam" "${bam}.bai"
			fi

			log_info "[STRINGTIE] Done processing $SRR (de novo)"
		done
	fi
	
	log_step "HISAT2 de novo pipeline completed for $fasta_tag"
}
