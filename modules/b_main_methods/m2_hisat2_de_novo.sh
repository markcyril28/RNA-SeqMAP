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
	index_prefix="$HISAT2_DE_NOVO_INDEX_DIR/${fasta_tag}_index"

	# BUILD HISAT2 INDEX
	mkdir -p "$HISAT2_DE_NOVO_INDEX_DIR"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1; then
		log_info "[HISAT2 INDEX] De novo index exists - skipping build"
	else
		log_step "Building HISAT2 de novo index from $fasta"
		log_file_size "$fasta" "Input FASTA for HISAT2 de novo index"
		run_with_space_time_log --input "$fasta" --output "$HISAT2_DE_NOVO_INDEX_DIR" \
			hisat2-build -p "${THREADS}" "$fasta" "$index_prefix"
		log_file_size "$HISAT2_DE_NOVO_INDEX_DIR" "HISAT2 de novo index output"
	fi

	# ALIGNMENT AND STRINGTIE ASSEMBLY
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		mkdir -p "$HISAT2_DIR"
		
		find_trimmed_fastq "$SRR"
		[[ -z "$trimmed1" || -z "$trimmed2" ]] && { log_warn "Trimmed FASTQ for $SRR not found - skipping"; continue; }

		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
		local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped.sam"
		
		if [[ -f "$bam" && -f "${bam}.bai" ]]; then
			log_info "[HISAT2 ALIGN] BAM exists for $SRR - skipping alignment"
		else
			log_step "Aligning $SRR using HISAT2 De Novo"
			
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$HISAT2_DIR" \
					hisat2 -p "${THREADS}" --dta -x "$index_prefix" -1 "$trimmed1" -2 "$trimmed2" -S "$sam"
			else
				run_with_space_time_log hisat2 -p "${THREADS}" --dta -x "$index_prefix" -U "$trimmed1" -S "$sam"
			fi
			
			log_info "[SAMTOOLS] Converting SAM to sorted BAM..."
			run_with_space_time_log --input "$sam" --output "$bam" samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			run_with_space_time_log samtools index -@ "${THREADS}" "$bam"
			rm -f "$sam"
		fi

		# StringTie assembly (de novo - no reference GTF)
		local out_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
		mkdir -p "$out_dir"
		
		if [[ -f "$out_gtf" ]]; then
			log_info "[STRINGTIE] De novo assembly exists for $SRR - skipping"
		else
			log_step "Assembling transcripts for $SRR (de novo)"
			run_with_space_time_log --input "$bam" --output "$out_dir" \
				stringtie -p "$THREADS" "$bam" -o "$out_gtf" \
					-A "$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo.tsv"
		fi
		
		# Cleanup BAM files if configured
		if [[ "$keep_bam_global" != "y" ]]; then
			rm -f "$bam" "${bam}.bai" "$sam"
		fi
		
		log_info "[STRINGTIE] Done processing $SRR (de novo)"
	done
	
	log_step "HISAT2 de novo pipeline completed for $fasta_tag"
}
