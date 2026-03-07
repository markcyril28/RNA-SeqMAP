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
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR"
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
					hisat2 -p "${THREADS}" --dta -x "$index_prefix" -1 "$trimmed1" -2 "$trimmed2" -S "$sam"
			else
				run_with_space_time_log hisat2 -p "${THREADS}" --dta -x "$index_prefix" -U "$trimmed1" -S "$sam"
			fi
			
			log_info "[SAMTOOLS] Converting to sorted BAM..."
			run_with_space_time_log --input "$sam" --output "$bam" samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			run_with_space_time_log samtools index -@ "${THREADS}" "$bam"
			rm -f "$sam"
		fi

		# StringTie assembly
		local out_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		local ballgown_dir="$out_dir/ballgown"
		mkdir -p "$out_dir" "$ballgown_dir"
		
		if [[ -f "$out_gtf" ]]; then
			log_info "[STRINGTIE] Assembly exists for $SRR/$fasta_tag - skipping"
		else
			log_step "Assembling transcripts: $SRR -> $fasta_tag"
			# -e: only estimate abundance for transcripts in reference GTF (no novel transcripts)
			# This preserves original gene_id values instead of creating MSTRG.* IDs
			run_with_space_time_log --input "$bam" --output "$out_dir" \
				stringtie -e -p "$THREADS" "$bam" -G "$gtf" -o "$out_gtf" \
					-A "$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv" \
					-B -C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf"
		fi
	done
	
	# MERGE GTF FILES
	log_step "Creating merged GTF file"
	local merge_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/merged"
	local merged_gtf="$merge_dir/${fasta_tag}_ref_guided_merged.gtf"
	local gtf_list="$merge_dir/gtf_list.txt"
	mkdir -p "$merge_dir"
	
	true > "$gtf_list"
	for SRR in "${rnaseq_list[@]}"; do
		local out_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		[[ -f "$out_gtf" ]] && echo "$out_gtf" >> "$gtf_list"
	done
	
	if [[ ! -f "$merged_gtf" ]]; then
		log_info "[STRINGTIE MERGE] Merging GTF files..."
		run_with_space_time_log stringtie --merge -p "$THREADS" -G "$gtf" -o "$merged_gtf" "$gtf_list"
	fi
	
	# RE-ESTIMATE ABUNDANCES WITH MERGED GTF
	for SRR in "${rnaseq_list[@]}"; do
		local bam="$HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
		local final_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/final"
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
	
	# PREPARE COUNT MATRICES
	log_step "Preparing count matrices for DESeq2"
	local deseq2_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/deseq2_input"
	local prepde_sample_list="$deseq2_dir/sample_list.txt"
	local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
	local transcript_count_matrix="$deseq2_dir/transcript_count_matrix.csv"
	mkdir -p "$deseq2_dir"
	
	true > "$prepde_sample_list"
	local samples_found=0
	for SRR in "${rnaseq_list[@]}"; do
		local final_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/final/${SRR}_${fasta_tag}_ref_guided_final.gtf"
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
