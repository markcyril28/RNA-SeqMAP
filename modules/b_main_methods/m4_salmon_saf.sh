#!/bin/bash
# ==============================================================================
# METHOD 4: SALMON SAF (SELECTIVE ALIGNMENT WITH DECOYS) PIPELINE
# ==============================================================================
# Quantify expression using decoy-aware Salmon (Selective Alignment)
# Fast and accurate pseudo-alignment
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M4_SALMON_SOURCED:-}" == "true" ]] && return 0
export M4_SALMON_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# SALMON CONFIGURATION - IMPORTANT PARAMETERS AT TOP
# ==============================================================================

# Salmon quantification settings
SALMON_KMER_SIZE="${SALMON_KMER_SIZE:-31}"
SALMON_NUM_BOOTSTRAPS="${SALMON_NUM_BOOTSTRAPS:-30}"

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
	local work="tmp_${tag}_gentrome"
	local idx_dir="$SALMON_INDEX_ROOT/${tag}_decoySAF"
	local quant_root="$SALMON_QUANT_ROOT/$tag"
	local matrix_dir="$SALMON_MATRIX_ROOT/$tag"

	mkdir -p "$SALMON_INDEX_ROOT" "$quant_root" "$matrix_dir" "$work"

	# BUILD DECOY-AWARE INDEX
	if [[ -f "$idx_dir/versionInfo.json" ]]; then
		log_info "[SALMON INDEX] Decoy index already exists. Skipping."
	else
		log_step "Building decoy-aware Salmon index for $tag"
		awk '/^>/{print substr($0,2); next}{next}' "$genome" > "$work/decoys.txt"
		cat "$fasta" "$genome" > "$work/gentrome.fa"
		log_file_size "$work/gentrome.fa" "Gentrome FASTA for Salmon - $tag"
		log_file_size "$work/decoys.txt" "Decoy list for Salmon - $tag"
		run_with_space_time_log --input "$work" --output "$idx_dir" salmon index \
			-t "$work/gentrome.fa" \
			-d "$work/decoys.txt" \
			-i "$idx_dir" \
			-k "$SALMON_KMER_SIZE" -p "$THREADS"
		log_file_size "$idx_dir" "Salmon index output - $tag"
		log_info "[CLEANUP] Removing temporary gentrome work directory"
		rm -rf "$work"
	fi

	# QUANTIFICATION PER SRR
	for SRR in "${rnaseq_list[@]}"; do
		local out_dir="$quant_root/$SRR"
		mkdir -p "$out_dir"
		
		[[ -f "$out_dir/quant.sf" ]] && { log_info "[SALMON QUANT] Quantification for $SRR already exists. Skipping."; continue; }

		find_trimmed_fastq "$SRR"
		[[ -z "$trimmed1" ]] && { log_warn "Missing trimmed reads for $SRR. Skipping."; continue; }

		log_step "Quantifying expression for $SRR with Salmon"
		
		if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
			# Paired-end reads
			log_info "[SALMON QUANT] Using paired-end reads for $SRR"
			run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$out_dir" salmon quant \
				-i "$idx_dir" -l A \
				-1 "$trimmed1" -2 "$trimmed2" \
				-p "$THREADS" \
				--validateMappings \
				--seqBias --gcBias --posBias \
				--numBootstraps "$SALMON_NUM_BOOTSTRAPS" \
				-o "$out_dir"
		else
			# Single-end reads
			log_info "[SALMON QUANT] Using single-end reads for $SRR"
			run_with_space_time_log salmon quant \
				-i "$idx_dir" -l A \
				-r "$trimmed1" \
				-p "$THREADS" \
				--validateMappings \
				--seqBias --gcBias --posBias \
				--numBootstraps "$SALMON_NUM_BOOTSTRAPS" \
				-o "$out_dir"
		fi
		log_file_size "$out_dir/quant.sf" "Salmon quantification output - $SRR"
	done

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
	local -n srr_list=$5

	log_step "Generating gene and transcript matrices (Salmon)"
	
	# Check or create gene_trans_map
	local gene_trans_map="${fasta}.gene_trans_map"
	if [[ ! -f "$gene_trans_map" ]]; then
		log_info "[SALMON MATRIX] Creating gene-transcript mapping file..."
		_create_gene_trans_map "$fasta" "$gene_trans_map"
	fi
	
	# Generate matrices using Trinity's script or manual creation
	if command -v abundance_estimates_to_matrix.pl >/dev/null 2>&1; then
		log_info "[SALMON MATRIX] Running abundance_estimates_to_matrix.pl..."
		run_with_space_time_log abundance_estimates_to_matrix.pl \
			--est_method salmon \
			--gene_trans_map "$gene_trans_map" \
			--out_prefix "$matrix_dir/genes" \
			--name_sample_by_basedir "$quant_root"/*/quant.sf
	else
		log_warn "abundance_estimates_to_matrix.pl not found. Creating manual count matrix..."
		_create_manual_salmon_matrix "$quant_root" "$matrix_dir" srr_list
	fi
	
	# Prepare DESeq2-compatible outputs
	_prepare_salmon_deseq2_output "$tag" "$quant_root" "$matrix_dir" srr_list
}

_create_gene_trans_map() {
	local fasta="$1"
	local output="$2"
	
	local is_trinity=false
	grep -q "^>TRINITY_" "$fasta" 2>/dev/null && is_trinity=true
	
	if [[ "$is_trinity" == "true" ]]; then
		log_info "[SALMON MATRIX] Detected Trinity assembly format"
		awk '/^>/ {
			trans = $1; gsub(/^>/, "", trans)
			gene = trans
			if (match(gene, /^(.+)_i[0-9]+$/, arr)) { gene = arr[1] }
			print gene "\t" trans
		}' "$fasta" > "$output"
	else
		grep "^>" "$fasta" | sed 's/^>//' | awk '{
			trans=$1
			if (match($0, /gene=([^ ]+)/, arr)) { gene=arr[1] }
			else if (match($0, /gene_id[=:]([^ ]+)/, arr)) { gene=arr[1] }
			else if (match(trans, /^([^|]+)\|/, arr)) { gene=arr[1] }
			else { gene=trans }
			print gene "\t" trans
		}' > "$output"
	fi
	
	local unique_genes=$(cut -f1 "$output" | sort -u | wc -l)
	local total_transcripts=$(wc -l < "$output")
	log_info "[SALMON MATRIX] Created gene-transcript map: $unique_genes genes, $total_transcripts transcripts"
}

_create_manual_salmon_matrix() {
	local quant_root="$1"
	local matrix_dir="$2"
	local -n srr_list=$3
	
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
				local num_genes=$(wc -l < "$temp_gene_ids")
				yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_counts.tmp"
			fi
		done
		
		paste "$temp_gene_ids" "$matrix_dir"/*_counts.tmp > "$temp_counts"
		
		echo -n "gene_id" > "$matrix_dir/genes.counts.matrix"
		for SRR in "${srr_list[@]}"; do
			echo -ne "\t$SRR" >> "$matrix_dir/genes.counts.matrix"
		done
		echo "" >> "$matrix_dir/genes.counts.matrix"
		cat "$temp_counts" >> "$matrix_dir/genes.counts.matrix"
		
		rm -f "$temp_gene_ids" "$temp_counts" "$matrix_dir"/*_counts.tmp
	fi
}

_prepare_salmon_deseq2_output() {
	local tag="$1"
	local quant_root="$2"
	local matrix_dir="$3"
	local -n srr_list=$4

	log_step "Preparing DESeq2-compatible count matrix for Salmon pipeline"
	
	local deseq2_dir="$matrix_dir/deseq2_input"
	local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
	local sample_metadata="$deseq2_dir/sample_metadata.csv"
	mkdir -p "$deseq2_dir"
	
	# Verify quantifications
	local quant_count=0
	for SRR in "${srr_list[@]}"; do
		[[ -f "$quant_root/$SRR/quant.sf" ]] && ((quant_count++))
	done
	
	[[ $quant_count -lt 2 ]] && { log_error "Insufficient Salmon quantifications (found: $quant_count, need: ≥2)"; return 1; }
	log_info "[SALMON] Found $quant_count samples with successful quantifications"
	
	# Convert to CSV
	if [[ -f "$matrix_dir/genes.counts.matrix" && ! -f "$gene_count_matrix" ]]; then
		log_info "[SALMON MATRIX] Converting count matrix to CSV format..."
		sed 's/\t/,/g' "$matrix_dir/genes.counts.matrix" | sed '1s/gene_id/Gene_ID/' > "$gene_count_matrix"
	fi
	
	# Create sample metadata
	[[ ! -f "$sample_metadata" ]] && create_sample_metadata "$sample_metadata" srr_list
	
	# Generate tximport script
	local tximport_script="$deseq2_dir/run_tximport_salmon.R"
	[[ ! -f "$tximport_script" ]] && generate_tximport_script "salmon" "$quant_root" "$tximport_script" "$sample_metadata"
	
	# Create TPM matrix
	if [[ -f "$matrix_dir/genes.TPM.not_cross_norm" ]]; then
		local tpm_matrix="$deseq2_dir/gene_tpm_matrix.csv"
		[[ ! -f "$tpm_matrix" ]] && sed 's/\t/,/g' "$matrix_dir/genes.TPM.not_cross_norm" | sed '1s/gene_id/Gene_ID/' > "$tpm_matrix"
	fi
	
	# Create summary
	_create_salmon_summary "$tag" "$quant_root" "$deseq2_dir" srr_list
	
	# Validate
	[[ -f "$gene_count_matrix" ]] && validate_count_matrix "$gene_count_matrix" "gene" 2
	
	log_info "DESeq2 input files:"
	log_info "  - Gene count matrix: $gene_count_matrix"
	log_info "  - Sample metadata: $sample_metadata"
}

_create_salmon_summary() {
	local tag="$1"
	local quant_root="$2"
	local deseq2_dir="$3"
	local -n srr_list=$4
	
	local summary_file="$deseq2_dir/salmon_summary.txt"
	[[ -f "$summary_file" ]] && return 0
	
	{
		echo "==================================================================="
		echo "Salmon SAF Quantification Summary for $tag"
		echo "==================================================================="
		echo "Date: $(date)"
		echo "Samples processed: ${#srr_list[@]}"
		echo "Method: Salmon Selective Alignment with decoy-aware indexing"
		echo ""
		echo "Per-sample statistics:"
		echo "-------------------------------------------------------------------"
		
		for SRR in "${srr_list[@]}"; do
			if [[ -f "$quant_root/$SRR/quant.sf" ]]; then
				local total=$(awk 'NR>1' "$quant_root/$SRR/quant.sf" | wc -l)
				local expressed=$(awk 'NR>1 && $5>0' "$quant_root/$SRR/quant.sf" | wc -l)
				local reads=$(awk 'NR>1 {sum+=$5} END {print int(sum)}' "$quant_root/$SRR/quant.sf")
				echo "$SRR: $expressed/$total expressed, $reads total counts"
			fi
		done
	} > "$summary_file"
	
	log_info "[SALMON SUMMARY] Summary saved to: $summary_file"
}
