#!/bin/bash
# ==============================================================================
# MAIN METHODS SHARED UTILITIES
# ==============================================================================
# Common helper functions used across all GEA analysis methods
# Sourced by: All method scripts in b_main_methods/
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${METHOD_SHARED_SOURCED:-}" == "true" ]] && return 0
export METHOD_SHARED_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/global_config_method.sh"
source "$SCRIPT_DIR/../logging/logging_utils.sh"
source "$SCRIPT_DIR/../a_preprocessing/shared_utils_preproc.sh"

# ==============================================================================
# VALIDATION FUNCTIONS
# ==============================================================================

# Validate count matrix quality
validate_count_matrix() {
	local matrix="$1"
	local matrix_type="${2:-gene}"
	local min_samples="${3:-2}"
	
	[[ ! -f "$matrix" ]] && { log_error "Matrix not found: $matrix"; return 1; }
	
	log_step "Validating $matrix_type matrix: $(basename "$matrix")"
	
	local delim="\t"
	[[ "$matrix" == *.csv ]] && delim=","
	
	local header=$(head -n1 "$matrix")
	local num_samples=$(echo "$header" | awk -F"$delim" '{print NF-1}')
	
	[[ $num_samples -lt $min_samples ]] && { log_error "Insufficient samples: $num_samples (need ≥$min_samples)"; return 1; }
	
	local total_rows=$(tail -n +2 "$matrix" | wc -l)
	local unique_ids=$(tail -n +2 "$matrix" | cut -d"${delim:0:1}" -f1 | sort -u | wc -l)
	
	[[ $total_rows -ne $unique_ids ]] && { log_error "Duplicate ${matrix_type} IDs detected!"; return 1; }
	
	log_info "[VALIDATION] Passed - Samples: $num_samples, ${matrix_type^}s: $unique_ids"
	return 0
}

# ==============================================================================
# METADATA FUNCTIONS
# ==============================================================================

# Load sample metadata from file
load_sample_metadata() {
	local metadata_file="${1:-${SAMPLE_CONDITIONS_FILE:-sample_conditions.txt}}"
	local -n metadata_array=$2
	
	[[ ! -f "$metadata_file" ]] && { log_warn "Metadata file not found: $metadata_file"; return 1; }
	
	# Count valid sample lines (exclude comments, blanks, and header)
	local line_count=$(grep -v '^#' "$metadata_file" | grep -v '^$' | grep -v '^SRR_ID' | wc -l)
	[[ $line_count -lt 2 ]] && { log_error "Metadata file must contain at least 2 samples (found: $line_count)"; return 1; }
	
	while IFS=$'\t' read -r srr condition batch; do
		[[ -z "$srr" || "$srr" =~ ^# || "$srr" == "SRR_ID" ]] && continue
		metadata_array["${srr}_condition"]="$condition"
		metadata_array["${srr}_batch"]="$batch"
	done < "$metadata_file"
	
	log_info "Loaded metadata: $(( ${#metadata_array[@]} / 2 )) samples from $metadata_file"
	return 0
}

# Create sample metadata CSV for DESeq2
create_sample_metadata() {
	local metadata_file="$1"
	local -a sample_list=("${!2}")
	local metadata_source="${3:-${SAMPLE_CONDITIONS_FILE:-sample_conditions.txt}}"
	
	local delim=","
	[[ "$metadata_file" == *.tsv ]] && delim=$'\t'
	
	declare -A sample_metadata
	local has_external=false
	load_sample_metadata "$metadata_source" sample_metadata 2>/dev/null && has_external=true
	
	echo -e "sample${delim}condition${delim}batch" > "$metadata_file"
	for SRR in "${sample_list[@]}"; do
		if [[ "$has_external" == "true" ]]; then
			local condition="${sample_metadata[${SRR}_condition]:-unknown}"
			local batch="${sample_metadata[${SRR}_batch]:-1}"
		else
			local condition="treatment"
			local batch="1"
		fi
		echo -e "$SRR${delim}$condition${delim}$batch" >> "$metadata_file"
	done
	
	[[ "$has_external" == "false" ]] && \
		log_warn "Sample conditions need manual specification in: $metadata_file"
}

# ==============================================================================
# TXIMPORT SCRIPT GENERATION
# ==============================================================================

# Get the helper scripts directory
HELPERS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/helpers" && pwd)"

# Run tximport using external R helper script
# Usage: run_tximport <method> <quant_dir> <metadata_file> [output_dir]
run_tximport() {
	local method="$1"
	local quant_dir="$2"
	local metadata_file="$3"
	local output_dir="${4:-$(dirname "$metadata_file")}"
	local helper_script="$HELPERS_DIR/tximport_helper.R"
	
	if [[ ! -f "$helper_script" ]]; then
		log_error "tximport_helper.R not found: $helper_script"
		return 1
	fi
	
	log_info "[TXIMPORT] Running $method import..."
	Rscript "$helper_script" "$method" "$quant_dir" "$metadata_file" "$output_dir"
}

# Generate tximport R script (legacy compatibility - copies helper)
generate_tximport_script() {
	local method="$1"
	local quant_dir="$2"
	local output_script="$3"
	local metadata_file="$4"
	local helper_script="$HELPERS_DIR/tximport_helper.R"
	
	if [[ -f "$helper_script" ]]; then
		cp "$helper_script" "$output_script"
		chmod +x "$output_script"
		log_info "[TXIMPORT] Copied helper to: $output_script"
	else
		log_error "tximport_helper.R not found: $helper_script"
		return 1
	fi
}

# ==============================================================================
# GENE-TRANSCRIPT MAPPING
# ==============================================================================

# Create gene-transcript mapping from FASTA
create_gene_trans_map() {
	local fasta="$1"
	local output_file="${2:-${fasta}.gene_trans_map}"
	
	if [[ -f "$output_file" ]]; then
		log_info "Gene-transcript map already exists: $output_file"
		return 0
	fi
	
	log_info "Creating gene-transcript mapping from FASTA..."
	
	# Detect if this is a Trinity assembly
	local is_trinity=false
	if grep -q "^>TRINITY_" "$fasta" 2>/dev/null; then
		is_trinity=true
		log_info "Detected Trinity assembly format"
	fi
	
	if [[ "$is_trinity" == "true" ]]; then
		awk '/^>/ {
			trans = $1
			gsub(/^>/, "", trans)
			gene = trans
			if (match(gene, /^(.+)_i[0-9]+$/, arr)) {
				gene = arr[1]
			}
			print gene "\t" trans
		}' "$fasta" > "$output_file"
	else
		grep "^>" "$fasta" | sed 's/^>//' | awk '{
			trans=$1
			if (match($0, /gene=([^ ]+)/, arr)) {
				gene=arr[1]
			} else if (match($0, /gene_id[=:]([^ ]+)/, arr)) {
				gene=arr[1]
			} else if (match(trans, /^([^|]+)\|/, arr)) {
				gene=arr[1]
			} else {
				gene=trans
			}
			print gene "\t" trans
		}' > "$output_file"
	fi
	
	local unique_genes=$(cut -f1 "$output_file" | sort -u | wc -l)
	local total_transcripts=$(wc -l < "$output_file")
	log_info "Created gene-transcript map: $unique_genes genes, $total_transcripts transcripts"
}

# ==============================================================================
# NORMALIZATION FUNCTIONS
# ==============================================================================

normalize_expression_data() {
	local matrix_dir="$1"
	local method="$2"
	
	if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
		log_step "Applying TMM normalization for $method"
		if command -v normalize_matrix.pl >/dev/null 2>&1; then
			run_with_space_time_log normalize_matrix.pl "$matrix_dir/genes.counts.matrix" \
				--est_method "$method" \
				--out_prefix "$matrix_dir/genes.TMM" || \
				log_warn "TMM normalization failed for $method"
		else
			log_warn "normalize_matrix.pl not found. Skipping TMM normalization."
		fi
	fi
}

# ==============================================================================
# CROSS-METHOD VALIDATION
# ==============================================================================

compare_methods_summary() {
	local fasta_tag="$1"
	
	log_step "Cross-Method Validation Summary for $fasta_tag"
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	
	# Check each method's output
	local methods=(
		"Method 1 (HISAT2 Ref-Guided):$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/deseq2_input/gene_count_matrix.csv"
		"Method 3 (STAR):$STAR_ALIGN_ROOT/$fasta_tag/6_matrices_from_stringtie/gene_counts_tximport.tsv"
		"Method 4 (Salmon SAF):$SALMON_SAF_MATRIX_ROOT/$fasta_tag/deseq2_input/gene_count_matrix.csv"
		"Method 5 (Bowtie2+RSEM):$RSEM_MATRIX_ROOT/$fasta_tag/deseq2_input/gene_count_matrix.csv"
	)
	
	for method_info in "${methods[@]}"; do
		local method_name="${method_info%%:*}"
		local matrix_file="${method_info#*:}"
		
		if [[ -f "$matrix_file" ]]; then
			local genes=$(tail -n +2 "$matrix_file" | wc -l)
			local delim=","
			[[ "$matrix_file" == *.tsv ]] && delim=$'\t'
			local samples=$(head -n1 "$matrix_file" | tr "$delim" '\n' | tail -n +2 | wc -l)
			log_info "✅ $method_name: $genes genes, $samples samples"
		else
			log_info "❌ $method_name: Not found"
		fi
	done
	
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
}
