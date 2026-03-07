#!/bin/bash
# ==============================================================================
# MAIN METHODS - MASTER LOADER
# ==============================================================================
# Sources all GEA analysis method scripts from b_main_methods/
# This file provides backward compatibility with older scripts
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${METHODS_MASTER_SOURCED:-}" == "true" ]] && return 0
export METHODS_MASTER_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source all method scripts
source "$SCRIPT_DIR/m1_hisat2_ref_guided.sh"
source "$SCRIPT_DIR/m2_hisat2_de_novo.sh"
source "$SCRIPT_DIR/m3_star_alignment.sh"
source "$SCRIPT_DIR/m4_salmon_saf.sh"
source "$SCRIPT_DIR/m5_bowtie2_rsem.sh"

# ==============================================================================
# METHOD COMPARISON AND VALIDATION FUNCTIONS
# ==============================================================================

# Cross-method validation and comparison
compare_methods_summary() {
	local fasta_tag="$1"
	
	log_step "Cross-Method Validation Summary for $fasta_tag"
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	
	# Method 1: HISAT2 Reference-Guided
	local m1_matrix="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/deseq2_input/gene_count_matrix.csv"
	if [[ -f "$m1_matrix" ]]; then
		local m1_genes=$(tail -n +2 "$m1_matrix" | wc -l)
		local m1_samples=$(head -n1 "$m1_matrix" | tr ',' '\n' | tail -n +2 | wc -l)
		log_info "✅ Method 1 (HISAT2 Ref-Guided): $m1_genes genes, $m1_samples samples"
		log_info "   Status: BEST for publication (true read counts via prepDE.py)"
	fi
	
	# Method 2: HISAT2 De Novo
	local m2_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$fasta_tag"
	if [[ -d "$m2_dir" ]]; then
		log_info "✅ Method 2 (HISAT2 De Novo): Transcript discovery mode"
		log_info "   Status: De novo assembly without reference GTF"
	fi
	
	# Method 3: STAR Alignment
	local m3_matrix="$STAR_ALIGN_ROOT/$fasta_tag/6_matrices_from_stringtie/gene_counts_tximport.tsv"
	if [[ -f "$m3_matrix" ]]; then
		local m3_genes=$(tail -n +2 "$m3_matrix" | wc -l)
		local m3_samples=$(head -n1 "$m3_matrix" | tr '\t' '\n' | tail -n +2 | wc -l)
		log_info "✅ Method 3 (STAR): $m3_genes genes, $m3_samples samples"
		log_info "   Status: Splice-aware alignment with Salmon quantification"
	fi
	
	# Method 4: Salmon SAF
	local m4_matrix="$SALMON_MATRIX_ROOT/*/deseq2_input/gene_count_matrix.csv"
	if compgen -G "$m4_matrix" >/dev/null 2>&1; then
		log_info "✅ Method 4 (Salmon SAF): Fast pseudo-alignment"
		log_info "   Status: EXCELLENT for publication (fast, accurate, modern)"
	fi
	
	# Method 5: Bowtie2 + RSEM
	local m5_matrix="$RSEM_MATRIX_ROOT/*/deseq2_input/gene_count_matrix.csv"
	if compgen -G "$m5_matrix" >/dev/null 2>&1; then
		log_info "⚠️  Method 5 (Bowtie2+RSEM): RSEM expected counts"
		log_info "   Status: Uses RSEM expected_count (statistical estimates)"
	fi
	
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
}

# Apply normalization to expression matrices
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
# IMPORTANT: Sample Metadata Configuration
# ==============================================================================
# All pipeline methods create a default sample_metadata.csv file with:
#   - All samples assigned to "treatment" condition
#   - All samples in batch "1"
# 
# THIS MUST BE CUSTOMIZED for differential expression analysis!
# 
# Option 1: Edit the sample_conditions.txt file at:
#   0_INPUT_FASTAs/sample_conditions.txt
#   Format (tab-separated):
#     SRR_ID	condition	batch
#     SRR3884597	Flowers	1
#     SRR3884653	Fruits	1
#
# Option 2: Set SAMPLE_CONDITIONS_FILE environment variable to custom path
#
# Option 3: Manually edit the generated sample_metadata.csv files
# ==============================================================================

log_info "[METHODS] All GEA analysis methods loaded successfully"
log_info "[METHODS] Sample conditions file: ${SAMPLE_CONDITIONS_FILE:-0_INPUT_FASTAs/sample_conditions.txt}"
