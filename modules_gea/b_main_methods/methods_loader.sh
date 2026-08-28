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
METHODS_MASTER_SOURCED="true"

# Source dependencies
# Use exported MODULES_DIR to avoid cd+dirname+pwd subshell fork; fallback for standalone sourcing
SCRIPT_DIR="${MODULES_DIR:+${MODULES_DIR}/b_main_methods}"
if [[ -z "$SCRIPT_DIR" ]]; then
	SCRIPT_DIR="${BASH_SOURCE[0]%/*}"; [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
	SCRIPT_DIR="$(cd "$SCRIPT_DIR" 2>/dev/null && pwd)"
fi

# Source all method scripts (log_error is available via logging_utils.sh sourced earlier in the chain)
source "$SCRIPT_DIR/m1_hisat2_ref_guided.sh" || { log_error "Failed to source m1_hisat2_ref_guided.sh"; return 1; }
source "$SCRIPT_DIR/m2_hisat2_de_novo.sh"    || { log_error "Failed to source m2_hisat2_de_novo.sh"; return 1; }
source "$SCRIPT_DIR/m3_star_alignment.sh"    || { log_error "Failed to source m3_star_alignment.sh"; return 1; }
source "$SCRIPT_DIR/m4_salmon_saf.sh"        || { log_error "Failed to source m4_salmon_saf.sh"; return 1; }
source "$SCRIPT_DIR/m5_bowtie2_rsem.sh"      || { log_error "Failed to source m5_bowtie2_rsem.sh"; return 1; }

# ==============================================================================
# METHOD COMPARISON AND VALIDATION FUNCTIONS
# ==============================================================================

# Cross-method validation and comparison
compare_methods_summary() {
	local fasta_tag="$1"
	
	log_step "Cross-Method Validation Summary for $fasta_tag"
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	
	# Method 1: HISAT2 Reference-Guided
	local m1_matrix="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/deseq2_input/gene_count_matrix.csv"
	if [[ -f "$m1_matrix" ]]; then
		local m1_genes m1_samples
		read -r m1_genes m1_samples < <(awk -F',' 'NR==1{s=NF-1} END{printf "%d %d", NR-1, s}' "$m1_matrix")
		log_info "✅ Method 1 (HISAT2 Ref-Guided): $m1_genes genes, $m1_samples samples"
		log_info "   Status: BEST for publication (true read counts via prepDE.py)"
	fi
	
	# Method 2: HISAT2 De Novo
	local m2_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT"
	if [[ -d "$m2_dir" ]]; then
		log_info "✅ Method 2 (HISAT2 De Novo): Transcript discovery mode"
		log_info "   Status: De novo assembly without reference GTF"
	fi
	
	# Method 3: STAR Alignment
	# Output filename follows the tximport naming convention set by tximport_star_helper.R
	local m3_matrix="$STAR_MATRIX_ROOT/gene_level/${fasta_tag}_NumReads_Gene_ID_from_${fasta_tag}_gene_level.tsv"
	if [[ -f "$m3_matrix" ]]; then
		local m3_genes m3_samples
		read -r m3_genes m3_samples < <(awk -F'\t' 'NR==1{s=NF-1} END{printf "%d %d", NR-1, s}' "$m3_matrix")
		log_info "✅ Method 3 (STAR): $m3_genes genes, $m3_samples samples"
		log_info "   Status: Splice-aware alignment with Salmon quantification"
	fi
	
	# Method 4: Salmon SAF
	local m4_matrix="$SALMON_MATRIX_ROOT/deseq2_input/gene_count_matrix.csv"
	if [[ -f "$m4_matrix" ]]; then
		log_info "✅ Method 4 (Salmon SAF): Fast pseudo-alignment"
		log_info "   Status: EXCELLENT for publication (fast, accurate, modern)"
	fi
	
	# Method 5: Bowtie2 + RSEM
	local m5_matrix="$RSEM_MATRIX_ROOT/deseq2_input/gene_count_matrix.csv"
	if [[ -f "$m5_matrix" ]]; then
		log_info "⚠️  Method 5 (Bowtie2+RSEM): RSEM expected counts"
		log_info "   Status: Uses RSEM expected_count (statistical estimates)"
	fi
	
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
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
#   I_INPUTS/inputs/eggplant/sample_conditions.txt
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
log_info "[METHODS] Sample conditions file: ${SAMPLE_CONDITIONS_FILE:-I_INPUTS/inputs/eggplant/sample_conditions.txt}"
