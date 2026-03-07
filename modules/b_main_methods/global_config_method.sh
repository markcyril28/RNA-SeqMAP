#!/bin/bash
# ==============================================================================
# MAIN METHODS GLOBAL CONFIGURATION
# ==============================================================================
# Centralized configuration for all GEA analysis methods
# Sourced by: All method scripts in b_main_methods/
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${METHOD_CONFIG_SOURCED:-}" == "true" ]] && return 0
export METHOD_CONFIG_SOURCED="true"

# ==============================================================================
# IMPORTANT PARAMETERS - MODIFY THESE AT TOP
# ==============================================================================

# Thread Configuration
THREADS="${THREADS:-12}"
JOBS="${JOBS:-2}"

# BAM File Retention
keep_bam_global="${keep_bam_global:-n}"  # y = keep BAM files, n = delete after processing

# ==============================================================================
# BOWTIE2/RSEM CONFIGURATION
# ==============================================================================

# Bowtie2 alignment mode (options: very-sensitive-local, sensitive-local, fast-local, very-fast-local)
BOWTIE2_MODE="${BOWTIE2_MODE:-very-sensitive-local}"

# ==============================================================================
# STAR CONFIGURATION
# ==============================================================================

STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-NoSharedMemory}"
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"
STAR_STRAND_SPECIFIC="${STAR_STRAND_SPECIFIC:-intronMotif}"

# ==============================================================================
# POST PROCESSING ROOT
# ==============================================================================
# Convert to absolute path if relative (prevents STAR/tool output file errors)
_POST_PROC_DEFAULT="${POST_PROCESSING_ROOT:-4_POST_PROC}"
if [[ "$_POST_PROC_DEFAULT" != /* ]]; then
	POST_PROCESSING_ROOT="$(pwd)/$_POST_PROC_DEFAULT"
else
	POST_PROCESSING_ROOT="$_POST_PROC_DEFAULT"
fi
unset _POST_PROC_DEFAULT

# ==============================================================================
# SAMPLE METADATA CONFIGURATION
# ==============================================================================
# Path to the sample conditions file (tab-separated: SRR_ID condition batch)
# Convert to absolute path if relative (same pattern as POST_PROCESSING_ROOT)
_SAMPLE_COND_DEFAULT="${SAMPLE_CONDITIONS_FILE:-0_INPUT_FASTAs/sample_conditions.txt}"
if [[ "$_SAMPLE_COND_DEFAULT" != /* ]]; then
	SAMPLE_CONDITIONS_FILE="$(pwd)/$_SAMPLE_COND_DEFAULT"
else
	SAMPLE_CONDITIONS_FILE="$_SAMPLE_COND_DEFAULT"
fi
unset _SAMPLE_COND_DEFAULT

# ==============================================================================
# METHOD 1: HISAT2 REFERENCE GUIDED DIRECTORIES
# ==============================================================================
HISAT2_REF_GUIDED_ROOT="$POST_PROCESSING_ROOT/M1_HISAT2_RefGuided/4_HISAT2_WD"
HISAT2_REF_GUIDED_INDEX_DIR="$HISAT2_REF_GUIDED_ROOT/index"
STRINGTIE_HISAT2_REF_GUIDED_ROOT="$POST_PROCESSING_ROOT/M1_HISAT2_RefGuided/5_stringtie_WD"
HISAT2_REF_GUIDED_MATRIX_ROOT="$POST_PROCESSING_ROOT/M1_HISAT2_RefGuided/6_matrices"

# ==============================================================================
# METHOD 2: HISAT2 DE NOVO DIRECTORIES
# ==============================================================================
HISAT2_DE_NOVO_ROOT="$POST_PROCESSING_ROOT/M2_HISAT2_DeNovo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="$HISAT2_DE_NOVO_ROOT/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="$POST_PROCESSING_ROOT/M2_HISAT2_DeNovo/5_stringtie_WD"
HISAT2_DE_NOVO_MATRIX_ROOT="$POST_PROCESSING_ROOT/M2_HISAT2_DeNovo/6_matrices"

# ==============================================================================
# METHOD 3: STAR ALIGNMENT DIRECTORIES
# ==============================================================================
STAR_ALIGN_ROOT="$POST_PROCESSING_ROOT/M3_STAR_Align"
STAR_INDEX_ROOT="$STAR_ALIGN_ROOT/star_index"
STAR_GENOME_DIR="$STAR_ALIGN_ROOT/alignments"

# ==============================================================================
# METHOD 4: SALMON SAF DIRECTORIES
# ==============================================================================
SALMON_SAF_ROOT="$POST_PROCESSING_ROOT/M4_Salmon_Saf"
SALMON_INDEX_ROOT="$SALMON_SAF_ROOT/4_Salmon_WD/index"
SALMON_QUANT_ROOT="$SALMON_SAF_ROOT/5_Salmon_Quant_WD"
SALMON_SAF_MATRIX_ROOT="$SALMON_SAF_ROOT/6_matrices_from_Salmon"
SALMON_MATRIX_ROOT="${SALMON_MATRIX_ROOT:-$SALMON_SAF_MATRIX_ROOT}"

# ==============================================================================
# METHOD 5: BOWTIE2 + RSEM DIRECTORIES
# ==============================================================================
BOWTIE2_RSEM_ROOT="$POST_PROCESSING_ROOT/M5_RSEM_Bowtie2"
RSEM_INDEX_ROOT="$BOWTIE2_RSEM_ROOT/4_Bowtie2_WD/index"
RSEM_QUANT_ROOT="$BOWTIE2_RSEM_ROOT/5_RSEM_Quant_WD"
RSEM_MATRIX_ROOT="$BOWTIE2_RSEM_ROOT/6_matrices_from_RSEM"

# ==============================================================================
# SRR SAMPLE ARRAYS (initialize if not set)
# ==============================================================================
if ! declare -p SRR_COMBINED_LIST &>/dev/null 2>&1; then
	declare -a SRR_COMBINED_LIST=()
fi

# ==============================================================================
# INITIALIZE ALL METHOD DIRECTORIES
# ==============================================================================
init_method_directories() {
	mkdir -p "$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
		"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
		"$STAR_ALIGN_ROOT" "$STAR_INDEX_ROOT" "$STAR_GENOME_DIR" \
		"$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" "$SALMON_SAF_MATRIX_ROOT" \
		"$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT" "$RSEM_MATRIX_ROOT"
}

# ==============================================================================
# DISPLAY CONFIGURATION
# ==============================================================================
show_method_configuration() {
	log_info "=== METHOD CONFIGURATION ==="
	log_info "Threads: $THREADS"
	log_info "Keep BAM: $keep_bam_global"
	log_info "Bowtie2 Mode: $BOWTIE2_MODE"
	log_info "STAR Read Length: $STAR_READ_LENGTH"
	log_info "STAR Genome Load: $STAR_GENOME_LOAD"
	if type -t log_gpu_status &>/dev/null; then
		log_gpu_status
	fi
	log_info "============================="
}
