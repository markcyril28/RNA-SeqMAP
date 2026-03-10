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
# IMPORTANT PARAMETERS (tweak here)
# ==============================================================================

# Total CPU threads available to the pipeline
THREADS="${THREADS:-12}"
# Base parallel job count; PARALLEL_JOBS inherits this if not set separately
JOBS="${JOBS:-2}"

# Concurrent sample jobs for GNU Parallel; each job gets THREADS/PARALLEL_JOBS threads
PARALLEL_JOBS="${PARALLEL_JOBS:-${JOBS:-2}}"

# BAM retention: "y" = keep BAM files after processing, "n" = delete to save disk space
keep_bam_global="${keep_bam_global:-n}"

# ==============================================================================
# BOWTIE2/RSEM CONFIGURATION
# ==============================================================================

# Alignment sensitivity for rsem-calculate-expression --bowtie2-sensitivity-level
# Valid options: sensitive (default), very-sensitive, fast, very-fast
BOWTIE2_MODE="${BOWTIE2_MODE:-sensitive}"

# ==============================================================================
# STAR CONFIGURATION
# ==============================================================================

# Genome loading: NoSharedMemory (safe default), LoadAndKeep (faster multi-run on HPC)
STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-NoSharedMemory}"
# Used to compute sjdbOverhang = read_length - 1; set to actual sequencing read length
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"
# Strandedness: None (unstranded), Forward, Reverse
STAR_STRAND_SPECIFIC="${STAR_STRAND_SPECIFIC:-None}"

# ==============================================================================
# OVERWRITE MODE
# ==============================================================================
# OVERWRITE_MODE controls whether alignment steps skip existing outputs.
# "overwrite" = re-run even if outputs exist; "skip" = skip existing (default).
# Derived from OVERWRITE_EXISTING (set/exported by run_all_post_processing.sh).
# Respects an existing OVERWRITE_MODE set by a_GEA_script_v12.sh or the environment.
if [[ "${OVERWRITE_EXISTING:-FALSE}" == "TRUE" ]]; then
	OVERWRITE_MODE="overwrite"
elif [[ -z "${OVERWRITE_MODE:-}" ]]; then
	OVERWRITE_MODE="skip"
fi
export OVERWRITE_MODE

# ==============================================================================
# POST PROCESSING ROOT
# ==============================================================================
# Convert to absolute path if relative (prevents STAR/tool output file errors)
_POST_PROC_DEFAULT="${POST_PROCESSING_ROOT:-3_POST_PROC}"
if [[ "$_POST_PROC_DEFAULT" != /* ]]; then
	POST_PROCESSING_ROOT="$(pwd)/$_POST_PROC_DEFAULT"
else
	POST_PROCESSING_ROOT="$_POST_PROC_DEFAULT"
fi
unset _POST_PROC_DEFAULT

# ==============================================================================
# ALIGNMENT RESULTS ROOT
# ==============================================================================
# Alignment outputs live separately from post-processing
_ALIGN_DEFAULT="${ALIGNMENT_RESULTS_ROOT:-2_ALIGNMENT_RESULTs}"
if [[ "$_ALIGN_DEFAULT" != /* ]]; then
	ALIGNMENT_RESULTS_ROOT="$(pwd)/$_ALIGN_DEFAULT"
else
	ALIGNMENT_RESULTS_ROOT="$_ALIGN_DEFAULT"
fi
unset _ALIGN_DEFAULT

# ==============================================================================
# SAMPLE METADATA CONFIGURATION
# ==============================================================================
# Path to the sample conditions file (tab-separated: SRR_ID condition batch)
# Convert to absolute path if relative (same pattern as POST_PROCESSING_ROOT)
_SAMPLE_COND_DEFAULT="${SAMPLE_CONDITIONS_FILE:-inputs/sample_conditions.txt}"
if [[ "$_SAMPLE_COND_DEFAULT" != /* ]]; then
	SAMPLE_CONDITIONS_FILE="$(pwd)/$_SAMPLE_COND_DEFAULT"
else
	SAMPLE_CONDITIONS_FILE="$_SAMPLE_COND_DEFAULT"
fi
unset _SAMPLE_COND_DEFAULT

# ==============================================================================
# METHOD 1: HISAT2 REFERENCE GUIDED DIRECTORIES
# ==============================================================================
HISAT2_REF_GUIDED_ROOT="$ALIGNMENT_RESULTS_ROOT/M1_HISAT2_RefGuided/HISAT2_WD"
HISAT2_REF_GUIDED_INDEX_DIR="$HISAT2_REF_GUIDED_ROOT/index"
STRINGTIE_HISAT2_REF_GUIDED_ROOT="$ALIGNMENT_RESULTS_ROOT/M1_HISAT2_RefGuided/stringtie_WD"
HISAT2_REF_GUIDED_MATRIX_ROOT="$POST_PROCESSING_ROOT/M1_HISAT2_RefGuided/count_matrices_from_stringtie"

# ==============================================================================
# METHOD 2: HISAT2 DE NOVO DIRECTORIES
# ==============================================================================
HISAT2_DE_NOVO_ROOT="$ALIGNMENT_RESULTS_ROOT/M2_HISAT2_DeNovo/HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="$HISAT2_DE_NOVO_ROOT/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="$ALIGNMENT_RESULTS_ROOT/M2_HISAT2_DeNovo/stringtie_WD"
HISAT2_DE_NOVO_MATRIX_ROOT="$POST_PROCESSING_ROOT/M2_HISAT2_DeNovo/count_matrices_from_stringtie"

# ==============================================================================
# METHOD 3: STAR ALIGNMENT DIRECTORIES
# ==============================================================================
STAR_ALIGN_ROOT="$ALIGNMENT_RESULTS_ROOT/M3_STAR_Align"
STAR_INDEX_ROOT="$STAR_ALIGN_ROOT/STAR_index"
STAR_MATRIX_ROOT="$POST_PROCESSING_ROOT/M3_STAR_Align/count_matrices_from_STAR"

# ==============================================================================
# METHOD 4: SALMON SAF DIRECTORIES
# ==============================================================================
SALMON_SAF_ROOT="$ALIGNMENT_RESULTS_ROOT/M4_Salmon_Saf"
SALMON_INDEX_ROOT="$SALMON_SAF_ROOT/Salmon_WD/index"
SALMON_QUANT_ROOT="$SALMON_SAF_ROOT/Salmon_Quant"
SALMON_SAF_MATRIX_ROOT="$POST_PROCESSING_ROOT/M4_Salmon_Saf/count_matrices_from_Salmon_Quant"
SALMON_MATRIX_ROOT="${SALMON_MATRIX_ROOT:-$SALMON_SAF_MATRIX_ROOT}"

# ==============================================================================
# METHOD 5: BOWTIE2 + RSEM DIRECTORIES
# ==============================================================================
BOWTIE2_RSEM_ROOT="$ALIGNMENT_RESULTS_ROOT/M5_RSEM_Bowtie2"
RSEM_INDEX_ROOT="$BOWTIE2_RSEM_ROOT/Bowtie2_WD/index"
RSEM_QUANT_ROOT="$BOWTIE2_RSEM_ROOT/RSEM_Quant_WD"
RSEM_MATRIX_ROOT="$POST_PROCESSING_ROOT/M5_RSEM_Bowtie2/count_matrices_from_RSEM_Quant"

# ==============================================================================
# SRR SAMPLE ARRAYS (initialize if not set)
# ==============================================================================
if ! declare -p SRR_COMBINED_LIST &>/dev/null; then
	declare -a SRR_COMBINED_LIST=()
fi

# ==============================================================================
# INITIALIZE ALL METHOD DIRECTORIES
# ==============================================================================
init_method_directories() {
	# Alignment directories
	mkdir -p "$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
		"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
		"$STAR_ALIGN_ROOT" "$STAR_INDEX_ROOT" \
		"$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" \
		"$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT"
	# Post-processing matrix directories
	mkdir -p "$HISAT2_REF_GUIDED_MATRIX_ROOT" "$HISAT2_DE_NOVO_MATRIX_ROOT" \
		"$STAR_MATRIX_ROOT" "$SALMON_SAF_MATRIX_ROOT" "$RSEM_MATRIX_ROOT"
}

# ==============================================================================
# DYNAMIC FASTA-BASED OUTPUT DIRECTORIES
# ==============================================================================
# Call before running any method to isolate outputs by input FASTA filename.
# Reconfigures all method directory variables to include fasta_tag as a
# subdirectory under ALIGNMENT_RESULTS_ROOT and POST_PROCESSING_ROOT.
# ==============================================================================

set_fasta_output_dirs() {
	local fasta_tag="$1"
	[[ -z "$fasta_tag" ]] && { echo "ERROR: set_fasta_output_dirs requires a fasta_tag argument" >&2; return 1; }

	# Method 1: HISAT2 Reference Guided
	HISAT2_REF_GUIDED_ROOT="$ALIGNMENT_RESULTS_ROOT/M1_HISAT2_RefGuided/HISAT2_WD/$fasta_tag"
	HISAT2_REF_GUIDED_INDEX_DIR="$HISAT2_REF_GUIDED_ROOT/index"
	STRINGTIE_HISAT2_REF_GUIDED_ROOT="$ALIGNMENT_RESULTS_ROOT/M1_HISAT2_RefGuided/stringtie_WD/$fasta_tag"
	HISAT2_REF_GUIDED_MATRIX_ROOT="$POST_PROCESSING_ROOT/M1_HISAT2_RefGuided/count_matrices_from_stringtie/$fasta_tag"

	# Method 2: HISAT2 De Novo
	HISAT2_DE_NOVO_ROOT="$ALIGNMENT_RESULTS_ROOT/M2_HISAT2_DeNovo/HISAT2_WD/$fasta_tag"
	HISAT2_DE_NOVO_INDEX_DIR="$HISAT2_DE_NOVO_ROOT/index"
	STRINGTIE_HISAT2_DE_NOVO_ROOT="$ALIGNMENT_RESULTS_ROOT/M2_HISAT2_DeNovo/stringtie_WD/$fasta_tag"
	HISAT2_DE_NOVO_MATRIX_ROOT="$POST_PROCESSING_ROOT/M2_HISAT2_DeNovo/count_matrices_from_stringtie/$fasta_tag"

	# Method 3: STAR Alignment
	STAR_ALIGN_ROOT="$ALIGNMENT_RESULTS_ROOT/M3_STAR_Align"
	STAR_INDEX_ROOT="$ALIGNMENT_RESULTS_ROOT/M3_STAR_Align/STAR_index/$fasta_tag"
	STAR_MATRIX_ROOT="$POST_PROCESSING_ROOT/M3_STAR_Align/count_matrices_from_STAR/$fasta_tag"

	# Method 4: Salmon SAF
	SALMON_SAF_ROOT="$ALIGNMENT_RESULTS_ROOT/M4_Salmon_Saf"
	SALMON_INDEX_ROOT="$ALIGNMENT_RESULTS_ROOT/M4_Salmon_Saf/Salmon_WD/$fasta_tag/index"
	SALMON_QUANT_ROOT="$ALIGNMENT_RESULTS_ROOT/M4_Salmon_Saf/Salmon_Quant/$fasta_tag"
	SALMON_SAF_MATRIX_ROOT="$POST_PROCESSING_ROOT/M4_Salmon_Saf/count_matrices_from_Salmon_Quant/$fasta_tag"
	SALMON_MATRIX_ROOT="$SALMON_SAF_MATRIX_ROOT"

	# Method 5: Bowtie2 + RSEM
	BOWTIE2_RSEM_ROOT="$ALIGNMENT_RESULTS_ROOT/M5_RSEM_Bowtie2"
	RSEM_INDEX_ROOT="$ALIGNMENT_RESULTS_ROOT/M5_RSEM_Bowtie2/Bowtie2_WD/$fasta_tag/index"
	RSEM_QUANT_ROOT="$ALIGNMENT_RESULTS_ROOT/M5_RSEM_Bowtie2/RSEM_Quant_WD/$fasta_tag"
	RSEM_MATRIX_ROOT="$POST_PROCESSING_ROOT/M5_RSEM_Bowtie2/count_matrices_from_RSEM_Quant/$fasta_tag"

	log_info "[CONFIG] Output directories configured for FASTA: $fasta_tag"
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
