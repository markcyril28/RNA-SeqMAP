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

# Total CPU threads available to the pipeline (auto-detect if not set)
THREADS="${THREADS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 12)}"
# Base parallel job count; set to "auto" to calculate from THREADS / OPTIMAL_THREADS_PER_JOB
JOBS="${JOBS:-2}"
# Default optimal threads per job for Stage 2 (alignment programs scale well up to ~16)
OPTIMAL_THREADS_PER_JOB="${OPTIMAL_THREADS_PER_JOB:-16}"
# Track original setting so method scripts can re-resolve with program-specific optimal
_JOBS_MODE="${JOBS}"
# Resolve JOBS="auto" → numeric value
if [[ "${JOBS}" == "auto" || "${JOBS}" == "AUTO" ]]; then
	JOBS=$(( THREADS / OPTIMAL_THREADS_PER_JOB ))
	(( JOBS < 1 )) && JOBS=1
fi

# Concurrent sample jobs for GNU Parallel; each job gets THREADS/PARALLEL_JOBS threads
PARALLEL_JOBS="${PARALLEL_JOBS:-${JOBS:-2}}"

# BAM retention: "y" = keep BAM files after processing, "n" = delete to save disk space
keep_bam_global="${keep_bam_global:-n}"

# ==============================================================================
# BOWTIE2/RSEM CONFIGURATION
# ==============================================================================

# Alignment sensitivity for rsem-calculate-expression --bowtie2-sensitivity-level
# Valid options: sensitive (default), very_sensitive, fast, very_fast
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
# Derived from OVERWRITE_EXISTING (set/exported by run_post_processing.sh).
# Respects an existing OVERWRITE_MODE set by a_GEA_script_v12.sh or the environment.
_ow="${OVERWRITE_EXISTING:-FALSE}"
if [[ "${_ow^^}" == "TRUE" ]]; then
	OVERWRITE_MODE="overwrite"
elif [[ -z "${OVERWRITE_MODE:-}" ]]; then
	OVERWRITE_MODE="skip"
fi
export OVERWRITE_MODE
unset _ow

# ==============================================================================
# PATH HELPER — O(1) absolute path resolution, no subprocess spawns
# ==============================================================================
# Converts relative paths to absolute using $PWD. Consolidates the repeated
# if [[ "$path" != /* ]] pattern (was 3 copies; now single source of truth).
_make_absolute_path() {
	local p="$1"
	[[ "$p" != /* ]] && p="$PWD/$p"
	printf '%s' "$p"
}

# ==============================================================================
# POST PROCESSING ROOT
# ==============================================================================
POST_PROCESSING_ROOT="$(_make_absolute_path "${POST_PROCESSING_ROOT:-3_POST_PROC}")"

# ==============================================================================
# ALIGNMENT RESULTS ROOT
# ==============================================================================
ALIGNMENT_RESULTS_ROOT="$(_make_absolute_path "${ALIGNMENT_RESULTS_ROOT:-2_ALIGNMENT_RESULTs}")"

# ==============================================================================
# SAMPLE METADATA CONFIGURATION
# ==============================================================================
# Path to the sample conditions file (tab-separated: SRR_ID condition batch)
SAMPLE_CONDITIONS_FILE="$(_make_absolute_path "${SAMPLE_CONDITIONS_FILE:-inputs/sample_conditions.txt}")"

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
	[[ -z "$fasta_tag" ]] && { log_error "set_fasta_output_dirs requires a fasta_tag argument"; return 1; }

	# Method 1: HISAT2 Reference Guided
	HISAT2_REF_GUIDED_ROOT="$ALIGNMENT_RESULTS_ROOT/M1_HISAT2_RefGuided/HISAT2_WD/$fasta_tag"
	HISAT2_REF_GUIDED_INDEX_DIR="$HISAT2_REF_GUIDED_ROOT/index"
	STRINGTIE_HISAT2_REF_GUIDED_ROOT="$ALIGNMENT_RESULTS_ROOT/M1_HISAT2_RefGuided/stringtie_WD/$fasta_tag"
	# Matrix root does NOT include fasta_tag — stringtie_matrix_builder.sh writes
	# gene-group folders directly under count_matrices_from_stringtie/
	HISAT2_REF_GUIDED_MATRIX_ROOT="$POST_PROCESSING_ROOT/M1_HISAT2_RefGuided/count_matrices_from_stringtie"

	# Method 2: HISAT2 De Novo
	HISAT2_DE_NOVO_ROOT="$ALIGNMENT_RESULTS_ROOT/M2_HISAT2_DeNovo/HISAT2_WD/$fasta_tag"
	HISAT2_DE_NOVO_INDEX_DIR="$HISAT2_DE_NOVO_ROOT/index"
	STRINGTIE_HISAT2_DE_NOVO_ROOT="$ALIGNMENT_RESULTS_ROOT/M2_HISAT2_DeNovo/stringtie_WD/$fasta_tag"
	# Matrix root does NOT include fasta_tag — stringtie_matrix_builder.sh writes
	# gene-group folders directly under count_matrices_from_stringtie/
	HISAT2_DE_NOVO_MATRIX_ROOT="$POST_PROCESSING_ROOT/M2_HISAT2_DeNovo/count_matrices_from_stringtie"

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

