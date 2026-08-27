#!/bin/bash

# ==============================================================================
# IMPORTANT PARAMETERS
# ==============================================================================

# Runtime Configuration
THREADS=64                              # Threads for parallel operations
JOBS=2									# Parallel jobs for GNU Parallel
USE_GNU_PARALLEL="TRUE"                 # TRUE/FALSE for GNU Parallel
keep_bam_global="n"                     # y=keep BAM files, n=delete after

# Pipeline Stages (comment/uncomment to enable/disable)
PIPELINE_STAGES=(
	#"MAMBA_INSTALLATION"

	# Option A: Separate download and trim (keeps raw files)
	#"DOWNLOAD_SRR"
	#"TRIM_SRR"

	# Option B: Combined download+trim+cleanup (auto-deletes raw after trim)
	#"DOWNLOAD_TRIM_and_DELETE_RAW_SRR"

	#"GZIP_TRIMMED_FILES"
	#"QUALITY_CONTROL"

	#"DELETE_RAW_SRR"				# Manually delete raw SRR files
	#"DELETE_TRIMMED_FASTQ_FILES"	# Manually delete trimmed files

	#"METHOD_1_HISAT2_REF_GUIDED"
	"METHOD_2_HISAT2_DE_NOVO"
	#"METHOD_3_STAR_ALIGNMENT"
	"METHOD_4_SALMON_SAF"
	"METHOD_5_BOWTIE2_RSEM"
)

# ==============================================================================
# SHARED RUNTIME (conda, modules, thread calculation)
# ==============================================================================

source "I_INPUTS/config/eggplant/shared/runtime_defaults.sh"

# ==============================================================================
# INPUT FILES AND DATA SOURCES
# ==============================================================================

decoy="inputs/fasta/experimental/TEST.fasta"

# FASTA Files for Analysis (transcripts, required for M2, M4, M5)
ALL_FASTA_FILES=(
	"inputs/fasta/reference_genome/GPE001970_transcripts.fa"
	"inputs/fasta/reference_genome/Eggplant_V4.1_transcripts.function.fa"
)

# ==============================================================================
# RNA-SEQ DATA SOURCES (SRR LISTS) — sourced from shared config
# ==============================================================================

source "config/shared/srr_datasets.sh"

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

# POST_PROCESSING_ROOT is set and converted to absolute by global_config_method.sh
# (sourced via modules_loader.sh above). Do NOT override it here — that would
# revert the absolute path back to a relative one, creating inconsistency with
# ALIGNMENT_RESULTS_ROOT which remains absolute.

# Create required directories
# NOTE: SALMON_INDEX_ROOT / SALMON_QUANT_ROOT / SALMON_SAF_MATRIX_ROOT and the RSEM equivalents
# are NOT pre-created here because they include the fasta_tag subdirectory (set by
# set_fasta_output_dirs() inside run_all()). Creating them now would produce stale base-level
# directories. Each method pipeline creates its own directories at the correct paths.
mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
	"$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
	"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
	"$STAR_ALIGN_ROOT" "$STAR_INDEX_ROOT" \
	"$SALMON_SAF_ROOT" \
	"$BOWTIE2_RSEM_ROOT" \
	"$HISAT2_REF_GUIDED_MATRIX_ROOT" "$HISAT2_DE_NOVO_MATRIX_ROOT" "$STAR_MATRIX_ROOT"

# ==============================================================================
# CLEANUP OPTIONS AND TESTING ESSENTIALS
# ==============================================================================
# Uncomment lines below to remove previous results before re-running.
# WARNING: These are destructive operations — verify before uncommenting.
# ==============================================================================

ACTIVATE_RM=FALSE

# --- Method 2: HISAT2 De Novo ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_ROOT"                   # HISAT2 de novo alignments
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_INDEX_DIR"              # HISAT2 de novo index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STRINGTIE_HISAT2_DE_NOVO_ROOT"         # StringTie (de novo) assemblies
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_MATRIX_ROOT"            # Count matrices (de novo)

# --- Method 4: Salmon SAF ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_SAF_ROOT"                       # Salmon SAF root
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_INDEX_ROOT"                     # Salmon index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_QUANT_ROOT"                     # Salmon quantification
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_SAF_MATRIX_ROOT"                # Count matrices (Salmon)

# --- Method 5: Bowtie2 + RSEM ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$BOWTIE2_RSEM_ROOT"                     # Bowtie2/RSEM root
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RSEM_INDEX_ROOT"                       # RSEM index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RSEM_QUANT_ROOT"                       # RSEM quantification
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RSEM_MATRIX_ROOT"                      # Count matrices (RSEM)
