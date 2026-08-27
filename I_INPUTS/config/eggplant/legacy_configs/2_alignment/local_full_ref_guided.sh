#!/bin/bash

# ==============================================================================
# IMPORTANT PARAMETERS
# ==============================================================================

# Runtime Configuration
THREADS=64                              # Threads for parallel operations
JOBS=2									# Parallel jobs for GNU Parallel
USE_GNU_PARALLEL="TRUE"                 # TRUE/FALSE for GNU Parallel
keep_bam_global="n"                     # y=keep BAM files, n=delete after
STAR_READ_LENGTH=89                      # Actual read length for PRJNA328564 (89 bp)
export STAR_READ_LENGTH

# Strandness for HISAT2/StringTie (M1/M2)
# Options: "FR" (ligation/forward), "RF" (dUTP/TruSeq/reverse), "" (unstranded/auto-detect)
# Leave empty to enable auto-detection via infer_experiment.py on the first aligned sample.
HISAT2_STRANDNESS=""

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

	"METHOD_1_HISAT2_REF_GUIDED"
	#"METHOD_2_HISAT2_DE_NOVO"
	"METHOD_3_STAR_ALIGNMENT"
	#"METHOD_4_SALMON_SAF"
	#"METHOD_5_BOWTIE2_RSEM"
)

# ==============================================================================
# SHARED RUNTIME (conda, modules, thread calculation)
# ==============================================================================

source "I_INPUTS/config/eggplant/shared/runtime_defaults.sh"

# ==============================================================================
# INPUT FILES AND DATA SOURCES
# ==============================================================================

# ------------------------------------------------------------------------------
# GENOME REFERENCE PAIRS  (M1: HISAT2 Ref-Guided  |  M3: STAR)
# Format: "GTF_FILE|FASTA_FILE|STAR_TRANSCRIPTOME_FASTA"
#   Field 3 is optional — leave empty ("GTF|FASTA|") to use auto-detect.
# The pipeline loops through all uncommented pairs.
# ------------------------------------------------------------------------------
GENOME_REF_PAIRS=(
	"inputs/gtf/reference/GPE001970_genome.gtf|inputs/fasta/reference_genome/GPE001970_genome.fa|inputs/fasta/reference_genome/GPE001970_transcripts.fa"                                       	# GPE001970
	"inputs/gtf/reference/Eggplant_V4.1_function_IPR_final_stringtie.gtf|inputs/fasta/reference_genome/Eggplant_V4.1.fa|inputs/fasta/reference_genome/Eggplant_V4.1_transcripts.function.fa"  	# Eggplant V4.1
)

# ==============================================================================
# RNA-SEQ DATA SOURCES (SRR LISTS) — sourced from shared config
# ==============================================================================

source "config/shared/srr_datasets.sh"

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

POST_PROCESSING_ROOT="3_POST_PROC"
export POST_PROCESSING_ROOT

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

# --- Method 1: HISAT2 Reference-Guided ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_ROOT"                # HISAT2 ref-guided alignments
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_INDEX_DIR"           # HISAT2 ref-guided index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STRINGTIE_HISAT2_REF_GUIDED_ROOT"      # StringTie (ref-guided) assemblies
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_MATRIX_ROOT"         # Count matrices (ref-guided)

# --- Method 3: STAR Alignment ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_ALIGN_ROOT"                       # STAR alignment results (removes index + BAMs)
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_INDEX_ROOT"                       # STAR genome index only
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_MATRIX_ROOT"                      # Count matrices (STAR)
