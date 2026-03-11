#!/bin/bash

# ==============================================================================
# TEST CONFIG: GENOME-BASED METHODS (M1 + M3) — 3 SRRs
# ==============================================================================
# Methods tested:
#   M1 — HISAT2 Reference-Guided (genome FASTA + GTF)
#   M3 — STAR Splice-Aware Alignment (genome FASTA + GTF)
#
# Pair with HPC_test_transcript_M2_M4_M5.sh to cover all M1-M5.
# ==============================================================================

# ==============================================================================
# IMPORTANT PARAMETERS
# ==============================================================================

# Runtime Configuration
THREADS=48                              # Threads for parallel operations
JOBS=3                                  # Parallel jobs for GNU Parallel
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

	"METHOD_1_HISAT2_REF_GUIDED"
	#"METHOD_2_HISAT2_DE_NOVO"
	"METHOD_3_STAR_ALIGNMENT"
	#"METHOD_4_SALMON_SAF"
	#"METHOD_5_BOWTIE2_RSEM"
)

# ==============================================================================
# CONDA ENVIRONMENT
# ==============================================================================

eval "$(conda shell.bash hook)"
conda activate gea

# ==============================================================================
# SOURCE MODULES
# ==============================================================================

source "modules/modules_loader.sh"
#bash init_setup.sh

# Calculate threads per job
if [[ "$USE_GNU_PARALLEL" == "TRUE" ]]; then
	THREADS_PER_JOB=$((THREADS / JOBS))
	[[ $THREADS_PER_JOB -lt 1 ]] && THREADS_PER_JOB=1
else
	THREADS_PER_JOB=$THREADS
fi
export THREADS JOBS USE_GNU_PARALLEL THREADS_PER_JOB keep_bam_global

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
	"inputs/gtf/reference/GPE001970.gtf|inputs/fasta/reference_genomes/GPE001970.fa|"                              # GPE001970
	#"inputs/gtf/reference/Eggplant_V4.1_function_IPR_final.gtf|inputs/fasta/reference_genomes/Eggplant_V4.1.fa|inputs/fasta/reference_genomes/Eggplant_V4.1_transcripts.function.fa"  # Eggplant V4.1
)

# ==============================================================================
# RNA-SEQ DATA SOURCES (SRR LISTS) — 3 samples for testing
# ==============================================================================

SRR_LIST_TEST=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa
	SRR3884686	# Buds_0.7cm (flower bud initiation)
	SRR3884687	# Opened_Buds (flower development)
	SRR3884597	# Flowers (anthesis)
)

SRR_COMBINED_LIST=(
	"${SRR_LIST_TEST[@]}"
)

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

POST_PROC_ROOT="3_POST_PROC"
export POST_PROC_ROOT

# ==============================================================================
# CLEANUP OPTIONS AND TESTING ESSENTIALS
# ==============================================================================
# Uncomment lines below to remove previous results before re-running.
# WARNING: These are destructive operations — verify before uncommenting.
# ==============================================================================

ACTIVATE_RM=FALSE

# --- Method 1: HISAT2 Reference-Guided ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_ROOT"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_INDEX_DIR"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STRINGTIE_HISAT2_REF_GUIDED_ROOT"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_MATRIX_ROOT"

# --- Method 3: STAR Alignment ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_ALIGN_ROOT"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_INDEX_ROOT"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_MATRIX_ROOT"
