#!/bin/bash

# ==============================================================================
# IMPORTANT PARAMETERS
# ==============================================================================

# Runtime Configuration
THREADS=48                              # Threads for parallel operations
JOBS=3									# Parallel jobs for GNU Parallel
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
# Structure: logging/, a_preprocessing/, b_main_methods/, 0_input_information/
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

gtf_file="inputs/gtf/reference/GPE001970_transcripts.gtf"

# FASTA Files for Analysis
ALL_FASTA_FILES=(
	"inputs/fasta/reference_genomes/GPE001970_transcripts.fa"
)

# ==============================================================================
# RNA-SEQ DATA SOURCES (SRR LISTS)
# ==============================================================================

SRR_LIST_PRJNA328564=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa
	SRR3884686	# Buds_0.7cm (flower bud initiation) [MAIN INTEREST]
	SRR3884687	# Opened_Buds (flower development) 	 [MAIN INTEREST]
	SRR3884597	# Flowers (anthesis)				 [MAIN INTEREST]
)

SRR_LIST_SAMN28540077=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540077&o=acc_s%3Aa
	SRR20722234	# Flowers
	SRR4243802 # Buds, Adopted Dataset from ID: PRJNA341784
)

SRR_LIST_SAMN28540068=(
	#Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540068&o=acc_s%3Aa
	#SRR3884597 	# Flower — duplicate: already in SRR_LIST_PRJNA328564
	SRR20722297 # flower_buds
)

SRR_COMBINED_LIST=(
	"${SRR_LIST_PRJNA328564[@]}"	# Main Dataset for GEA.
	#"${SRR_LIST_SAMN28540077[@]}"	# Chinese Dataset for replicability.
	#"${SRR_LIST_SAMN28540068[@]}"	# Chinese Dataset for replicability.
)

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

POST_PROCESSING_ROOT="3_POST_PROC"
export POST_PROCESSING_ROOT

# Preprocessing directories are created by the main script after sourcing this config.

# ==============================================================================
# CLEANUP OPTIONS AND TESTING ESSENTIALS
# ==============================================================================
# Uncomment lines below to remove previous results before re-running.
# WARNING: These are destructive operations — verify before uncommenting.
# ==============================================================================

ACTIVATE_RM=FALSE

# --- Method 2: HISAT2 De Novo ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_ROOT"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_INDEX_DIR"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STRINGTIE_HISAT2_DE_NOVO_ROOT"
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_MATRIX_ROOT"
