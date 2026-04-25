#!/bin/bash

# ==============================================================================
# DOWNLOAD & TRIM CONFIG — Full SRR List
# ==============================================================================
# Purpose: Download and quality-trim all RNA-seq SRR files before alignment.
# Run this config BEFORE any alignment config (2_alignment/).
#
# Two download options (enable only one):
#   Option A — Keep raw files:  enable DOWNLOAD_SRR + TRIM_SRR
#   Option B — Auto-cleanup:    enable DOWNLOAD_TRIM_and_DELETE_RAW_SRR
# ==============================================================================

# ==============================================================================
# IMPORTANT PARAMETERS
# ==============================================================================

# Optimal Thread for TrimGalore is 4.
# Runtime Configuration
THREADS=32                              # Threads for parallel operations
JOBS=8                                  # Parallel jobs for GNU Parallel
USE_GNU_PARALLEL="TRUE"                 # TRUE/FALSE for GNU Parallel
keep_bam_global="n"                     # y=keep BAM files, n=delete after

# Pipeline Stages (comment/uncomment to enable/disable)
PIPELINE_STAGES=(
	#"MAMBA_INSTALLATION"

	# Option A: Separate download and trim (keeps raw files)
	#"DOWNLOAD_SRR"
	#"TRIM_SRR"

	# Option B: Combined download+trim+cleanup (auto-deletes raw after trim)
	"DOWNLOAD_TRIM_and_DELETE_RAW_SRR"

	"GZIP_TRIMMED_FILES"
	"QUALITY_CONTROL"

	#"DELETE_RAW_SRR"				# Manually delete raw SRR files
	#"DELETE_TRIMMED_FASTQ_FILES"	# Manually delete trimmed files

	# Alignment methods — disabled for this download-only config
	#"METHOD_1_HISAT2_REF_GUIDED"
	#"METHOD_2_HISAT2_DE_NOVO"
	#"METHOD_3_STAR_ALIGNMENT"
	#"METHOD_4_SALMON_SAF"
	#"METHOD_5_BOWTIE2_RSEM"
)

# ==============================================================================
# SHARED RUNTIME (conda, modules, thread calculation)
# ==============================================================================

source "config/shared/runtime_defaults.sh"

# ==============================================================================
# INPUT FILES AND DATA SOURCES
# ==============================================================================

# NOTE: FASTA is not used by the download/trim stages.
# A single placeholder entry is required for the pipeline loop to execute.
ALL_FASTA_FILES=(
	"inputs/fasta/reference_genome/GPE001970_transcripts.fa"
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

# Create required preprocessing directories
mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT"

# ==============================================================================
# CLEANUP OPTIONS
# ==============================================================================
# Uncomment lines below to remove previous results before re-running.
# WARNING: These are destructive operations — verify before uncommenting.
# ==============================================================================

ACTIVATE_RM=FALSE

# --- Preprocessing ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RAW_DIR_ROOT"                          # Raw SRR downloads
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$TRIM_DIR_ROOT"                         # Trimmed FASTQ files
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$FASTQC_ROOT"                           # FastQC reports
