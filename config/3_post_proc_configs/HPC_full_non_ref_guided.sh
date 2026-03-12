#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: TRANSCRIPT-BASED METHODS (M2 + M4 + M5) — FULL DATASET
# ==============================================================================
# Purpose: Production post-processing run for transcript-alignment methods.
#          Uses the full eggplant tissue atlas and the V4.1 transcript reference.
#
# Methods:
#   M2 — HISAT2 De Novo         →  uses Eggplant_V4.1_transcripts.function
#   M4 — Salmon SAF             →  uses Eggplant_V4.1_transcripts.function
#   M5 — RSEM + Bowtie2        →  uses Eggplant_V4.1_transcripts.function
#
# Usage:
#   bash run_all_post_processing.sh --config config/3_post_proc_configs/HPC_full_non_ref_guided.sh
# ==============================================================================

# ==============================================================================
# SYSTEM RESOURCES
# ==============================================================================


THREADS=96
ENABLE_GPU="FALSE"
ENABLE_GNU_PARALLEL="TRUE"
DESIRED_CPU_PER_JOB=2
AVAILABLE_RAM_GB=64
GPU_VRAM_GB=8


# ==============================================================================
# LOGGING AND OUTPUT
# ==============================================================================

CLEAR_LOGS="FALSE"
CLEAR_OUTPUT_FOLDER="TRUE"

# ==============================================================================
# MASTER REFERENCE (uncomment ONE)
# ==============================================================================

MASTER_REFERENCES=(
    #"All_Smel_Genes"
    "Eggplant_V4.1_transcripts.function"
    #"GPE001970_genome"            # M1_HISAT2_RefGuided / M3_STAR_Align
    #"GPE001970_transcripts"       # M2_HISAT2_DeNovo
)

# ==============================================================================
# GENE GROUPS
# ==============================================================================

GENE_GROUPS=(
    "SmelDMPs_with_1_18s_rRNA"
    #"SmelDMPs"
    #"SmelDMPs_with_2_18s_rRNA"
    #"SmelDMPs_with_SmelCyclo"
    #"SmelGRF-GIF"
    #"SmelGRF-GIF_with_1_18s_rRNA"
    #"SmelGRFs"
    #"SmelGIF"
    #"Selected_GRF_GIF_Genes_vAll_GIF"
)

# ==============================================================================
# SRR DATASETS
# ==============================================================================

SRR_DATASETS=(
    "PRJNA328564"             # Main Dataset — Eggplant tissue atlas (full)
    #"PRJNA328564_selected"
    #"PRJNA865018"            # SET_1: Good Dataset for SmelDMP GEA
    #"PRJNA941250"            # SET_2: Good Dataset for SmelDMP GEA
    #"PRJNA865018_and_PRJNA941250"
    #"SAMN28540077"           # Chinese Dataset 1 — replicability
    #"SAMN28540068"           # Chinese Dataset 2 — replicability
    #"OTHER_SRR_LIST"
)

# ==============================================================================
# ALIGNMENT METHODS (comment/uncomment to enable/disable)
# ==============================================================================

METHODS=(
    #"M1_HISAT2_RefGuided"
    "M2_HISAT2_DeNovo"
    #"M3_STAR_Align"
    "M4_Salmon_Saf"
    "M5_RSEM_Bowtie2"
)

# ==============================================================================
# ANALYSES (comment/uncomment to enable/disable)
# ==============================================================================

ANALYSES=(
    #"Stringtie_Matrix"         # M2 only — skip if matrices already built

    #"Tximport_STAR"            # M3 only
    "Tximport_Salmon"           # M4 standalone (no-op if Matrix_Creation active)
    "Tximport_RSEM"             # M5 standalone (no-op if Matrix_Creation active)
    "Matrix_Creation"           # M3/M4/M5 only — NOT for HISAT2

    "Basic_Heatmap"
    "Heatmap_with_CV"
    #"BarGraph"
    #"Differential_Expression"
    #"Gene_Set_Enrichment"
    "PCA_Dimensionality_Reduction"
    "Sample_Correlation_Clustering"
    "Tissue_Specificity"
    #"Coexpression_using_WGCNA"
)
