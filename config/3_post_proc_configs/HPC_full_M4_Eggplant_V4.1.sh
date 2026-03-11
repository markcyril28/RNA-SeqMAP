#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: M4 — Salmon SAF (Eggplant_V4.1_transcripts.function)
# ==============================================================================
# Purpose: Production post-processing for M4 (Salmon quasi-mapping).
#          Uses the full eggplant tissue atlas.
#
# Method:
#   M4 — Salmon SAF  →  uses Eggplant_V4.1_transcripts.function
#
# Preprocessing:  Tximport_Salmon → Matrix_Creation (Tximport is a no-op when
#                 Matrix_Creation is active; keep both for flexibility).
# DE support:     Possible via tximport counts, but not enabled by default.
# NOT applicable: Stringtie_Matrix, Tximport_STAR, Tximport_RSEM
#
# Usage:
#   Uncomment this config in PIPELINE_CONFIGS inside run_all_post_processing.sh
# ==============================================================================

# ==============================================================================
# SYSTEM RESOURCES
# ==============================================================================

THREADS=64
ENABLE_GPU="FALSE"
ENABLE_GNU_PARALLEL="TRUE"
DESIRED_CPU_PER_JOB=4
AVAILABLE_RAM_GB=64
GPU_VRAM_GB=8

# ==============================================================================
# LOGGING AND OUTPUT
# ==============================================================================

CLEAR_LOGS="FALSE"
CLEAR_OUTPUT_FOLDER="TRUE"

# ==============================================================================
# MASTER REFERENCE
# ==============================================================================
#
# M4 uses transcript-level references.

MASTER_REFERENCES=(
    "Eggplant_V4.1_transcripts.function"
)

# ==============================================================================
# GENE GROUPS
# ==============================================================================

GENE_GROUPS=(
    "SmelDMPs_with_1_18s_rRNA"
    #"SmelDMPs"
    #"SmelDMPs_with_2_18s_rRNA"
    #"SmelDMPs_with_SmelCyclo"
    #"SmelGRF-GIFs"              # Needs stringtie_matrix_builder.sh first
    #"SmelGRF-GIFs_with_1_18s_rRNA"
    #"SmelGRFs"
    #"SmelGIFs"
    #"Selected_GRF_GIF_Genes_vAll_GIFs"
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
# ALIGNMENT METHODS
# ==============================================================================

METHODS=(
    "M4_Salmon_Saf"
)

# ==============================================================================
# ANALYSES
# ==============================================================================
#
# Tximport_Salmon is a no-op when Matrix_Creation is also active (keep for
# standalone use).  Stringtie_Matrix and Tximport_STAR/RSEM are NOT applicable.
# Tissue_Specificity is applicable to transcript-level quantification.

ANALYSES=(
    # ---- M4 preprocessing ----
    "Tximport_Salmon"           # standalone; no-op when Matrix_Creation is active
    "Matrix_Creation"           # builds count matrix from Salmon quant output

    # ---- Visualisation ----
    "Basic_Heatmap"
    "Heatmap_with_CV"
    #"BarGraph"

    # ---- Statistics ----
    #"Differential_Expression"   # Possible for M4 (tximport counts)
    "PCA_Dimensionality_Reduction"
    "Sample_Correlation_Clustering"
    "Tissue_Specificity"        # Applicable — transcript-level quantification
    #"Coexpression_using_WGCNA"
    #"Gene_Set_Enrichment"
)
